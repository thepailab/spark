import pandas as pd
import random
import string
import argparse
import gzip
import os
import numpy as np
import ast
import hashlib
from spark_reads import pool_paths, write_units, read_pool, remove_paths, sample_rows, gene_references
from spark_labeling import parse_conversion_map

def reverse_complement(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def process_sequence(sequence, sequencing_type):
    if sequencing_type == "RNA":
        return sequence[::-1]
    elif sequencing_type == "cDNA":
        return sequence if random.choice([True, False]) else reverse_complement(sequence)
    return sequence

def id_generator(size=12, chars=string.ascii_uppercase + string.digits):
    return "@" + ''.join(random.choice(chars) for _ in range(size - 1))   

def reconstruct_main(ref_arr, stop_pos, conv_pos_str, err_pos_str, spliced_introns_str, conv_map, polyA_len):
    mol_arr = ref_arr[:int(stop_pos)].copy()

    if pd.notna(conv_pos_str) and conv_pos_str not in ("[]", "NA", "") and conv_map:
        try:
            conv_pos = ast.literal_eval(conv_pos_str)
            valid = [p for p in conv_pos if p < len(mol_arr)]
            if valid:
                orig_bases = mol_arr[valid]
                mol_arr[valid] = [conv_map.get(b, b) for b in orig_bases]
        except: pass
            
    if pd.notna(err_pos_str) and err_pos_str not in ("[]", "NA", ""):
        try:
            err_pos = ast.literal_eval(err_pos_str)
            valid = [p for p in err_pos if p < len(mol_arr)]
            for p in valid:
                orig = mol_arr[p]
                poss = [b for b in ['A','C','G','T'] if b != orig]
                mol_arr[p] = random.choice(poss)
        except: pass
            
    keep_mask = np.ones(len(mol_arr), dtype=bool)
    
    keep_mask[int(stop_pos):] = False
    
    if pd.notna(spliced_introns_str) and spliced_introns_str not in ("[]", "NA", ""):
        try:
            spliced_introns = ast.literal_eval(spliced_introns_str)
            for start, end in spliced_introns:
                keep_mask[start:end] = False
        except: pass
            
    seq_str = "".join(mol_arr[keep_mask])
    
    if int(polyA_len) > 0:
        seq_str += "A" * int(polyA_len)
            
    return seq_str

def reconstruct_bg(bg_ref_arr, err_pos_str, polyA_len):
    mol_arr = bg_ref_arr.copy()
    if pd.notna(err_pos_str) and err_pos_str not in ("[]", "NA", ""):
        try:
            err_pos = ast.literal_eval(err_pos_str)
            valid = [p for p in err_pos if p < len(mol_arr)]
            for p in valid:
                orig = mol_arr[p]
                poss = [b for b in ['A','C','G','T'] if b != orig]
                mol_arr[p] = random.choice(poss)
        except:
            pass
            
    seq_str = "".join(mol_arr)
    if int(polyA_len) > 0:
        seq_str += "A" * int(polyA_len)
        
    return seq_str

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")   
    parser.add_argument("--seq_type", choices=["RNA", "cDNA"], default='RNA')
    parser.add_argument('-o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--bkg_molecules", type=float, default=0.0)
    parser.add_argument("--labeling_base", type=str, default='T')
    parser.add_argument("--sub_base", type=str, default='C')
    parser.add_argument("--experiment_type", type=str, default="metaboliclabeling")
    parser.add_argument("--prepare", action="store_true", help="count the molecule pools instead of writing reads")
    parser.add_argument("--nascent_reads", type=int, default=0)
    parser.add_argument("--bkg_reads", type=int, default=0)

    args = parser.parse_args()

    _from_bases, conv_map = parse_conversion_map(args.labeling_base, args.sub_base)

    filename = os.path.basename(args.input_df)
    base_filename = filename.replace(".tsv.gz", "")

    
    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % 4294967295 
        random.seed(unique_seed)
        np.random.seed(unique_seed)
    

    output_filename = os.path.join(args.o, 'reads', f"{base_filename}_{args.seq_type}.fastq.gz")
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    path_to_BGmRNAs = os.path.join(args.o.rstrip('/'), "mRNA", base_filename + "_background.tsv.gz")
    _, _, units_path = pool_paths(args.o, base_filename)

    gtf_path = os.path.join(args.o, "gtf", f"{base_filename}.tsv.gz")
    if not os.path.exists(gtf_path):
        print(f"[Warning] GTF not found at {gtf_path}. Cannot assemble sequences.")
        exit(1)

    df = read_pool(args.input_df)
    if not df.empty and args.experiment_type != 'chromatin_associated':
        df = df[df['incorporated_positions'].notna() & (df['incorporated_positions'] != '[]')]
    df_bg = read_pool(path_to_BGmRNAs) if args.bkg_molecules > 0.0 else pd.DataFrame()

    if args.prepare:
        write_units(units_path, base_filename, len(df), len(df_bg), len(df_bg))
        exit(0)

    ref_arr, bg_ref_arr = gene_references(gtf_path, path_to_BGmRNAs)

    print(f"[LongRead] Target: {args.nascent_reads} gene reads, {args.bkg_reads} background reads.")

    with gzip.open(output_filename, 'wt') as f_out:
        for row in sample_rows(df, args.nascent_reads).itertuples():
            molecule_id = str(getattr(row, 'molecule_id', id_generator()))
            read_name = f"{molecule_id}_{base_filename}"
            
            raw_seq = reconstruct_main(
                ref_arr,
                getattr(row, 'stop_label_pos', 0), 
                getattr(row, 'converted_positions', '[]'), 
                getattr(row, 'seq_err_positions', '[]'), 
                getattr(row, 'spliced_introns', '[]'),
                conv_map,
                getattr(row, 'polyA_length', 0)
            )
            
            sequence = process_sequence(raw_seq, args.seq_type)
            quality_scores = "I" * len(sequence)

            f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")

        for row in sample_rows(df_bg, args.bkg_reads).itertuples():
            molecule_id = str(getattr(row, 'molecule_id', id_generator()))
            read_name = f"{molecule_id}_{base_filename}_BG"
            
            raw_seq = reconstruct_bg(bg_ref_arr.copy(), getattr(row, 'seq_err_positions', '[]'), getattr(row, 'polyA_length', 0))
            
            sequence = process_sequence(raw_seq, args.seq_type)
            quality_scores = "I" * len(sequence)

            f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")

    remove_paths([units_path])
