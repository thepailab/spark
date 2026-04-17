import pandas as pd
import random
import string
import argparse
import gzip
import os
import numpy as np
import ast

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

def reconstruct_main(ref_arr, stop_pos, conv_pos_str, err_pos_str, spliced_introns_str, sub_base):
    mol_arr = ref_arr[:int(stop_pos)].copy()
    
    if pd.notna(conv_pos_str) and conv_pos_str not in ("[]", "NA", "") and sub_base:
        try:
            conv_pos = ast.literal_eval(conv_pos_str)
            valid = [p for p in conv_pos if p < stop_pos]
            if valid:
                mol_arr[valid] = sub_base
        except:
            pass
            
    if pd.notna(err_pos_str) and err_pos_str not in ("[]", "NA", ""):
        try:
            err_pos = ast.literal_eval(err_pos_str)
            valid = [p for p in err_pos if p < stop_pos]
            for p in valid:
                orig = mol_arr[p]
                poss = [b for b in ['A','C','G','T'] if b != orig]
                mol_arr[p] = random.choice(poss)
        except:
            pass
            
    seq_str = "".join(mol_arr)
    
    if pd.notna(spliced_introns_str) and spliced_introns_str not in ("[]", "NA", ""):
        try:
            spliced_introns = ast.literal_eval(spliced_introns_str)
            spliced_introns.sort(key=lambda x: x[0], reverse=True)
            for start, end in spliced_introns:
                seq_str = seq_str[:start] + seq_str[end:]
        except:
            pass
            
    return seq_str

def reconstruct_bg(bg_ref_arr, err_pos_str):
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
    return "".join(mol_arr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")    
    parser.add_argument("--seq_type", choices=["RNA", "cDNA"], default='RNA')
    parser.add_argument('--seq_depth', type=int, default=20000000)
    parser.add_argument('--tpm_lower_limit', type=int, default=5)
    parser.add_argument('--tpm_upper_limit', type=int, default=200)
    parser.add_argument('-o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--bkg_molecules", type=float, default=0.0)
    parser.add_argument("--sub_base", type=str, default='C')

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, 'reads', f"{base_filename}_{args.seq_type}.fastq.gz")
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    gtf_path = os.path.join(args.o, "gtf", f"{base_filename}.tsv.gz")
    if os.path.exists(gtf_path):
        df_gtf = pd.read_csv(gtf_path, sep="\t")
        df_gtf = df_gtf.sort_values('position')
        
        ref_seq_full = "".join(df_gtf['sequence'].astype(str).tolist())
        ref_arr = np.array(list(ref_seq_full), dtype='U1')
        
        is_bg_spliced = True
        path_to_BGmRNAs = os.path.join(args.o, "mRNA", f"{base_filename}_background.tsv.gz")
        
        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            try:
                peek_df = pd.read_csv(path_to_BGmRNAs, delimiter="\t", nrows=1)
                if 'sequence_length' in peek_df.columns and peek_df.iloc[0]['sequence_length'] == len(ref_seq_full):
                    is_bg_spliced = False
            except:
                pass
                
        if is_bg_spliced:
            exon_gtf = df_gtf[df_gtf['feature'] == 'exon']
            bg_ref_seq_spliced = "".join(exon_gtf['sequence'].astype(str).tolist())
            bg_ref_arr = np.array(list(bg_ref_seq_spliced), dtype='U1')
        else:
            bg_ref_arr = ref_arr.copy()
            
    else:
        print(f"[Warning] GTF not found at {gtf_path}. Cannot assemble sequences.")
        exit(1)

    cpm_gene = random.randint(args.tpm_lower_limit, args.tpm_upper_limit)
    reads_to_get_gene = int(round(args.seq_depth * cpm_gene / 10**6))
    
    reads_to_get_bg = 0
    if args.bkg_molecules > 0.0:
        if args.bkg_molecules >= 1.0:
            reads_to_get_bg = 1000 
            reads_to_get_gene = 0
        elif reads_to_get_gene > 0:
            ratio = args.bkg_molecules / (1.0 - args.bkg_molecules)
            reads_to_get_bg = int(reads_to_get_gene * ratio)

    print(f"[LongRead] Target: {reads_to_get_gene} gene reads, {reads_to_get_bg} background reads.")

    with gzip.open(output_filename, 'wt') as f_out:
        
        if reads_to_get_gene > 0:
            df = pd.read_csv(args.input_df, delimiter="\t")
            
            if not df.empty:
                if len(df) >= reads_to_get_gene:
                    sampled_df = df.sample(n=reads_to_get_gene, replace=False)
                else:
                    extra_needed = reads_to_get_gene - len(df)
                    additional_samples = df.sample(n=extra_needed, replace=True)
                    sampled_df = pd.concat([df, additional_samples], ignore_index=True)
                
                for row in sampled_df.itertuples():
                    molecule_id = id_generator()
                    read_name = f"{molecule_id}_{base_filename}"
                    
                    raw_seq = reconstruct_main(
                        ref_arr, 
                        getattr(row, 'stop_label_pos', 0), 
                        getattr(row, 'converted_positions', '[]'), 
                        getattr(row, 'seq_err_positions', '[]'), 
                        getattr(row, 'spliced_introns', '[]'), 
                        args.sub_base
                    )
                    
                    sequence = process_sequence(raw_seq, args.seq_type)
                    quality_scores = "I" * len(sequence)

                    f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")
            del df 

        if reads_to_get_bg > 0:
            path_to_BGmRNAs = os.path.join(args.o, "mRNA", f"{base_filename}_background.tsv.gz")
            
            if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
                try:
                    df_bg = pd.read_csv(path_to_BGmRNAs, delimiter="\t")
                    
                    if not df_bg.empty:
                        bg_arr = bg_ref_arr.copy()

                        if len(df_bg) >= reads_to_get_bg:
                            sampled_bg = df_bg.sample(n=reads_to_get_bg, replace=False)
                        else:
                            extra_needed = reads_to_get_bg - len(df_bg)
                            additional_samples = df_bg.sample(n=extra_needed, replace=True)
                            sampled_bg = pd.concat([df_bg, additional_samples], ignore_index=True)

                        for row in sampled_bg.itertuples():
                            molecule_id = id_generator()
                            read_name = f"{molecule_id}_{base_filename}_BG"
                            
                            raw_seq = reconstruct_bg(bg_arr, getattr(row, 'seq_err_positions', '[]'))
                            
                            sequence = process_sequence(raw_seq, args.seq_type)
                            quality_scores = "I" * len(sequence)

                            f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")
                    del df_bg
                except Exception as e:
                    print(f"Warning: Failed to process background file: {e}")
