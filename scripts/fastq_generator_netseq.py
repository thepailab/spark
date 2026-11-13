import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import string
import numpy as np
import hashlib
from spark_reads import pool_paths, write_units, read_pool, remove_paths, sample_rows, gene_references, coordinate_mapper

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=True)

def make_random_suffix(length=5):
    return ''.join(random.choices(string.ascii_lowercase, k=length))

def reverse_complement(seq):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def process_sequence(sequence, reverse_comp=False):
    return reverse_complement(sequence) if reverse_comp else sequence

def process_and_write_fastq(df, output_prefix, sequencing_type, strandedness, read_length, rates_df, ref_arr, bg_ref_arr):
    to_absolute = coordinate_mapper(rates_df)
    if sequencing_type == "SE":
        output_file = f"{output_prefix}.fastq.gz"
        with gzip.open(output_file, 'wt') as f:
            for row in df.itertuples():
                start_offset = int(getattr(row, 'start_offset', 0))
                frag_start = int(row.fragment_start_pos)
                frag_end = int(row.fragment_end_pos)
                
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                frag_end = min(frag_end, len(current_ref))
                frag_start = min(frag_start, frag_end)

                if strandedness == "rf":
                    
                    d_start = max(frag_start, frag_end - read_length)
                    seq = "".join(current_ref[d_start:frag_end])
                    rev = True
                    abs_coords = to_absolute(d_start, frag_end - 1)
                elif strandedness == "fr":
                    u_end = min(frag_start + read_length, frag_end)
                    seq = "".join(current_ref[frag_start:u_end])
                    rev = False
                    abs_coords = to_absolute(frag_start, u_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        d_start = max(frag_start, frag_end - read_length)
                        seq = "".join(current_ref[d_start:frag_end])
                        rev = True
                        abs_coords = to_absolute(d_start, frag_end - 1)
                    else:
                        u_end = min(frag_start + read_length, frag_end)
                        seq = "".join(current_ref[frag_start:u_end])
                        rev = False
                        abs_coords = to_absolute(frag_start, u_end - 1)

                seq = process_sequence(seq, reverse_comp=rev)
                if not seq: continue

                random_suffix = make_random_suffix()
                read_name = f"{row.molecule_id}{random_suffix}_{abs_coords}"
                quality_scores = "I" * len(seq)
                
                f.write(f"{read_name}\n{seq}\n+\n{quality_scores}\n")

    elif sequencing_type == "PE":
        output_file_r1 = f"{output_prefix}_R1.fastq.gz"
        output_file_r2 = f"{output_prefix}_R2.fastq.gz"
        
        with gzip.open(output_file_r1, 'wt') as f1, gzip.open(output_file_r2, 'wt') as f2:
            for row in df.itertuples():
                start_offset = int(getattr(row, 'start_offset', 0))
                frag_start = int(row.fragment_start_pos)
                frag_end = int(row.fragment_end_pos)
                
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                frag_end = min(frag_end, len(current_ref))
                frag_start = min(frag_start, frag_end)
                
                u_end = min(frag_start + read_length, frag_end)
                
                
                d_start = max(frag_start, frag_end - read_length)
                
                seq_upstream = "".join(current_ref[frag_start:u_end])
                seq_downstream = "".join(current_ref[d_start:frag_end])
                
                if strandedness == "rf":
                    seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                    seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                    abs_coords_r1 = to_absolute(d_start, frag_end - 1)
                    abs_coords_r2 = to_absolute(frag_start, u_end - 1)
                elif strandedness == "fr":
                    seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                    seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                    abs_coords_r1 = to_absolute(frag_start, u_end - 1)
                    abs_coords_r2 = to_absolute(d_start, frag_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                        seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                        abs_coords_r1 = to_absolute(d_start, frag_end - 1)
                        abs_coords_r2 = to_absolute(frag_start, u_end - 1)
                    else:
                        seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                        seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                        abs_coords_r1 = to_absolute(frag_start, u_end - 1)
                        abs_coords_r2 = to_absolute(d_start, frag_end - 1)

                if not seq_r1 or not seq_r2: continue

                random_suffix = make_random_suffix()
                read_name_r1 = f"{row.molecule_id}{random_suffix}_{abs_coords_r1}"
                read_name_r2 = f"{row.molecule_id}{random_suffix}_{abs_coords_r2}"
                
                f1.write(f"{read_name_r1}\n{seq_r1}\n+\n{'I'*len(seq_r1)}\n")
                f2.write(f"{read_name_r2}\n{seq_r2}\n+\n{'I'*len(seq_r2)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")
    parser.add_argument("--insert_size", type=str, default="200,300")
    parser.add_argument("--read_length", type=int, default=100)
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE')
    parser.add_argument("--s", choices=["rf", "fr", "unstranded"], type=str, default='rf')
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--prepare", action="store_true", help="count the footprint pools instead of writing reads")
    parser.add_argument("--nascent_reads", type=int, default=0)
    parser.add_argument("--bkg_reads", type=int, default=0)

    args = parser.parse_args()

    filename = os.path.basename(args.input_df)
    base_filename = filename.replace("_netseq.tsv.gz", "")

    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % 4294967295 
        random.seed(unique_seed)
        np.random.seed(unique_seed)

    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/netseqmRNAs/{base_filename}_netseq.tsv.gz"
    path_to_BGmRNAs = os.path.join(args.o.rstrip('/'), "mRNA", base_filename + "_background.tsv.gz")
    _, _, units_path = pool_paths(args.o, base_filename)

    gtf_path = f"{args.o}/gtf/{base_filename}.tsv.gz"
    if not os.path.exists(gtf_path):
        exit(1)

    df = read_pool(chopped_coordinates_file_path)
    is_bg = df['molecule_id'].astype(str).str.contains('_BG_') if not df.empty else pd.Series(dtype=bool)
    df_main = df[~is_bg] if not df.empty else df
    df_bg = df[is_bg] if not df.empty else df

    if args.prepare:
        write_units(units_path, base_filename, len(df_main), len(df_bg), len(df_bg))
        exit(0)

    df = pd.concat([sample_rows(df_main, args.nascent_reads), sample_rows(df_bg, args.bkg_reads)], ignore_index=True)

    if not df.empty:
        ref_arr, bg_ref_arr = gene_references(gtf_path, path_to_BGmRNAs)

        rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
        if os.path.exists(rates_for_gene):
            rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
        else:
            rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])

        process_and_write_fastq(
            df, 
            output_prefix, 
            args.seq_type, 
            args.s, 
            args.read_length, 
            rates_df, 
            ref_arr, 
            bg_ref_arr
        )

    remove_paths([units_path])
