import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import string
import ast
import numpy as np
import hashlib
import sys
from spark_reads import pool_paths, write_units, read_pool, remove_paths, sample_rows, count_rows, gene_references, coordinate_mapper

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
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                try:
                    r_start = int(row.read_start) - 1
                    r_end = int(row.read_end)
                except:
                    continue

                r_start = max(0, r_start)
                r_end = min(r_end, len(current_ref))
                if r_start >= r_end: continue

                d_start = max(r_start, r_end - read_length)
                u_end = min(r_start + read_length, r_end)

                if strandedness == "rf":
                    seq = "".join(current_ref[d_start:r_end])
                    rev = True
                    abs_coords = to_absolute(d_start, r_end - 1)
                elif strandedness == "fr":
                    seq = "".join(current_ref[r_start:u_end])
                    rev = False
                    abs_coords = to_absolute(r_start, u_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq = "".join(current_ref[d_start:r_end])
                        rev = True
                        abs_coords = to_absolute(d_start, r_end - 1)
                    else:
                        seq = "".join(current_ref[r_start:u_end])
                        rev = False
                        abs_coords = to_absolute(r_start, u_end - 1)

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
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                try:
                    r_start = int(row.read_start) - 1
                    r_end = int(row.read_end)
                except:
                    continue

                r_start = max(0, r_start)
                r_end = min(r_end, len(current_ref))
                if r_start >= r_end: continue

                d_start = max(r_start, r_end - read_length)
                u_end = min(r_start + read_length, r_end)

                seq_upstream = "".join(current_ref[r_start:u_end])
                seq_downstream = "".join(current_ref[d_start:r_end])

                if strandedness == "rf":
                    seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                    seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                    abs_coords_r1 = to_absolute(d_start, r_end - 1)
                    abs_coords_r2 = to_absolute(r_start, u_end - 1)
                elif strandedness == "fr":
                    seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                    seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                    abs_coords_r1 = to_absolute(r_start, u_end - 1)
                    abs_coords_r2 = to_absolute(d_start, r_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                        seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                        abs_coords_r1 = to_absolute(d_start, r_end - 1)
                        abs_coords_r2 = to_absolute(r_start, u_end - 1)
                    else:
                        seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                        seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                        abs_coords_r1 = to_absolute(r_start, u_end - 1)
                        abs_coords_r2 = to_absolute(d_start, r_end - 1)

                if not seq_r1 or not seq_r2: continue

                random_suffix = make_random_suffix()
                read_name_r1 = f"{row.molecule_id}{random_suffix}_{abs_coords_r1}"
                read_name_r2 = f"{row.molecule_id}{random_suffix}_{abs_coords_r2}"
                
                f1.write(f"{read_name_r1}\n{seq_r1}\n+\n{'I'*len(seq_r1)}\n")
                f2.write(f"{read_name_r2}\n{seq_r2}\n+\n{'I'*len(seq_r2)}\n")

def run_chopper(tsv_path, args):
    cmd_chop = [
        f"python {os.path.join(SCRIPT_DIR, 'short_read_chopper.py')}",
        "--mode proseq",
        f"--tsv {tsv_path}",
        f"--insert_size {args.insert_size}",
        f"--read_length {args.read_length}",
        f"--threads {args.threads}",
        f"-o {args.o}"
    ]
    if args.seed is not None:
        cmd_chop += ['--seed', str(args.seed)]
    run_cmd(" ".join(cmd_chop))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")
    parser.add_argument("--insert_size", type=str, default="200,300")
    parser.add_argument("--read_length", type=int, default=100)
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE')
    parser.add_argument("--s", choices=["rf", "fr", "unstranded"], type=str, default='rf')
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--bkg_molecules', type=float, default=0)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--prepare", action="store_true", help="build the fragment pools and count them instead of writing reads")
    parser.add_argument("--nascent_reads", type=int, default=0)
    parser.add_argument("--bkg_reads", type=int, default=0)

    args = parser.parse_args()

    filename = os.path.basename(args.input_df)
    base_filename = filename.replace("_proseq.tsv.gz", "")

    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % 4294967295
        random.seed(unique_seed)
        np.random.seed(unique_seed)

    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"
    chopped_coordinates_BG_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
    path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
    nascent_pool_path, bg_pool_path, units_path = pool_paths(args.o, base_filename)

    gtf_path = f"{args.o}/gtf/{base_filename}.tsv.gz"
    if not os.path.exists(gtf_path):
        exit(1)

    if args.prepare:
        run_chopper(args.input_df, args)
        nascent_pool = read_pool(chopped_coordinates_file_path)

        bg_pool = pd.DataFrame()
        bg_molecules = 0
        if args.bkg_molecules > 0.0 and count_rows(path_to_BGmRNAs) > 0:
            run_chopper(path_to_BGmRNAs, args)
            bg_pool = read_pool(chopped_coordinates_BG_file_path)
            if not bg_pool.empty:
                bg_pool['molecule_id'] = bg_pool['molecule_id'].astype(str) + '_BG_'
            bg_molecules = count_rows(path_to_BGmRNAs)

        nascent_pool.to_csv(nascent_pool_path, sep="\t", index=False)
        bg_pool.to_csv(bg_pool_path, sep="\t", index=False)
        write_units(units_path, base_filename, len(nascent_pool), len(bg_pool), bg_molecules)
        remove_paths([chopped_coordinates_file_path, chopped_coordinates_BG_file_path])
        sys.exit(0)

    result_df = pd.concat([
        sample_rows(read_pool(nascent_pool_path), args.nascent_reads),
        sample_rows(read_pool(bg_pool_path), args.bkg_reads),
    ], ignore_index=True)

    if not result_df.empty:
        ref_arr, bg_ref_arr = gene_references(gtf_path, path_to_BGmRNAs)

        rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
        if os.path.exists(rates_for_gene):
            rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
        else:
            rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])
            
        process_and_write_fastq(
            result_df, 
            output_prefix, 
            args.seq_type, 
            args.s,
            args.read_length, 
            rates_df, 
            ref_arr, 
            bg_ref_arr
        )

    remove_paths([nascent_pool_path, bg_pool_path, units_path])
