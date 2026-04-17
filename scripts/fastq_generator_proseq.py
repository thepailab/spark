import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import string
import ast
import numpy as np

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

def get_absolute_coords(rates_df, read_start, read_end):
    if rates_df.empty:
        return "NA"
    subset = rates_df[
        (rates_df['nucleotide_coord'] >= read_start) &
        (rates_df['nucleotide_coord'] <= read_end)
    ]
    if not subset.empty:
        chrom = subset['chromosome'].iloc[0]
        abs_start = subset['absolute_position'].iloc[0]
        abs_end = subset['absolute_position'].iloc[-1]
        strand = subset['strand'].iloc[0]
        return f"{chrom}:{abs_start}-{abs_end}:{strand}"
    return "NA"

def process_and_write_fastq(df, output_prefix, sequencing_type, strandedness, read_length, rates_df, ref_arr, bg_ref_arr):
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
                    abs_coords = get_absolute_coords(rates_df, d_start, r_end - 1)
                elif strandedness == "fr":
                    seq = "".join(current_ref[r_start:u_end])
                    rev = False
                    abs_coords = get_absolute_coords(rates_df, r_start, u_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq = "".join(current_ref[d_start:r_end])
                        rev = True
                        abs_coords = get_absolute_coords(rates_df, d_start, r_end - 1)
                    else:
                        seq = "".join(current_ref[r_start:u_end])
                        rev = False
                        abs_coords = get_absolute_coords(rates_df, r_start, u_end - 1)

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
                    abs_coords_r1 = get_absolute_coords(rates_df, d_start, r_end - 1)
                    abs_coords_r2 = get_absolute_coords(rates_df, r_start, u_end - 1)
                elif strandedness == "fr":
                    seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                    seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                    abs_coords_r1 = get_absolute_coords(rates_df, r_start, u_end - 1)
                    abs_coords_r2 = get_absolute_coords(rates_df, d_start, r_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                        seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                        abs_coords_r1 = get_absolute_coords(rates_df, d_start, r_end - 1)
                        abs_coords_r2 = get_absolute_coords(rates_df, r_start, u_end - 1)
                    else:
                        seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                        seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                        abs_coords_r1 = get_absolute_coords(rates_df, r_start, u_end - 1)
                        abs_coords_r2 = get_absolute_coords(rates_df, d_start, r_end - 1)

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
    parser.add_argument('--seq_depth', type=int, default=20000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--tpm_lower_limit', type=int, default=5)
    parser.add_argument('--tpm_upper_limit', type=int, default=200)
    parser.add_argument('--bkg_molecules', type=float, default=0)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    cmd_chop = [
        f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper_proseq.R')}",
        f"--tsv {args.input_df}",
        f"--insert_size {args.insert_size}",
        f"--read_length {args.read_length}",
        f"--threads {args.threads}",
        f"--seq_depth {args.seq_depth}",
        f"--tpm_lower_limit {args.tpm_lower_limit}",
        f"--tpm_upper_limit {args.tpm_upper_limit}",
        f"-o {args.o}"
    ]
    if args.seed:
        cmd_chop += ['--seed', str(args.seed)]
    run_cmd(" ".join(cmd_chop))

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split("_")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"

    gtf_path = f"{args.o}/gtf/{base_filename}.tsv.gz"
    if os.path.exists(gtf_path):
        df_gtf = pd.read_csv(gtf_path, delimiter="\t")
        df_gtf = df_gtf.sort_values('position')
        
        ref_seq_full = "".join(df_gtf['sequence'].astype(str).tolist())
        ref_arr = np.array(list(ref_seq_full), dtype='U1')
        
        is_bg_spliced = True
        path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            try:
                peek_df = pd.read_csv(path_to_BGmRNAs, delimiter="\t", nrows=1)
                if 'sequence_length' in peek_df.columns and peek_df.iloc[0]['sequence_length'] == len(ref_seq_full):
                    is_bg_spliced = False
            except:
                pass
        
        if is_bg_spliced:
            exon_df = df_gtf[df_gtf['feature'] == 'exon'].copy()
            bg_ref_seq_spliced = "".join(exon_df['sequence'].astype(str).tolist())
            bg_ref_arr = np.array(list(bg_ref_seq_spliced), dtype='U1')
            gene_length = exon_df["sequence"].str.len().sum() / 1000
        else:
            bg_ref_arr = ref_arr.copy()
            gene_length = len(ref_seq_full) / 1000
    else:
        exit(1)

    try:
        df = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
        result_df = df.copy()

        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            cmd_chop_bg = [
                f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper_TTseq.R')}",
                f"--tsv {path_to_BGmRNAs}",
                f"--insert_size {args.insert_size}",
                f"--read_length {args.read_length}",
                f"--threads {args.threads}",
                f"-o {args.o}"
            ]
            if args.seed:
                cmd_chop_bg += ['--seed', str(args.seed)]

            run_cmd(" ".join(cmd_chop_bg))
            chopped_coordinates_BG_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
            
            try:
                df_bg = pd.read_csv(chopped_coordinates_BG_file_path, delimiter="\t")
                df_bg['molecule_id'] = df_bg['molecule_id'].astype(str) + '_BG_'
                
                if 'read_coordinates' in df_bg.columns:
                    df_bg['read_coordinates'] = df_bg['read_coordinates'].astype(str).str.split(',')
                    df_bg = df_bg.explode('read_coordinates').reset_index(drop=True)
                    coord_split = df_bg['read_coordinates'].str.split('-', expand=True)
                    df_bg['read_start'] = pd.to_numeric(coord_split[0], errors='coerce')
                    df_bg['read_end'] = pd.to_numeric(coord_split[1], errors='coerce')
                    df_bg = df_bg.dropna(subset=['read_start', 'read_end'])

                if not df_bg.empty and not result_df.empty:
                    n_bg = int(len(result_df) * args.bkg_molecules)
                    if n_bg > 0:
                        if len(df_bg) >= n_bg:
                            df_bg_sampled = df_bg.sample(n=n_bg, replace=False, random_state=42)
                        else:
                            df_bg_sampled = df_bg.sample(n=n_bg, replace=True, random_state=42)
                    
                        n_main = int(len(result_df) * (1 - args.bkg_molecules))
                        result_df = result_df.sample(n=n_main, replace=False, random_state=42)
                        result_df = pd.concat([result_df, df_bg_sampled]).reset_index(drop=True)
            except FileNotFoundError:
                pass

        gene_tpm = np.random.uniform(low=int(args.tpm_lower_limit), high=int(args.tpm_upper_limit))
        seq_depth = args.seq_depth / 1e6
        reads_to_get = int(gene_length * gene_tpm * seq_depth)

        if len(result_df) >= reads_to_get:
            result_df = result_df.sample(n=reads_to_get, replace=False, random_state=None)
        elif not result_df.empty:
            sampled_reads = result_df.sample(n=reads_to_get - len(result_df), replace=True, random_state=None)
            result_df = pd.concat([result_df, sampled_reads], ignore_index=True)

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

        try:
            if os.path.exists(chopped_coordinates_file_path):
                os.remove(chopped_coordinates_file_path)
            if 'chopped_coordinates_BG_file_path' in locals() and os.path.exists(chopped_coordinates_BG_file_path):
                os.remove(chopped_coordinates_BG_file_path)
        except OSError:
            pass

    except FileNotFoundError:
        pass
