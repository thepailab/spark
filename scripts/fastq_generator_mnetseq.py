import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import string
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
                frag_start = int(row.fragment_start_pos)
                frag_end = int(row.fragment_end_pos)
                
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                frag_end = min(frag_end, len(current_ref))
                frag_start = min(frag_start, frag_end)

                if strandedness == "rf":
                    d_start = max(0, frag_end - read_length)
                    seq = "".join(current_ref[d_start:frag_end])
                    rev = True
                    abs_coords = get_absolute_coords(rates_df, d_start, frag_end - 1)
                elif strandedness == "fr":
                    u_end = min(frag_start + read_length, frag_end)
                    seq = "".join(current_ref[frag_start:u_end])
                    rev = False
                    abs_coords = get_absolute_coords(rates_df, frag_start, u_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        d_start = max(0, frag_end - read_length)
                        seq = "".join(current_ref[d_start:frag_end])
                        rev = True
                        abs_coords = get_absolute_coords(rates_df, d_start, frag_end - 1)
                    else:
                        u_end = min(frag_start + read_length, frag_end)
                        seq = "".join(current_ref[frag_start:u_end])
                        rev = False
                        abs_coords = get_absolute_coords(rates_df, frag_start, u_end - 1)

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
                frag_start = int(row.fragment_start_pos)
                frag_end = int(row.fragment_end_pos)
                
                is_bg = "_BG_" in str(row.molecule_id)
                current_ref = bg_ref_arr if is_bg else ref_arr
                
                frag_end = min(frag_end, len(current_ref))
                frag_start = min(frag_start, frag_end)
                
                u_end = min(frag_start + read_length, frag_end)
                d_start = max(0, frag_end - read_length)
                
                seq_upstream = "".join(current_ref[frag_start:u_end])
                seq_downstream = "".join(current_ref[d_start:frag_end])
                
                if strandedness == "rf":
                    seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                    seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                    abs_coords_r1 = get_absolute_coords(rates_df, d_start, frag_end - 1)
                    abs_coords_r2 = get_absolute_coords(rates_df, frag_start, u_end - 1)
                elif strandedness == "fr":
                    seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                    seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                    abs_coords_r1 = get_absolute_coords(rates_df, frag_start, u_end - 1)
                    abs_coords_r2 = get_absolute_coords(rates_df, d_start, frag_end - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seq_r1 = process_sequence(seq_downstream, reverse_comp=True)
                        seq_r2 = process_sequence(seq_upstream, reverse_comp=False)
                        abs_coords_r1 = get_absolute_coords(rates_df, d_start, frag_end - 1)
                        abs_coords_r2 = get_absolute_coords(rates_df, frag_start, u_end - 1)
                    else:
                        seq_r1 = process_sequence(seq_upstream, reverse_comp=False)
                        seq_r2 = process_sequence(seq_downstream, reverse_comp=True)
                        abs_coords_r1 = get_absolute_coords(rates_df, frag_start, u_end - 1)
                        abs_coords_r2 = get_absolute_coords(rates_df, d_start, frag_end - 1)

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
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split("_")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/mnetseqmRNAs/{base_filename}_mnetseq.tsv.gz"

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
        
        gene_tpm = np.random.uniform(low=int(args.tpm_lower_limit), high=int(args.tpm_upper_limit))
        seq_depth = args.seq_depth / 1e6
        reads_to_get = gene_length * gene_tpm * seq_depth

        if len(df) >= reads_to_get:
            df = df.sample(n=int(reads_to_get), replace=False, random_state=None)
        else:
            sampled_reads = df.sample(n=int(reads_to_get) - len(df), replace=True, random_state=None)
            df = pd.concat([df, sampled_reads], ignore_index=True)

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
            
    except FileNotFoundError:
        pass
