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
    if rates_df.empty: return "NA"
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
    else:
        return "NA"

def process_and_write_fastq(
    df_sampled, lookup_dict, output_prefix, sequencing_type, 
    strandedness, read_length, rates_df, ref_arr, bg_ref_arr, sub_base='C', ttseq=False
):
    if sequencing_type == "SE":
        f_out = gzip.open(f"{output_prefix}.fastq.gz", 'wt')
    else:
        f1_out = gzip.open(f"{output_prefix}_R1.fastq.gz", 'wt')
        f2_out = gzip.open(f"{output_prefix}_R2.fastq.gz", 'wt')

    try:
        for row in df_sampled.itertuples():
            ref_data = lookup_dict.get(row.transcript_id)
            if not ref_data: continue 
            
            is_bg = "_BG_" in ref_data['mol_id']
            mol_arr = bg_ref_arr.copy() if is_bg else ref_arr.copy()
            orig_idx_arr = np.arange(len(mol_arr))
            
            is_inc = np.zeros(len(mol_arr), dtype=bool)
            is_conv = np.zeros(len(mol_arr), dtype=bool)
            
            if not is_bg:
                inc_pos = [p for p in ref_data.get('inc', []) if p < len(mol_arr)]
                conv_pos = [p for p in ref_data.get('conv', []) if p < len(mol_arr)]
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
                
                is_inc[inc_pos] = True
                is_conv[conv_pos] = True
                if sub_base: mol_arr[conv_pos] = sub_base
                
                for p in err_pos:
                    mol_arr[p] = random.choice([b for b in ['A','C','G','T'] if b != mol_arr[p]])
                    
                spliced_introns = ref_data.get('spliced_introns', [])
                if spliced_introns:
                    keep_mask = np.ones(len(mol_arr), dtype=bool)
                    for start, end in spliced_introns:
                        keep_mask[start:end] = False
                        
                    mol_arr = mol_arr[keep_mask]
                    orig_idx_arr = orig_idx_arr[keep_mask]
                    is_inc = is_inc[keep_mask]
                    is_conv = is_conv[keep_mask]
            else:
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
                for p in err_pos:
                    mol_arr[p] = random.choice([b for b in ['A','C','G','T'] if b != mol_arr[p]])

            try:
                r_start = max(0, int(row.read_start) - 1) 
                r_end = min(len(mol_arr), int(row.read_end))
            except: continue
            
            if r_start >= r_end: continue

            upstream_start = r_start
            upstream_end = min(r_end, r_start + read_length)
            downstream_start = max(r_start, r_end - read_length)
            downstream_end = r_end
            
            if ttseq and not is_bg:
                if not (np.any(is_inc[upstream_start:upstream_end]) or np.any(is_inc[downstream_start:downstream_end])):
                    continue

            seq_upstream = "".join(mol_arr[upstream_start:upstream_end])
            seq_downstream = "".join(mol_arr[downstream_start:downstream_end])

            inc_upstream_rel = np.where(is_inc[upstream_start:upstream_end])[0].tolist()
            conv_upstream_rel = np.where(is_conv[upstream_start:upstream_end])[0].tolist()
            inc_downstream_rel = np.where(is_inc[downstream_start:downstream_end])[0].tolist()
            conv_downstream_rel = np.where(is_conv[downstream_start:downstream_end])[0].tolist()

            abs_u_start = orig_idx_arr[upstream_start] if upstream_start < len(orig_idx_arr) else 0
            abs_u_end = orig_idx_arr[upstream_end - 1] if upstream_end > upstream_start else 0
            coords_up = get_absolute_coords(rates_df, abs_u_start, abs_u_end)
            
            abs_d_start = orig_idx_arr[downstream_start] if downstream_start < len(orig_idx_arr) else 0
            abs_d_end = orig_idx_arr[downstream_end - 1] if downstream_end > downstream_start else 0
            coords_down = get_absolute_coords(rates_df, abs_d_start, abs_d_end)

            random_suffix = make_random_suffix()

            if sequencing_type == "SE":
                target_read = None
                if strandedness == "rf":
                    target_read = "downstream"
                    rev = True
                elif strandedness == "fr":
                    target_read = "upstream"
                    rev = False
                elif strandedness == "unstranded":
                    target_read = "downstream" if random.choice([True, False]) else "upstream"
                    rev = target_read == "downstream"
                
                if target_read == "upstream":
                    final_seq, n_inc, p_inc, n_conv, p_conv, abs_coords = seq_upstream, len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream_rel), conv_upstream_rel, coords_up
                else:
                    final_seq, n_inc, p_inc, n_conv, p_conv, abs_coords = seq_downstream, len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream_rel), conv_downstream_rel, coords_down

                final_seq_proc = process_sequence(final_seq, reverse_comp=rev)
                conv_str = ','.join(str(pos + 1) for pos in p_conv)
                inc_str = ','.join(str(pos + 1) for pos in p_inc)

                read_name = f"{ref_data['mol_id']}{random_suffix}_{abs_coords}_ninc{n_inc}:{inc_str}_nsubs{n_conv}:{conv_str}"
                f_out.write(f"{read_name}\n{final_seq_proc}\n+\n{'I' * len(final_seq_proc)}\n")

            elif sequencing_type == "PE":
                if strandedness == "rf":
                    s1, rev1, c1, m1 = seq_downstream, True, coords_down, (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream_rel), conv_downstream_rel)
                    s2, rev2, c2, m2 = seq_upstream, False, coords_up, (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream_rel), conv_upstream_rel)
                elif strandedness == "fr":
                    s1, rev1, c1, m1 = seq_upstream, False, coords_up, (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream_rel), conv_upstream_rel)
                    s2, rev2, c2, m2 = seq_downstream, True, coords_down, (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream_rel), conv_downstream_rel)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        s1, rev1, c1, m1 = seq_downstream, True, coords_down, (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream_rel), conv_downstream_rel)
                        s2, rev2, c2, m2 = seq_upstream, False, coords_up, (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream_rel), conv_upstream_rel)
                    else:
                        s1, rev1, c1, m1 = seq_upstream, False, coords_up, (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream_rel), conv_upstream_rel)
                        s2, rev2, c2, m2 = seq_downstream, True, coords_down, (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream_rel), conv_downstream_rel)

                seq_r1 = process_sequence(s1, reverse_comp=rev1)
                seq_r2 = process_sequence(s2, reverse_comp=rev2)

                name_r1 = f"{ref_data['mol_id']}{random_suffix}_{c1}_ninc{m1[0]}:{','.join(str(p+1) for p in m1[1])}_nsubs{m1[2]}:{','.join(str(p+1) for p in m1[3])}"
                name_r2 = f"{ref_data['mol_id']}{random_suffix}_{c2}_ninc{m2[0]}:{','.join(str(p+1) for p in m2[1])}_nsubs{m2[2]}:{','.join(str(p+1) for p in m2[3])}"

                f1_out.write(f"{name_r1}\n{seq_r1}\n+\n{'I'*len(seq_r1)}\n")
                f2_out.write(f"{name_r2}\n{seq_r2}\n+\n{'I'*len(seq_r2)}\n")

    finally:
        if sequencing_type == "SE":
            f_out.close()
        else:
            f1_out.close()
            f2_out.close()

def safe_eval(val):
    if pd.isna(val) or val == 'NA': return []
    if isinstance(val, str):
        try: return ast.literal_eval(val)
        except: return []
    if isinstance(val, list): return val
    return []

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")
    parser.add_argument("--insert_size", type=str, default="200,300")
    parser.add_argument("--read_length", type=int, default=100)
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE')
    parser.add_argument("--s", choices=["rf", "fr","unstranded"], type=str, default='rf')
    parser.add_argument('--seq_depth', type=int, default=20000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--tpm_lower_limit', type=int, default=5)
    parser.add_argument('--tpm_upper_limit', type=int, default=200)
    parser.add_argument("--fragments", action="store_true")
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--experiment_time", type=int, default=15)
    parser.add_argument("--bkg_molecules", type=float, default=0.0)
    parser.add_argument("--sizeselectiontype", choices=["none", "hardcut", "probabilistic"], default="probabilistic")
    parser.add_argument("--no_fragmentation", action="store_true")
    parser.add_argument("--sub_base", type=str, default='C')

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    cmd_chop = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper.R')}",
            f"--tsv {args.input_df}",
            f"--insert_size {args.insert_size}",
            f"--read_length {args.read_length}",
            f"--threads {args.threads}",
            f"--seq_depth {args.seq_depth}",
            f"--tpm_lower_limit {args.tpm_lower_limit}",
            f"--tpm_upper_limit {args.tpm_upper_limit}",
            f"-o {args.o}",
            f"--sizeselectiontype {args.sizeselectiontype}"
        ]
    if args.seed: cmd_chop += ['--seed', str(args.seed)]
    if args.fragments: cmd_chop += ['--fragments with_ground_truth']
    if args.no_fragmentation: cmd_chop += ['--no_fragmentation']

    run_cmd(" ".join(cmd_chop))
    
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"
    rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
    
    gtf_path = os.path.join(args.o, "gtf", f"{base_filename}.tsv.gz")
    if os.path.exists(gtf_path):
        df_gtf = pd.read_csv(gtf_path, sep="\t")
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
            exon_gtf = df_gtf[df_gtf['feature'] == 'exon']
            bg_ref_seq_spliced = "".join(exon_gtf['sequence'].astype(str).tolist())
            bg_ref_arr = np.array(list(bg_ref_seq_spliced), dtype='U1')
            gene_length = exon_gtf["sequence"].str.len().sum() / 1000
        else:
            bg_ref_arr = ref_arr.copy()
            gene_length = len(ref_seq_full) / 1000
    else:
        exit(1)

    df_ref_raw = pd.read_csv(args.input_df, delimiter="\t")
    clean_ref_dict = {}
    
    for col in ['converted_positions', 'incorporated_positions', 'seq_err_positions', 'spliced_introns']:
        if col in df_ref_raw.columns:
            df_ref_raw[col] = df_ref_raw[col].apply(safe_eval)
        else:
            df_ref_raw[col] = [[] for _ in range(len(df_ref_raw))]
            
    df_ref_raw['transcript_id'] = range(1, len(df_ref_raw) + 1)
    
    for _, row in df_ref_raw.iterrows():
        clean_ref_dict[row['transcript_id']] = {
            'conv': row['converted_positions'],
            'inc': row['incorporated_positions'],
            'err': row['seq_err_positions'],
            'spliced_introns': row['spliced_introns'],
            'mol_id': row['molecule_id']
        }
    del df_ref_raw 

    try:
        df_exploded = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        df_exploded = pd.DataFrame(columns=['transcript_id', 'read_start', 'read_end'])

    gene_tpm = np.random.uniform(low=int(args.tpm_lower_limit), high=int(args.tpm_upper_limit))
    seq_depth_million = args.seq_depth / 1e6
    
    total_reads_budget = int(gene_length * gene_tpm * seq_depth_million)
    
    if args.bkg_molecules > 0.0:
        reads_to_get_bg = int(total_reads_budget * args.bkg_molecules)
        reads_to_get_gene = total_reads_budget - reads_to_get_bg
        if reads_to_get_gene < 0: 
            reads_to_get_gene = 0
            reads_to_get_bg = total_reads_budget
    else:
        reads_to_get_bg = 0
        reads_to_get_gene = total_reads_budget

    if not df_exploded.empty and reads_to_get_gene > 0:
        if len(df_exploded) >= reads_to_get_gene:
            df_gene_final = df_exploded.sample(n=reads_to_get_gene, replace=False)
        else:
            remainder = reads_to_get_gene - len(df_exploded)
            if remainder > 0:
                sampled_remainder = df_exploded.sample(n=remainder, replace=True)
                df_gene_final = pd.concat([df_exploded, sampled_remainder], ignore_index=True)
            else:
                df_gene_final = df_exploded
    else:
        df_gene_final = pd.DataFrame()

    df_bg_final = pd.DataFrame()
    
    if reads_to_get_bg > 0:
        path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
        
        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            
            cmd_chop_bg = [
                f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper.R')}",
                f"--tsv {path_to_BGmRNAs}",
                f"--insert_size {args.insert_size}",
                f"--read_length {args.read_length}",
                f"--threads {args.threads}",
                f"--seq_depth {args.seq_depth}",
                f"--tpm_lower_limit {args.tpm_lower_limit}",
                f"--tpm_upper_limit {args.tpm_upper_limit}",
                f"-o {args.o}",
                f"--sizeselectiontype {args.sizeselectiontype}"
            ]
            if args.seed: cmd_chop_bg += ['--seed', str(args.seed)]
            if args.no_fragmentation: cmd_chop_bg += ['--no_fragmentation']
            
            run_cmd(" ".join(cmd_chop_bg))
            
            offset_bg = 10_000_000 
            
            df_bg_ref = pd.read_csv(path_to_BGmRNAs, delimiter="\t")
            for col in ['seq_err_positions', 'converted_positions', 'incorporated_positions', 'spliced_introns']:
                if col in df_bg_ref.columns:
                    df_bg_ref[col] = df_bg_ref[col].apply(safe_eval)
                else:
                    df_bg_ref[col] = [[] for _ in range(len(df_bg_ref))]
            
            df_bg_ref['transcript_id'] = range(1, len(df_bg_ref) + 1)
            for _, row in df_bg_ref.iterrows():
                clean_ref_dict[row['transcript_id'] + offset_bg] = {
                    'conv': row['converted_positions'],
                    'inc': row['incorporated_positions'],
                    'err': row['seq_err_positions'],
                    'spliced_introns': row['spliced_introns'],
                    'mol_id': str(row['molecule_id']) + "_BG_"
                }
            del df_bg_ref

            chopped_bg_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
            try:
                df_bg_exploded = pd.read_csv(chopped_bg_path, delimiter="\t")
                if not df_bg_exploded.empty:
                    df_bg_exploded['transcript_id'] = pd.to_numeric(df_bg_exploded['transcript_id'], errors='coerce')
                    df_bg_exploded = df_bg_exploded.dropna(subset=['transcript_id'])
                    df_bg_exploded['transcript_id'] = df_bg_exploded['transcript_id'].astype(int) + offset_bg

                    if len(df_bg_exploded) >= reads_to_get_bg:
                        df_bg_final = df_bg_exploded.sample(n=reads_to_get_bg, replace=False)
                    else:
                        remainder = reads_to_get_bg - len(df_bg_exploded)
                        sampled_remainder = df_bg_exploded.sample(n=remainder, replace=True)
                        df_bg_final = pd.concat([df_bg_exploded, sampled_remainder], ignore_index=True)
            except (FileNotFoundError, pd.errors.EmptyDataError):
                pass

    df_total = pd.concat([df_gene_final, df_bg_final], ignore_index=True)
    df_total = df_total.sample(frac=1).reset_index(drop=True)

    if os.path.exists(rates_for_gene):
        rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
    else:
        rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])

    process_and_write_fastq(
        df_total, 
        clean_ref_dict, 
        output_prefix, 
        args.seq_type, 
        args.s, 
        args.read_length, 
        rates_df, 
        ref_arr,
        bg_ref_arr,
        sub_base=args.sub_base,
        ttseq=False
    )

    try:
        if os.path.exists(chopped_coordinates_file_path):
            os.remove(chopped_coordinates_file_path)
        if 'chopped_bg_path' in locals() and os.path.exists(chopped_bg_path):
            os.remove(chopped_bg_path)
    except OSError: pass
