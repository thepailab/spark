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

def safe_eval(val):
    if pd.isna(val) or val == 'NA': return []
    if isinstance(val, str):
        try: return ast.literal_eval(val)
        except: return []
    if isinstance(val, list): return val
    return []

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
    return "NA"

def get_valid_ttseq_fragments(df_frags, lookup_dict, read_length):
    valid_indices = []
    for row in df_frags.itertuples():
        ref = lookup_dict.get(row.transcript_id)
        if not ref: 
            continue
        
        rel_inc = ref['rel_inc']
        if len(rel_inc) == 0:
            continue
            
        try:
            r_start, r_end = map(int, str(row.read_coordinate_split).split('-'))
            r_start = max(0, r_start - 1)
        except: continue
        
        u_start = r_start
        u_end = r_start + read_length
        d_start = max(r_start, r_end - read_length)
        d_end = r_end
        
        has_inc = np.any((rel_inc >= u_start) & (rel_inc < u_end)) or \
                  np.any((rel_inc >= d_start) & (rel_inc < d_end))
        
        if has_inc:
            valid_indices.append(row.Index)
            
    return df_frags.loc[valid_indices].copy()

def process_and_write_fastq(
    df_sampled, lookup_dict, output_prefix, sequencing_type, 
    strandedness, read_length, rates_df, ref_arr, bg_ref_arr, sub_base='C'
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
            
            if not is_bg:
                conv_pos = [p for p in ref_data['conv'] if p < len(mol_arr)]
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
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
            else:
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
                for p in err_pos:
                    mol_arr[p] = random.choice([b for b in ['A','C','G','T'] if b != mol_arr[p]])

            try:
                r_start, r_end = map(int, str(row.read_coordinate_split).split('-'))
                r_start = max(0, r_start - 1) 
                r_end = min(r_end, len(mol_arr))
            except: continue
            if r_start >= r_end: continue

            upstream_start = r_start
            upstream_end = r_start + read_length
            downstream_start = max(r_start, r_end - read_length)
            downstream_end = r_end
            
            seq_upstream = "".join(mol_arr[upstream_start:min(upstream_end, len(mol_arr))])
            if len(seq_upstream) < read_length:
                seq_upstream += "N" * (read_length - len(seq_upstream))

            seq_downstream = "".join(mol_arr[downstream_start:downstream_end])
            if len(seq_downstream) < read_length:
                padding = "N" * (read_length - len(seq_downstream))
                if downstream_start == 0:
                    seq_downstream = padding + seq_downstream
                else:
                    seq_downstream += padding

            abs_u_start = orig_idx_arr[upstream_start] if upstream_start < len(orig_idx_arr) else 0
            abs_u_end = orig_idx_arr[min(upstream_end, len(orig_idx_arr)) - 1] if upstream_end > upstream_start else 0
            coords_up = get_absolute_coords(rates_df, abs_u_start, abs_u_end)
            
            abs_d_start = orig_idx_arr[downstream_start] if downstream_start < len(orig_idx_arr) else 0
            abs_d_end = orig_idx_arr[min(downstream_end, len(orig_idx_arr)) - 1] if downstream_end > downstream_start else 0
            coords_down = get_absolute_coords(rates_df, abs_d_start, abs_d_end)

            random_suffix = make_random_suffix()

            if sequencing_type == "SE":
                if strandedness == "rf":
                    final_seq, rev, coords = seq_downstream, True, coords_down
                elif strandedness == "fr":
                    final_seq, rev, coords = seq_upstream, False, coords_up
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        final_seq, rev, coords = seq_downstream, True, coords_down
                    else:
                        final_seq, rev, coords = seq_upstream, False, coords_up
                
                final_seq_proc = process_sequence(final_seq, reverse_comp=rev)
                read_name = f"{ref_data['mol_id']}{random_suffix}_{coords}"
                f_out.write(f"{read_name}\n{final_seq_proc}\n+\n{'I' * len(final_seq_proc)}\n")

            elif sequencing_type == "PE":
                if strandedness == "rf":
                    s1, rev1, c1 = seq_downstream, True, coords_down
                    s2, rev2, c2 = seq_upstream, False, coords_up
                elif strandedness == "fr":
                    s1, rev1, c1 = seq_upstream, False, coords_up
                    s2, rev2, c2 = seq_downstream, True, coords_down
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        s1, rev1, c1 = seq_downstream, True, coords_down
                        s2, rev2, c2 = seq_upstream, False, coords_up
                    else:
                        s1, rev1, c1 = seq_upstream, False, coords_up
                        s2, rev2, c2 = seq_downstream, True, coords_down

                seq_r1 = process_sequence(s1, reverse_comp=rev1)
                seq_r2 = process_sequence(s2, reverse_comp=rev2)

                name_r1 = f"{ref_data['mol_id']}{random_suffix}_{c1}"
                name_r2 = f"{ref_data['mol_id']}{random_suffix}_{c2}"

                f1_out.write(f"{name_r1}\n{seq_r1}\n+\n{'I'*len(seq_r1)}\n")
                f2_out.write(f"{name_r2}\n{seq_r2}\n+\n{'I'*len(seq_r2)}\n")

    finally:
        if sequencing_type == "SE":
            f_out.close()
        else:
            f1_out.close()
            f2_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")
    parser.add_argument("--insert_size", type=str, default="200,300")
    parser.add_argument("--read_length", type=int,default=100)
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE')
    parser.add_argument("--s", choices=["rf", "fr","unstranded"], type=str, default='rf')
    parser.add_argument('--seq_depth', type=int, default=20000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--tpm_lower_limit', type=int, default=5)
    parser.add_argument('--tpm_upper_limit', type=int, default=200)
    parser.add_argument("--fragments", action="store_true")
    parser.add_argument('--bkg_molecules', type=float, default=0)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--sizeselectiontype", choices=["none", "hardcut", "probabilistic"], default="probabilistic")
    parser.add_argument("--no_fragmentation", action="store_true")
    parser.add_argument("--sub_base", type=str, default='C')

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    cmd_chop = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper_TTseq.R')}",
            f"--tsv {args.input_df}",
            f"--insert_size {args.insert_size}",
            f"--read_length {args.read_length}",
            f"--threads {args.threads}",
            f"-o {args.o}",
            f"--sizeselectiontype {args.sizeselectiontype}"
        ]
    if args.seed: cmd_chop += ['--seed', str(args.seed)]
    if args.no_fragmentation: cmd_chop += ['--no_fragmentation']

    run_cmd(" ".join(cmd_chop))

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"

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
        is_inc = np.zeros(len(ref_arr), dtype=bool)
        valid_incs = [p for p in row['incorporated_positions'] if p < len(ref_arr)]
        is_inc[valid_incs] = True
        
        spliced_introns = row.get('spliced_introns', [])
        if spliced_introns:
            for start, end in sorted(spliced_introns, key=lambda x: x[0], reverse=True):
                is_inc = np.delete(is_inc, slice(start, end))
                
        clean_ref_dict[row['transcript_id']] = {
            'conv': row['converted_positions'],
            'inc': row['incorporated_positions'],
            'err': row['seq_err_positions'],
            'spliced_introns': spliced_introns,
            'mol_id': row['molecule_id'],
            'rel_inc': np.where(is_inc)[0]
        }
    del df_ref_raw 

    try: 
        df = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
        df['read_coordinates'] = df['read_coordinates'].astype(str).str.split(',')
        df_exploded = df.explode('read_coordinates').reset_index(drop=True)
        df_exploded.rename(columns={'read_coordinates': 'read_coordinate_split'}, inplace=True)
        
        coord_split = df_exploded['read_coordinate_split'].str.split('-', expand=True)
        df_exploded['read_start'] = pd.to_numeric(coord_split[0], errors='coerce')
        df_exploded['read_end'] = pd.to_numeric(coord_split[1], errors='coerce')
        df_exploded = df_exploded.dropna(subset=['read_start', 'read_end'])
        
        result_df = get_valid_ttseq_fragments(df_exploded, clean_ref_dict, args.read_length)

        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            cmd_chop_bg = [
                f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper_TTseq.R')}",
                f"--tsv {path_to_BGmRNAs}",
                f"--insert_size {args.insert_size}",
                f"--read_length {args.read_length}",
                f"--threads {args.threads}",
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
                df_bg = pd.read_csv(chopped_bg_path, delimiter="\t")
                if not df_bg.empty:
                    df_bg['transcript_id'] = pd.to_numeric(df_bg['transcript_id'], errors='coerce')
                    df_bg = df_bg.dropna(subset=['transcript_id'])
                    df_bg['transcript_id'] = df_bg['transcript_id'].astype(int) + offset_bg
                    df_bg['read_coordinates'] = df_bg['read_coordinates'].astype(str).str.split(',')
                    result_df_bg = df_bg.explode('read_coordinates').reset_index(drop=True)
                    result_df_bg.rename(columns={'read_coordinates': 'read_coordinate_split'}, inplace=True)
                    
                    coord_split_bg = result_df_bg['read_coordinate_split'].str.split('-', expand=True)
                    result_df_bg['read_start'] = pd.to_numeric(coord_split_bg[0], errors='coerce')
                    result_df_bg['read_end'] = pd.to_numeric(coord_split_bg[1], errors='coerce')
                    result_df_bg = result_df_bg.dropna(subset=['read_start', 'read_end'])
                    
                    n_bg = int(len(result_df) * args.bkg_molecules)
                    if n_bg > 0:
                        if len(result_df_bg) >= n_bg:
                            result_df_bg = result_df_bg.sample(n=n_bg, random_state=42)
                        else:
                            remainder = n_bg - len(result_df_bg)
                            sampled_remainder = result_df_bg.sample(n=remainder, replace=True)
                            result_df_bg = pd.concat([result_df_bg, sampled_remainder], ignore_index=True)
                else:
                    result_df_bg = pd.DataFrame()
            except:
                result_df_bg = pd.DataFrame()
            
            n_main = int(len(result_df) * (1 - args.bkg_molecules))
            result_df = result_df.sample(n=n_main, random_state=42) 
            result_df = pd.concat([result_df, result_df_bg]).reset_index(drop=True)

        result_df["insert_size"] = result_df["read_end"] - result_df["read_start"]
        insert_values = args.insert_size.split(",")
        min_insert = int(insert_values[0])
        max_insert = int(insert_values[1])

        result_df_pre_selection = result_df
        
        if args.sizeselectiontype == 'hardcut':
            result_df = result_df[(result_df['insert_size'] >= min_insert) & (result_df['insert_size'] <= max_insert)]

        gene_tpm = np.random.uniform(low=int(args.tpm_lower_limit), high=int(args.tpm_upper_limit))
        seq_depth = args.seq_depth / 1e6
        reads_to_get = gene_length * gene_tpm * seq_depth

        if args.fragments:
            if not result_df_pre_selection.empty:
                n_sample_pre = min(len(result_df_pre_selection), int(reads_to_get))
                replace_pre = n_sample_pre > len(result_df_pre_selection)
                result_df_pre_selection_sampled = result_df_pre_selection.sample(n=int(reads_to_get), replace=True, random_state=None)
                
                output_dir_pre = os.path.join(args.o, 'ground_truth_pre_size_selection')
                os.makedirs(output_dir_pre, exist_ok=True)
                max_position_pre = result_df_pre_selection_sampled['read_end'].max()
                if not pd.isna(max_position_pre):
                    coverage_pre = np.zeros(int(max_position_pre) + 2, dtype=int)
                    np.add.at(coverage_pre, result_df_pre_selection_sampled['read_start'].astype(int), 1)
                    np.add.at(coverage_pre, result_df_pre_selection_sampled['read_end'].astype(int), -1)
                    coverage_pre = np.cumsum(coverage_pre[:-1])
                    coverage_df_pre = pd.DataFrame({'position': np.arange(len(coverage_pre)),'frequency': coverage_pre})
                    coverage_output_path_pre = os.path.join(output_dir_pre, f"{base_filename}.tsv.gz")
                    coverage_df_pre.to_csv(coverage_output_path_pre, sep='\t', index=False, compression='gzip')


            if not result_df.empty:
                result_df_sampled = result_df.sample(n=int(reads_to_get), replace=True, random_state=None)
                output_dir_post = os.path.join(args.o, 'ground_truth_after_size_selection')
                os.makedirs(output_dir_post, exist_ok=True)
                max_position_post = result_df_sampled['read_end'].max()
                if not pd.isna(max_position_post):
                    coverage_post = np.zeros(int(max_position_post) + 2, dtype=int)
                    np.add.at(coverage_post, result_df_sampled['read_start'].astype(int), 1)
                    np.add.at(coverage_post, result_df_sampled['read_end'].astype(int), -1)
                    coverage_post = np.cumsum(coverage_post[:-1])
                    coverage_df_post = pd.DataFrame({
                        'position': np.arange(len(coverage_post)),
                        'frequency': coverage_post
                    })
                    coverage_output_path_post = os.path.join(output_dir_post, f"{base_filename}.tsv.gz")
                    coverage_df_post.to_csv(coverage_output_path_post, sep='\t', index=False, compression='gzip')


        if not result_df.empty:
            if len(result_df) >= reads_to_get:
                result_df = result_df.sample(
                    n=int(reads_to_get), replace=False, random_state=None
                    )
            else:
                sampled_reads = result_df.sample(
                    n=int(reads_to_get) - len(result_df), replace=True, random_state=None
                    )
                result_df = pd.concat([result_df, sampled_reads], ignore_index=True)

            rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
            if os.path.exists(rates_for_gene):
                rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
            else:
                rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])

            process_and_write_fastq(
                result_df, clean_ref_dict, output_prefix, args.seq_type, args.s,
                args.read_length, rates_df, ref_arr, bg_ref_arr, sub_base=args.sub_base
            )

        try:
            if os.path.exists(chopped_coordinates_file_path):
                os.remove(chopped_coordinates_file_path)
            if 'chopped_bg_path' in locals() and os.path.exists(chopped_bg_path):
                os.remove(chopped_bg_path)
        except OSError:
            pass

    except FileNotFoundError:
        pass
