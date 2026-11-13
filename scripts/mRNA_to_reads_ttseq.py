import pandas as pd
import gzip
import random
import argparse
import os
import sys
import subprocess
import string
import ast
import numpy as np
import hashlib
from spark_labeling import parse_conversion_map
from spark_reads import pool_paths, write_units, read_pool, remove_paths, sample_rows, count_rows, explode_fragments, gene_references, coordinate_mapper

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

def get_valid_ttseq_fragments(df_frags, lookup_dict, read_length):
    valid_indices = []
    for row in df_frags.itertuples():
        ref = lookup_dict.get(row.transcript_id)
        if not ref: 
            continue
        
        rel_inc = ref['rel_inc']
        if len(rel_inc) == 0:
            continue
            
        r_start = max(0, int(row.read_start) - 1)
        r_end = int(row.read_end)
        
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
    strandedness, read_length, rates_df, ref_arr, bg_ref_arr, conv_map=None
):
    to_absolute = coordinate_mapper(rates_df)
    conv_map = conv_map or {}
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
                if conv_map and conv_pos:
                    orig_bases = mol_arr[conv_pos]
                    mol_arr[conv_pos] = [conv_map.get(b, b) for b in orig_bases]
                
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

            r_start = max(0, int(row.read_start) - 1)
            r_end = min(int(row.read_end), len(mol_arr))
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
            coords_up = to_absolute(abs_u_start, abs_u_end)
            
            abs_d_start = orig_idx_arr[downstream_start] if downstream_start < len(orig_idx_arr) else 0
            abs_d_end = orig_idx_arr[min(downstream_end, len(orig_idx_arr)) - 1] if downstream_end > downstream_start else 0
            coords_down = to_absolute(abs_d_start, abs_d_end)

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

def run_chopper(tsv_path, args):
    cmd_chop = [
        f"python {os.path.join(SCRIPT_DIR, 'short_read_chopper.py')}",
        "--mode ttseq",
        f"--tsv {tsv_path}",
        f"--insert_size {args.insert_size}",
        f"--read_length {args.read_length}",
        f"--threads {args.threads}",
        f"-o {args.o}",
        f"--sizeselectiontype {args.sizeselectiontype}"
    ]
    if args.seed is not None: cmd_chop += ['--seed', str(args.seed)]
    if args.no_fragmentation: cmd_chop += ['--no_fragmentation']
    run_cmd(" ".join(cmd_chop))

def build_reference_dict(input_path, ref_length, with_label_index, id_offset=0, mol_suffix=""):
    df_ref_raw = read_pool(input_path)
    for col in ['converted_positions', 'incorporated_positions', 'seq_err_positions', 'spliced_introns']:
        if col in df_ref_raw.columns:
            df_ref_raw[col] = df_ref_raw[col].apply(safe_eval)
        else:
            df_ref_raw[col] = [[] for _ in range(len(df_ref_raw))]

    df_ref_raw['transcript_id'] = range(1, len(df_ref_raw) + 1)

    clean_ref_dict = {}
    for _, row in df_ref_raw.iterrows():
        entry = {
            'conv': row['converted_positions'],
            'inc': row['incorporated_positions'],
            'err': row['seq_err_positions'],
            'spliced_introns': row['spliced_introns'],
            'mol_id': str(row['molecule_id']) + mol_suffix,
        }
        if with_label_index:
            is_inc = np.zeros(ref_length, dtype=bool)
            valid_incs = [p for p in row['incorporated_positions'] if p < ref_length]
            is_inc[valid_incs] = True
            spliced_introns = row.get('spliced_introns', [])
            if spliced_introns:
                for start, end in sorted(spliced_introns, key=lambda x: x[0], reverse=True):
                    is_inc = np.delete(is_inc, slice(start, end))
            entry['rel_inc'] = np.where(is_inc)[0]
        clean_ref_dict[row['transcript_id'] + id_offset] = entry
    return clean_ref_dict

def size_selected(df, args):
    if args.sizeselectiontype != 'hardcut' or df.empty:
        return df
    min_insert, max_insert = map(int, args.insert_size.split(","))
    insert_size = df['read_end'] - df['read_start']
    return df[(insert_size >= min_insert) & (insert_size <= max_insert)]

def write_coverage(df, output_dir, base_filename):
    if df.empty:
        return
    os.makedirs(output_dir, exist_ok=True)
    max_position = df['read_end'].max()
    coverage = np.zeros(int(max_position) + 2, dtype=int)
    np.add.at(coverage, df['read_start'].astype(int), 1)
    np.add.at(coverage, df['read_end'].astype(int), -1)
    coverage = np.cumsum(coverage[:-1])
    coverage_df = pd.DataFrame({'position': np.arange(len(coverage)), 'frequency': coverage})
    coverage_df.to_csv(os.path.join(output_dir, f"{base_filename}.tsv.gz"), sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df")
    parser.add_argument("--insert_size", type=str, default="200,300")
    parser.add_argument("--read_length", type=int,default=100)
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE')
    parser.add_argument("--s", choices=["rf", "fr","unstranded"], type=str, default='rf')
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--fragments", action="store_true")
    parser.add_argument('--bkg_molecules', type=float, default=0)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--sizeselectiontype", choices=["none", "hardcut", "probabilistic"], default="probabilistic")
    parser.add_argument("--no_fragmentation", action="store_true")
    parser.add_argument("--labeling_base", type=str, default='T')
    parser.add_argument("--sub_base", type=str, default='C')
    parser.add_argument("--prepare", action="store_true", help="build the fragment pools and count them instead of writing reads")
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

    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"
    chopped_bg_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
    path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
    nascent_pool_path, bg_pool_path, units_path = pool_paths(args.o, base_filename)
    offset_bg = 10_000_000

    gtf_path = os.path.join(args.o, "gtf", f"{base_filename}.tsv.gz")
    if not os.path.exists(gtf_path):
        exit(1)
    ref_arr, bg_ref_arr = gene_references(gtf_path, path_to_BGmRNAs)

    if args.prepare:
        run_chopper(args.input_df, args)
        nascent_pool = explode_fragments(chopped_coordinates_file_path)
        if not nascent_pool.empty:
            nascent_pool = get_valid_ttseq_fragments(
                nascent_pool, build_reference_dict(args.input_df, len(ref_arr), True), args.read_length
            )

        bg_pool = pd.DataFrame(columns=nascent_pool.columns)
        bg_molecules = 0
        if args.bkg_molecules > 0.0 and count_rows(path_to_BGmRNAs) > 0:
            run_chopper(path_to_BGmRNAs, args)
            bg_pool = explode_fragments(chopped_bg_path)
            bg_molecules = count_rows(path_to_BGmRNAs)

        nascent_pool.to_csv(nascent_pool_path, sep="\t", index=False)
        bg_pool.to_csv(bg_pool_path, sep="\t", index=False)
        write_units(units_path, base_filename, len(size_selected(nascent_pool, args)), len(size_selected(bg_pool, args)), bg_molecules)
        remove_paths([chopped_coordinates_file_path, chopped_bg_path])
        sys.exit(0)

    clean_ref_dict = build_reference_dict(args.input_df, len(ref_arr), False)
    nascent_pool = read_pool(nascent_pool_path)
    bg_pool = read_pool(bg_pool_path)
    if not bg_pool.empty:
        bg_pool = bg_pool.assign(transcript_id=bg_pool['transcript_id'] + offset_bg)
        clean_ref_dict.update(build_reference_dict(path_to_BGmRNAs, len(bg_ref_arr), False, offset_bg, "_BG_"))

    result_df = pd.concat([
        sample_rows(size_selected(nascent_pool, args), args.nascent_reads),
        sample_rows(size_selected(bg_pool, args), args.bkg_reads),
    ], ignore_index=True)

    if args.fragments:
        result_df_pre_selection = pd.concat([
            sample_rows(nascent_pool, args.nascent_reads),
            sample_rows(bg_pool, args.bkg_reads),
        ], ignore_index=True)
        write_coverage(result_df_pre_selection, os.path.join(args.o, 'ground_truth_pre_size_selection'), base_filename)
        write_coverage(result_df, os.path.join(args.o, 'ground_truth_after_size_selection'), base_filename)

    if not result_df.empty:
        rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
        if os.path.exists(rates_for_gene):
            rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
        else:
            rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])

        process_and_write_fastq(
            result_df, clean_ref_dict, output_prefix, args.seq_type, args.s,
            args.read_length, rates_df, ref_arr, bg_ref_arr, conv_map=conv_map
        )

    remove_paths([nascent_pool_path, bg_pool_path, units_path])
