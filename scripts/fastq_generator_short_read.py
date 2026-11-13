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
from spark_labeling import parse_conversion_map
from spark_reads import pool_paths, write_units, read_pool, remove_paths, sample_rows, count_rows, explode_fragments, gene_references, coordinate_mapper

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

def process_and_write_fastq(
    df_sampled, lookup_dict, output_prefix, sequencing_type, 
    strandedness, read_length, rates_df, ref_arr, bg_ref_arr, conv_map=None, ttseq=False
):
    to_absolute = coordinate_mapper(rates_df, extrapolate=True)
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
            
            is_inc = np.zeros(len(mol_arr), dtype=bool)
            is_conv = np.zeros(len(mol_arr), dtype=bool)
            
            if not is_bg:
                inc_pos = [p for p in ref_data.get('inc', []) if p < len(mol_arr)]
                conv_pos = [p for p in ref_data.get('conv', []) if p < len(mol_arr)]
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
                
                is_inc[inc_pos] = True
                is_conv[conv_pos] = True
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
                    is_inc = is_inc[keep_mask]
                    is_conv = is_conv[keep_mask]
            else:
                err_pos = [p for p in ref_data.get('err', []) if p < len(mol_arr)]
                for p in err_pos:
                    mol_arr[p] = random.choice([b for b in ['A','C','G','T'] if b != mol_arr[p]])

            try:
                r_start = max(0, int(row.read_start) - 1)
                r_end = int(row.read_end)
                r_end = min(len(mol_arr), r_end)
            except: continue
            
            if r_start >= r_end: continue

            upstream_start = r_start
            upstream_end = min(r_end, r_start + read_length)
            downstream_start = max(r_start, r_end - read_length)
            downstream_end = r_end
            
            seq_upstream = "".join(mol_arr[upstream_start:upstream_end])
            seq_downstream = "".join(mol_arr[downstream_start:downstream_end])

            inc_upstream_rel = np.where(is_inc[upstream_start:upstream_end])[0].tolist()
            conv_upstream_rel = np.where(is_conv[upstream_start:upstream_end])[0].tolist()
            inc_downstream_rel = np.where(is_inc[downstream_start:downstream_end])[0].tolist()
            conv_downstream_rel = np.where(is_conv[downstream_start:downstream_end])[0].tolist()

            abs_u_start = orig_idx_arr[upstream_start] if upstream_start < len(orig_idx_arr) else 0
            abs_u_end = orig_idx_arr[upstream_end - 1] if upstream_end > upstream_start else 0
            coords_up = to_absolute(abs_u_start, abs_u_end)
            
            abs_d_start = orig_idx_arr[downstream_start] if downstream_start < len(orig_idx_arr) else 0
            abs_d_end = orig_idx_arr[downstream_end - 1] if downstream_end > downstream_start else 0
            coords_down = to_absolute(abs_d_start, abs_d_end)

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

def run_chopper(tsv_path, args, with_ground_truth):
    cmd_chop = [
        f"python {os.path.join(SCRIPT_DIR, 'short_read_chopper.py')}",
        "--mode normal",
        f"--tsv {tsv_path}",
        f"--insert_size {args.insert_size}",
        f"--read_length {args.read_length}",
        f"--threads {args.threads}",
        f"-o {args.o}",
        f"--sizeselectiontype {args.sizeselectiontype}"
    ]
    if args.seed is not None: cmd_chop += ['--seed', str(args.seed)]
    if with_ground_truth: cmd_chop += ['--fragments with_ground_truth']
    if args.no_fragmentation: cmd_chop += ['--no_fragmentation']
    run_cmd(" ".join(cmd_chop))

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
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--fragments", action="store_true")
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--seed", type=int)
    parser.add_argument("--experiment_time", type=int, default=15)
    parser.add_argument("--bkg_molecules", type=float, default=0.0)
    parser.add_argument("--sizeselectiontype", choices=["none", "hardcut", "probabilistic"], default="probabilistic")
    parser.add_argument("--no_fragmentation", action="store_true")
    parser.add_argument("--labeling_base", type=str, default='T')
    parser.add_argument("--sub_base", type=str, default='C')
    parser.add_argument("--experiment_type", type=str, default="metaboliclabeling")
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
    rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
    path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
    nascent_pool_path, bg_pool_path, units_path = pool_paths(args.o, base_filename)
    offset_bg = 10_000_000

    gtf_path = os.path.join(args.o, "gtf", f"{base_filename}.tsv.gz")
    if not os.path.exists(gtf_path):
        exit(1)

    if args.prepare:
        run_chopper(args.input_df, args, args.fragments)
        nascent_pool = explode_fragments(chopped_coordinates_file_path)
        if args.experiment_type != 'chromatin_associated' and not nascent_pool.empty:
            incorporated = pd.read_csv(args.input_df, delimiter="\t", usecols=['incorporated_positions'])['incorporated_positions'].apply(safe_eval)
            labeled_ids = np.flatnonzero(incorporated.map(len).values > 0) + 1
            nascent_pool = nascent_pool[nascent_pool['transcript_id'].isin(labeled_ids)]

        bg_pool = pd.DataFrame(columns=nascent_pool.columns)
        bg_molecules = 0
        if args.bkg_molecules > 0.0 and count_rows(path_to_BGmRNAs) > 0:
            run_chopper(path_to_BGmRNAs, args, False)
            bg_pool = explode_fragments(chopped_bg_path)
            bg_molecules = count_rows(path_to_BGmRNAs)

        nascent_pool.to_csv(nascent_pool_path, sep="\t", index=False)
        bg_pool.to_csv(bg_pool_path, sep="\t", index=False)
        write_units(units_path, base_filename, len(nascent_pool), len(bg_pool), bg_molecules)
        remove_paths([chopped_coordinates_file_path, chopped_bg_path])
        sys.exit(0)

    ref_arr, bg_ref_arr = gene_references(gtf_path, path_to_BGmRNAs)

    df_ref_raw = read_pool(args.input_df)
    clean_ref_dict = {}
    
    for col in ['converted_positions', 'incorporated_positions', 'seq_err_positions', 'spliced_introns']:
        if col in df_ref_raw.columns:
            df_ref_raw[col] = df_ref_raw[col].apply(safe_eval)
        else:
            df_ref_raw[col] = [[] for _ in range(len(df_ref_raw))]
            
    df_ref_raw['transcript_id'] = range(1, len(df_ref_raw) + 1)
    
    for _, row in df_ref_raw.iterrows():
        #if the molecule incorporated zero analogs, it does not bind streptavidin beas or whatever the method is doing
        #if the molecules come from chromatin fractionation there is no need for nucleotide analogs
        if args.experiment_type != 'chromatin_associated':
            if not row['incorporated_positions']: 
                continue
            
        clean_ref_dict[row['transcript_id']] = {
            'conv': row['converted_positions'],
            'inc': row['incorporated_positions'],
            'err': row['seq_err_positions'],
            'spliced_introns': row['spliced_introns'],
            'mol_id': row['molecule_id'],
        }
    del df_ref_raw 

    df_gene_final = sample_rows(read_pool(nascent_pool_path), args.nascent_reads)
    df_bg_final = sample_rows(read_pool(bg_pool_path), args.bkg_reads)

    if not df_bg_final.empty:
        df_bg_final = df_bg_final.assign(transcript_id=df_bg_final['transcript_id'] + offset_bg)
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
                'mol_id': str(row['molecule_id']) + "_BG_",
                'start_offset': row['start_label_pos'] if 'start_label_pos' in row else 0
            }
        del df_bg_ref

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
        conv_map=conv_map,
        ttseq=False
    )

    remove_paths([nascent_pool_path, bg_pool_path, units_path])
