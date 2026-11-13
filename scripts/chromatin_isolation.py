import pandas as pd
import numpy as np
import argparse
import os
import hashlib
import shutil
import sys

def main():
    parser = argparse.ArgumentParser(description="keep transcripts that are chromatin associated")
    parser.add_argument('--input_df', required=True, help="path to the mRNA dataframe")
    parser.add_argument('--release_half_life', required=True, type=str, help="range of chromatin release half life in minutes (e.g., '2.0,5.0')")
    parser.add_argument('--seed', type=int, default=None, help="random seed")
    parser.add_argument('--o', type=str, default='./', help="output directory")
    parser.add_argument('--reuse_parameters', action='store_true', help="reuse the gene's stored chromatin release half-life instead of drawing a new one")
    
    args = parser.parse_args()
    
    os.makedirs(args.o, exist_ok=True)
    out_file = os.path.join(args.o, os.path.basename(args.input_df))
    
    #if the mRNA file is background molecules, just copy, I don't want to filter them
    if "_background.tsv.gz" in args.input_df:
        shutil.copy2(args.input_df, out_file)
        sys.exit(0)
    
    base_filename = os.path.basename(args.input_df).replace(".tsv.gz", "")
    
    #generate seed based on gene name
    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % (2**32)
        np.random.seed(unique_seed)
        
    df = pd.read_csv(args.input_df, sep='\t')
    
    parent_dir = os.path.dirname(os.path.abspath(args.o.rstrip('/')))
    parameters_path = os.path.join(parent_dir, 'rate_per_gene', f"{base_filename}_parameters.tsv")
    stored_parameters = pd.read_csv(parameters_path, sep='\t') if os.path.exists(parameters_path) else pd.DataFrame([{'gene_id': base_filename}])

    #get gene-specific chromatin release half life
    if args.reuse_parameters and 'chromatin_release_half_life_min' in stored_parameters.columns:
        gene_hl = float(stored_parameters['chromatin_release_half_life_min'].iloc[0])
    else:
        hl_min, hl_max = map(float, args.release_half_life.split(','))
        gene_hl = np.random.uniform(hl_min, hl_max)
        stored_parameters['chromatin_release_half_life_min'] = gene_hl
        os.makedirs(os.path.dirname(parameters_path), exist_ok=True)
        stored_parameters.to_csv(parameters_path, sep='\t', index=False)
    
    
    is_uncleaved = df['molecule_id'].str.endswith('_uncleaved')
    is_tail = df['molecule_id'].str.endswith('_readthrough_tail')
    is_body = df['molecule_id'].str.endswith('_cleaved_body')
    
    uncleaved_df = df[is_uncleaved].copy()
    tail_df = df[is_tail].copy()
    body_df = df[is_body].copy()
    
    #apply the half life ONLY to the _cleaved_body molecules
    if gene_hl > 0 and not body_df.empty:
        mean_life = gene_hl / np.log(2)
        p_retained = np.exp(-body_df['time_since_cleaved_min'].values / mean_life)
        retained_bool = np.random.rand(len(body_df)) < p_retained
        retained_body_df = body_df[retained_bool].copy()
    else:
        retained_body_df = pd.DataFrame(columns=df.columns)
        

    final_caRNA_df = pd.concat([uncleaved_df, tail_df, retained_body_df], ignore_index=True)
    
    final_caRNA_df.to_csv(out_file, sep='\t', index=False, compression='gzip')
    
    #export ground truth to the central temp folder for aggregation
    temp_dir = os.path.join(parent_dir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    log_file = os.path.join(temp_dir, f"temp_{base_filename}_chromatin_halflife.tsv")
    pd.DataFrame([{
        'gene_id': base_filename, 
        'chromatin_release_halflife_min': gene_hl
    }]).to_csv(log_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
