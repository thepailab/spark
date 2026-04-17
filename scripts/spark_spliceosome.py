import math
import numpy as np
import pandas as pd
import argparse
import os
from scipy.stats import linregress
import hashlib

def calculate_time_sums_optimized(mRNA_df, intron_df, mRNA_rates, strand):
    mRNA_rates_sorted = mRNA_rates.sort_values('nucleotide_coord').reset_index(drop=True)
    coords = mRNA_rates_sorted['nucleotide_coord'].values
    times = mRNA_rates_sorted['time_for_this_nt'].values
    
    # Pre-calculate cumulative sum of times for O(1) range summation
    cum_times = np.concatenate(([0], np.cumsum(times)))
    
    results = []
    initiation_times = mRNA_df['initiation_time'].values
    start_positions = mRNA_df['start_label_pos'].values
    stop_positions = mRNA_df['stop_label_pos'].values
    
    if strand == '+':
        splice_sites_3 = intron_df['end'].values
    else:
        splice_sites_3 = intron_df['start'].values

    for i in range(len(mRNA_df)):
        init_time = initiation_times[i]
        start_pos = start_positions[i]
        stop_pos = stop_positions[i]
        
        row_res = []
        
        for ss3 in splice_sites_3:
            if stop_pos < ss3:
                row_res.append(np.nan)
                continue

            idx_ss3 = np.searchsorted(coords, ss3, side='left')
            idx_start = np.searchsorted(coords, start_pos, side='left')
            
            if ss3 < start_pos:
                if idx_ss3 < idx_start:
                    extra_time = cum_times[idx_start] - cum_times[idx_ss3]
                    row_res.append(init_time + extra_time)
                else:
                    row_res.append(init_time)
            else:
                idx_ss3_right = np.searchsorted(coords, ss3, side='right')
                if idx_start < idx_ss3_right:
                    time_to_splice = cum_times[idx_ss3_right] - cum_times[idx_start]
                    row_res.append(init_time - time_to_splice)
                else:
                    row_res.append(init_time)

        results.append(row_res)
    return results

def determine_splicing_constant_half_life(t_available, half_life_value):
    n = len(t_available)
    valid_mask = ~np.isnan(t_available)
    
    spliced = np.zeros(n, dtype=bool)
    used_half_lives = np.full(n, np.nan)
    
    if np.any(valid_mask):
        valid_t = t_available[valid_mask]
        n_valid = len(valid_t)
        
        used_half_lives[valid_mask] = half_life_value
        mean_life = half_life_value / np.log(2)
        p_spliced = 1 - np.exp(-valid_t / mean_life)
        
        spliced[valid_mask] = np.random.rand(n_valid) < p_spliced
        
    return spliced, used_half_lives

def run_splicing_per_exon_vary_introns(time_sums_df, half_life_range, mean_half_life_dict=None, seed=None):
    spliced_results = {}
    half_lives_results = {}
    this_run_half_life_dict = {}
    
    rng = np.random.default_rng(seed)
    
    for exon_col in time_sums_df.columns:
        if mean_half_life_dict is not None:
            half_life_val = mean_half_life_dict[exon_col]
        else:
            half_life_val = rng.uniform(half_life_range[0], half_life_range[1])
            
        this_run_half_life_dict[exon_col] = half_life_val
        t_available = time_sums_df[exon_col].values
        
        spliced, used_half_lives = determine_splicing_constant_half_life(t_available, half_life_val)
            
        spliced_results[exon_col] = spliced
        half_lives_results[exon_col] = used_half_lives
        
    spliced_df = pd.DataFrame(spliced_results, index=time_sums_df.index)
    half_lives_df = pd.DataFrame(half_lives_results, index=time_sums_df.index)
    return spliced_df, half_lives_df, this_run_half_life_dict

def record_splicing_events(mRNA_df, spliced_df, intron_df, strand):
    """
    Replaces string manipulation by mathematically tracking the new sequence length 
    and recording the 0-based relative intervals of the introns that were spliced out.
    """
    spliced_matrix = spliced_df.values
    intron_starts = intron_df['start'].values
    intron_ends = intron_df['end'].values
    stop_positions = mRNA_df['stop_label_pos'].values
    
    sequence_lengths = np.zeros(len(mRNA_df), dtype=int)
    spliced_introns_list = []
    
    for i in range(len(mRNA_df)):
        stop_pos = stop_positions[i]
        length = stop_pos
        removed_intervals = []
        
        row_spliced = spliced_matrix[i]
        for j, is_spliced in enumerate(row_spliced):
            if is_spliced:
                start = intron_starts[j]
                end = intron_ends[j]
                
                # Map genomic offsets to 0-based sequence indices
                if strand == '+':
                    rs, re = start - 1, end
                else:
                    rs, re = end, start + 1
                    
                # Cap the removed interval by the transcription stop position
                rs = max(0, min(rs, stop_pos))
                re = max(0, min(re, stop_pos))
                
                removed_len = re - rs
                if removed_len > 0:
                    length -= removed_len
                    removed_intervals.append((int(rs), int(re)))
                    
        sequence_lengths[i] = length
        spliced_introns_list.append(str(removed_intervals))
        
    mRNA_df_spliced = mRNA_df.copy()
    mRNA_df_spliced['sequence_length'] = sequence_lengths
    mRNA_df_spliced['spliced_introns'] = spliced_introns_list
    return mRNA_df_spliced

def get_surviving_features(mRNA_df, spliced_df, mRNA_gtf, strand, TSS_coord):
    features = mRNA_gtf[mRNA_gtf['feature'].isin(['exon', 'intron'])].copy()
    
    if strand == '+':
        features['rel_start'] = (features['start'] - TSS_coord) + 1
        features['rel_end'] = (features['end'] - TSS_coord) + 1
    else:
        features['rel_start'] = (TSS_coord - features['end']) + 1
        features['rel_end'] = (TSS_coord - features['start']) + 1

    exons = features[features['feature'] == 'exon'].sort_values('position')
    introns = features[features['feature'] == 'intron'].sort_values('position')

    exon_starts = exons['rel_start'].values
    exon_ends = exons['rel_end'].values
    exon_nums = exons['position'].values
    
    intron_starts = introns['rel_start'].values
    intron_ends = introns['rel_end'].values
    intron_nums = introns['position'].values

    mol_stops = mRNA_df['stop_label_pos'].values[:, np.newaxis]
    mol_indices = mRNA_df.index.values

    exon_mask = exon_starts[np.newaxis, :] <= mol_stops
    kept_exon_rows, kept_exon_cols = np.nonzero(exon_mask)
    
    current_mol_stops_ex = mol_stops[kept_exon_rows, 0]
    current_exon_ends = exon_ends[kept_exon_cols]
    final_exon_ends = np.minimum(current_exon_ends, current_mol_stops_ex)
    
    df_exons = pd.DataFrame({
        'molecule_index': mol_indices[kept_exon_rows],
        'start': exon_starts[kept_exon_cols], 
        'end': final_exon_ends,
        'feature': 'exon',
        'feature_number': exon_nums[kept_exon_cols]
    })

    is_spliced_matrix = spliced_df.values
    intron_transcribed_mask = intron_starts[np.newaxis, :] <= mol_stops
    intron_mask = intron_transcribed_mask & (~is_spliced_matrix)
    
    kept_intron_rows, kept_intron_cols = np.nonzero(intron_mask)
    
    current_mol_stops_int = mol_stops[kept_intron_rows, 0]
    current_intron_ends = intron_ends[kept_intron_cols]
    final_intron_ends = np.minimum(current_intron_ends, current_mol_stops_int)
    
    df_introns = pd.DataFrame({
        'molecule_index': mol_indices[kept_intron_rows],
        'start': intron_starts[kept_intron_cols],
        'end': final_intron_ends,
        'feature': 'intron',
        'feature_number': intron_nums[kept_intron_cols]
    })

    all_features = pd.concat([df_exons, df_introns], ignore_index=True)
    if not all_features.empty:
        all_features = all_features.sort_values(['molecule_index', 'start'])
    
    return all_features

def estimate_half_lives_from_simulation(time_sums_df, spliced_df, min_points=3, n_bins=15):
    estimated_results = []
    
    for col in time_sums_df.columns:
        times = time_sums_df[col].values
        is_spliced = spliced_df[col].values
        
        mask = ~np.isnan(times)
        times = times[mask]
        is_spliced = is_spliced[mask]
        
        if len(times) == 0:
            continue

        try:
            if np.max(times) == 0:
                bins = [0, 0]
            else:
                bins = np.linspace(0, np.max(times), n_bins + 1)
            
            bin_times = []
            log_frac_unspliced = []
            
            bin_indices = np.digitize(times, bins)
            
            for i in range(1, len(bins)):
                bin_mask = (bin_indices == i)
                if np.sum(bin_mask) < 5: 
                    continue
                
                mean_time = np.mean(times[bin_mask])
                total = np.sum(bin_mask)
                spliced_count = np.sum(is_spliced[bin_mask])
                unspliced_count = total - spliced_count
                
                frac_unspliced = unspliced_count / total
                
                if frac_unspliced > 0:
                    bin_times.append(mean_time)
                    log_frac_unspliced.append(np.log(frac_unspliced))
            
            if len(bin_times) >= min_points:
                slope, intercept, r_value, p_value, std_err = linregress(bin_times, log_frac_unspliced)
                if slope < 0:
                    est_half_life = np.log(2) / -slope
                else:
                    est_half_life = np.inf 
                
                estimated_results.append({
                    'intron': col,
                    'estimated_half_life': est_half_life,
                    'r_squared': r_value**2
                })
            else:
                estimated_results.append({
                    'intron': col,
                    'estimated_half_life': np.nan,
                    'r_squared': np.nan
                })

        except Exception as e:
            estimated_results.append({
                'intron': col,
                'estimated_half_life': np.nan,
                'r_squared': np.nan
            })

    return pd.DataFrame(estimated_results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mRNA_df', help="full path to the mRNA dataframe")
    parser.add_argument('--intron_half_life', default='0.8,1.2', help="range of intron half life in minutes for the introns to simulate splicing")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")
    parser.add_argument("--nosplicing", action="store_true", help="If set, skip all splicing simulation steps and just output the mRNA_df")
    parser.add_argument("--mRNA_coordinates", action="store_true", help="If set, output a TSV with the coordinates of the surviving features (exons + retained introns)")

    args = parser.parse_args()
    
    filename = os.path.splitext(os.path.basename(args.mRNA_df))[0]
    base_filename = filename.split(".")[0]
    
    unique_seed = None
    if args.seed is not None:
        file_hash = int(hashlib.md5(filename.encode('utf-8')).hexdigest(), 16) % (10**8)
        unique_seed = args.seed + file_hash
        np.random.seed(unique_seed) 

    path_to_gtf = os.path.join(args.o, 'gtf', base_filename + ".tsv.gz")
    path_to_rates = os.path.join(args.o, 'rate_per_gene', base_filename + "_RatesandTraversalTimes.gtf")

    mRNA_df = pd.read_csv(args.mRNA_df, sep="\t", comment="#")
    output_filename_spliced_mRNAs = os.path.join(args.o, 'mRNA', base_filename + ".tsv.gz")

   
    if args.nosplicing or "_background" in filename or mRNA_df.empty:
        if not mRNA_df.empty:
            if "_background" in filename:
                mRNA_gtf = pd.read_csv(path_to_gtf, sep="\t", comment="#")
                if args.nosplicing:
                    bg_len = mRNA_gtf['sequence'].str.len().sum()
                else:
                    bg_len = mRNA_gtf[mRNA_gtf['feature'] == 'exon']['sequence'].str.len().sum()
                mRNA_df['sequence_length'] = bg_len
                mRNA_df['spliced_introns'] = "[]"
            else:
                mRNA_df['sequence_length'] = mRNA_df['stop_label_pos']
                mRNA_df['spliced_introns'] = "[]"
        mRNA_df.to_csv(output_filename_spliced_mRNAs, sep='\t', index=False)
        exit(0)

    mRNA_gtf = pd.read_csv(path_to_gtf, sep="\t", comment="#")
    mRNA_rates = pd.read_csv(path_to_rates, sep="\t", comment="#")

    strand = mRNA_gtf['strand'].unique()
    strand = strand[0] if len(strand) == 1 else strand[0]
    
    intron_df = mRNA_gtf[mRNA_gtf['feature'] == 'intron'].copy()

    if intron_df.shape[0] == 0:
        mRNA_df['sequence_length'] = mRNA_df['stop_label_pos']
        mRNA_df['spliced_introns'] = "[]"
        mRNA_df.to_csv(output_filename_spliced_mRNAs, sep='\t', index=False)
        exit(0)

    exon_row = mRNA_gtf[(mRNA_gtf['feature'] == 'exon') & (mRNA_gtf['position'] == 1)]
    if not exon_row.empty:
        TSS_coord = int(exon_row['start'].iloc[0]) if strand == '+' else int(exon_row['end'].iloc[0])
    else:
        sorted_gtf = mRNA_gtf.sort_values('start')
        if strand == '+':
            TSS_coord = int(sorted_gtf['start'].iloc[0])
        else:
            TSS_coord = int(sorted_gtf['end'].iloc[-1])
    
    if strand == '+':
        intron_df['start'] = intron_df['start'] - TSS_coord
        intron_df['end'] = intron_df['end'] - TSS_coord
    else:
        intron_df['start'] = TSS_coord - intron_df['start']
        intron_df['end'] = TSS_coord - intron_df['end']
    
    intron_df = intron_df[['start', 'end', 'position']]

    time_sums = calculate_time_sums_optimized(mRNA_df, intron_df, mRNA_rates, strand)
    time_sums_df = pd.DataFrame(time_sums, columns=[f'intron_{i+1}' for i in range(len(intron_df))])
    
    half_life_range = tuple(float(v) for v in args.intron_half_life.split(","))
    mean_half_life_dict = None

    output_filename_halflives = os.path.join(args.o, 'intron_half_lives', base_filename + ".tsv.gz")
    if os.path.exists(output_filename_halflives):
        mean_half_life_dict = pd.read_csv(output_filename_halflives, sep='\t', compression='gzip', index_col=0)['half_life'].to_dict()

    spliced_df, half_lives_df, run_mean_half_life_dict = run_splicing_per_exon_vary_introns(
        time_sums_df, half_life_range, mean_half_life_dict=mean_half_life_dict, seed=unique_seed)

    print(f"\n--- Splicing Validation: Estimated vs Input Half-Lives ({base_filename}) ---")
    validation_df = estimate_half_lives_from_simulation(time_sums_df, spliced_df)
    
    if not validation_df.empty:
        validation_df['input_half_life'] = validation_df['intron'].map(run_mean_half_life_dict)
        validation_df['intron_num'] = validation_df['intron'].str.extract(r'(\d+)').astype(float)
        validation_df = validation_df.sort_values('intron_num')
        
        for _, row in validation_df.iterrows():
            intron_name = row['intron']
            input_hl = row['input_half_life']
            est_hl = row['estimated_half_life']
            r2 = row['r_squared']
            
            if pd.isna(est_hl):
                print(f"{intron_name}: Not enough points for regression.")
            else:
                print(f"{intron_name}: Input={input_hl:.2f} min | Estimated={est_hl:.2f} min | R2={r2:.2f}")

        output_validation = os.path.join(args.o, 'intron_half_lives', base_filename + "_validation.tsv")
        validation_df.to_csv(output_validation, sep='\t', index=False)
    else:
        print("No validation data available (possibly no introns or no data points).")
    
    if args.mRNA_coordinates:
        print("Generating surviving feature coordinates...", flush=True)
        surviving_features_df = get_surviving_features(mRNA_df, spliced_df, mRNA_gtf, strand, TSS_coord)
        output_filename_coords = os.path.join(args.o, 'intron_half_lives', base_filename + "_mRNAfeature_coords.tsv.gz")
        surviving_features_df.to_csv(output_filename_coords, sep='\t', index=False, compression='gzip')
    
   
    mRNA_df_spliced = record_splicing_events(mRNA_df, spliced_df, intron_df, strand)
    mRNA_df_spliced.to_csv(output_filename_spliced_mRNAs, sep='\t', index=False)

    if mean_half_life_dict is None:
        half_life_df = pd.DataFrame(run_mean_half_life_dict.items(), columns=['intron', 'half_life'])
        half_life_df.to_csv(output_filename_halflives, sep='\t', compression='gzip', index=False)

        intron_df_out = mRNA_gtf[mRNA_gtf['feature'] == 'intron'].copy()
        half_life_map = {int(k.replace('intron_', '')): v for k, v in run_mean_half_life_dict.items()}
        intron_df_out['half_life'] = intron_df_out['position'].map(half_life_map)
        output_intron_abscoords = os.path.join(args.o, 'intron_half_lives', base_filename + "_abscoords.tsv.gz")
        if 'sequence' in intron_df_out.columns:
            intron_df_to_save = intron_df_out.drop(columns=['sequence'])
        else:
            intron_df_to_save = intron_df_out
        intron_df_to_save.to_csv(output_intron_abscoords, sep='\t', compression='gzip', index=False)
