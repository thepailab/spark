import math
import numpy as np
import pandas as pd
import argparse
import os
from scipy.stats import linregress
import hashlib


def atomic_to_csv(df, path, **kwargs):
    """Write a DataFrame to `path` atomically.

    Writes to a unique temp file in the same directory, then os.replace()s it
    into place. os.replace is atomic on the same filesystem, so a concurrent or
    re-queued writer can never leave a half-written file (or a valid gzip stream
    with trailing garbage appended) behind. If `path` ends in .gz, gzip is set
    explicitly rather than relying on pandas' extension inference.
    """
    if path.endswith(".gz") and "compression" not in kwargs:
        kwargs["compression"] = "gzip"
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        df.to_csv(tmp, **kwargs)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

def calculate_time_sums_optimized(mRNA_df, intron_df, mRNA_rates, strand):
    mRNA_rates_sorted = mRNA_rates.sort_values('nucleotide_coord').reset_index(drop=True)
    coords = mRNA_rates_sorted['nucleotide_coord'].values
    times  = mRNA_rates_sorted['time_for_this_nt'].values
    cum_times = np.concatenate(([0], np.cumsum(times)))

    initiation_times = mRNA_df['initiation_time'].values
    stop_positions   = mRNA_df['stop_label_pos'].values

    if strand == '+':
        splice_sites_3 = intron_df['end'].values
    else:
        splice_sites_3 = intron_df['start'].values

    ss3_idx = np.searchsorted(coords, splice_sites_3, side='right')
    ss3_idx = np.clip(ss3_idx, 0, len(cum_times) - 1)
    cum_at_ss3 = cum_times[ss3_idx]

    t_available = initiation_times[:, np.newaxis] - cum_at_ss3[np.newaxis, :]
    t_available = np.maximum(0.0, t_available)

    not_reached = stop_positions[:, np.newaxis] < splice_sites_3[np.newaxis, :]
    t_available = np.where(not_reached, np.nan, t_available)

    return t_available.tolist()

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
        valid_t = np.maximum(valid_t, 0.0)
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
    spliced_matrix = spliced_df.values
    intron_starts = intron_df['start'].values
    intron_ends = intron_df['end'].values
    stop_positions = mRNA_df['stop_label_pos'].values
    start_positions = mRNA_df['start_label_pos'].values # ADD THIS
    
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
                
                if strand == '+':
                    rs, re = start - 1, end
                else:
                    rs, re = end, start + 1
                    
                
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
    
    #extract molecule IDs
    mol_ids = mRNA_df['molecule_id'].values

    exon_mask = exon_starts[np.newaxis, :] <= mol_stops
    kept_exon_rows, kept_exon_cols = np.nonzero(exon_mask)
    
    current_mol_stops_ex = mol_stops[kept_exon_rows, 0]
    current_exon_ends = exon_ends[kept_exon_cols]
    final_exon_ends = np.minimum(current_exon_ends, current_mol_stops_ex)
    
    df_exons = pd.DataFrame({
        'molecule_id': mol_ids[kept_exon_rows],
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
        'molecule_id': mol_ids[kept_intron_rows], 
        'start': intron_starts[kept_intron_cols],
        'end': final_intron_ends,
        'feature': 'intron',
        'feature_number': intron_nums[kept_intron_cols]
    })

    all_features = pd.concat([df_exons, df_introns], ignore_index=True)
    if not all_features.empty:
        all_features = all_features.sort_values(['molecule_id', 'start'])
    
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
                
                correct_gtf_path = path_to_gtf.replace("_background.tsv.gz", ".tsv.gz")
                mRNA_gtf = pd.read_csv(correct_gtf_path, sep="\t", comment="#")
                
                if args.nosplicing:
                    bg_len = mRNA_gtf['sequence'].str.len().sum()
                else:
                    bg_len = mRNA_gtf[mRNA_gtf['feature'] == 'exon']['sequence'].str.len().sum()
                mRNA_df['sequence_length'] = bg_len
                mRNA_df['spliced_introns'] = "[]"
            else:
                mRNA_df['sequence_length'] = mRNA_df['stop_label_pos']
                mRNA_df['spliced_introns'] = "[]"
        atomic_to_csv(mRNA_df, output_filename_spliced_mRNAs, sep='\t', index=False)
        exit(0)

    mRNA_gtf = pd.read_csv(path_to_gtf, sep="\t", comment="#")
    mRNA_rates = pd.read_csv(path_to_rates, sep="\t", comment="#")

    strand = mRNA_gtf['strand'].unique()
    strand = strand[0] if len(strand) == 1 else strand[0]
    
    intron_df = mRNA_gtf[mRNA_gtf['feature'] == 'intron'].copy()

    if intron_df.shape[0] == 0:
        mRNA_df['sequence_length'] = mRNA_df['stop_label_pos']
        mRNA_df['spliced_introns'] = "[]"
        atomic_to_csv(mRNA_df, output_filename_spliced_mRNAs, sep='\t', index=False)
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
        atomic_to_csv(validation_df, output_validation, sep='\t', index=False)
    else:
        print("No validation data available (possibly no introns or no data points).")
    
    if args.mRNA_coordinates:
        print("Generating surviving feature coordinates...", flush=True)
        surviving_features_df = get_surviving_features(mRNA_df, spliced_df, mRNA_gtf, strand, TSS_coord)
        output_filename_coords = os.path.join(args.o, 'intron_half_lives', base_filename + "_mRNAfeature_coords.tsv.gz")
        atomic_to_csv(surviving_features_df, output_filename_coords, sep='\t', index=False)
    
   
    mRNA_df_spliced = record_splicing_events(mRNA_df, spliced_df, intron_df, strand)
    atomic_to_csv(mRNA_df_spliced, output_filename_spliced_mRNAs, sep='\t', index=False)

    if mean_half_life_dict is None:
        half_life_df = pd.DataFrame(run_mean_half_life_dict.items(), columns=['intron', 'half_life'])
        atomic_to_csv(half_life_df, output_filename_halflives, sep='\t', index=False)

        intron_df_out = mRNA_gtf[mRNA_gtf['feature'] == 'intron'].copy()
        
        #create a sequential 1-based index for the introns to map against
        intron_df_out['intron_idx'] = range(1, len(intron_df_out) + 1)
        half_life_map = {int(k.replace('intron_', '')): v for k, v in run_mean_half_life_dict.items()}
        
        #nmap using the new index, NOT the GTF position
        intron_df_out['half_life'] = intron_df_out['intron_idx'].map(half_life_map)
        
        output_intron_abscoords = os.path.join(args.o, 'intron_half_lives', base_filename + "_abscoords.tsv.gz")
        
        cols_to_drop = ['intron_idx']
        if 'sequence' in intron_df_out.columns:
            cols_to_drop.append('sequence')
            
        intron_df_to_save = intron_df_out.drop(columns=cols_to_drop)
        atomic_to_csv(intron_df_to_save, output_intron_abscoords, sep='\t', index=False)
