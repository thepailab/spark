import pandas as pd
import argparse
import math
import os
import string
import numpy as np
import zlib
import hashlib  

# 1: Create the DNA sequence as a numpy array for fast slicing/mutation
def getDNAseq_array(dna_sequence):
    # Convert to numpy array of characters for fast coordinate finding
    return np.array(list(''.join(str(i) for i in dna_sequence)), dtype='U1')

# 4: Optimized findStopnt
def findStopnt_optimized(times_array, cumsum_times_padded, labeling_duration, startsite, rng):
    idx_start = max(0, startsite - 1)
    
    start_time_mean = cumsum_times_padded[idx_start]
    target_time_mean = start_time_mean + labeling_duration
    
    idx_approx = np.searchsorted(cumsum_times_padded, target_time_mean, side='left')
    
    needed = idx_approx - idx_start
    idx_end_safe = min(len(times_array), idx_start + int(needed * 1.2) + 50)
    
    times_slice = times_array[idx_start : idx_end_safe]
    jitters = rng.uniform(0.9, 1.1, size=times_slice.size)
    total_times = np.cumsum(times_slice * jitters)
    
    idx_in_slice = np.searchsorted(total_times, labeling_duration, side="right")
    
    return startsite + idx_in_slice

# 5 & 6: Vectorized Coordinate Extractor (No String Manipulation)
def generate_mutated_read(seq_array_ref, read_stop_idx, label_start_idx, 
                          nt_inc_range, sub_rate_range, seq_err_range, 
                          from_base, to_base, rng):
    """
    Finds the coordinates for labeling and sequencing errors without 
    materializing or mutating the string sequence.
    """
    stop_safe = min(len(seq_array_ref), read_stop_idx)
    lbl_start_safe = max(0, label_start_idx)
    
    nt_inc_rate_val = 0.0
    incorporated_absolute = []
    converted_absolute = []
    
    # 1. LABELING (Incorporation & Conversion)
    if lbl_start_safe < stop_safe:
        # Look at the segment transcribed during the pulse
        label_view = seq_array_ref[lbl_start_safe:stop_safe]
        
        # Find local indices of the analog base (e.g., 'T')
        from_indices = np.where(label_view == from_base)[0]
        n_candidates = len(from_indices)
        
        nt_inc_rate_val = rng.uniform(nt_inc_range[0], nt_inc_range[1])
        n_inc = int(np.ceil(n_candidates * nt_inc_rate_val))
        
        if n_inc > 0:
            chosen_indices_rel = rng.choice(from_indices, size=n_inc, replace=False)
            incorporated_absolute = (lbl_start_safe + chosen_indices_rel).tolist()
            
            # Conversion (e.g., T to C)
            if to_base and to_base != "":
                conv_frac = rng.uniform(sub_rate_range[0], sub_rate_range[1])
                n_conv = int(np.floor(n_inc * conv_frac))
                
                if n_conv > 0:
                    conv_indices_rel = rng.choice(chosen_indices_rel, size=n_conv, replace=False)
                    converted_absolute = (lbl_start_safe + conv_indices_rel).tolist()

    # 2. SEQUENCING ERRORS (Determine random coordinates across the whole read)
    read_len = stop_safe
    seq_err_rate = rng.uniform(seq_err_range[0], seq_err_range[1])
    n_err = rng.binomial(read_len, seq_err_rate)
    
    seq_err_absolute = []
    if n_err > 0:
        err_locs = rng.choice(read_len, size=n_err, replace=False)
        seq_err_absolute = [int(x) for x in err_locs]
    
    return nt_inc_rate_val, converted_absolute, incorporated_absolute, seq_err_rate, seq_err_absolute

# 8: Bulk ID Generator
def generate_ids_bulk(n, size=12, rng=None):
    chars = list(string.ascii_uppercase + string.digits)
    rand_chars = rng.choice(chars, (n, size-1))
    return ["@" + "".join(row) for row in rand_chars]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--region_file', default='tmp.gtf', help='Input region specific elongation rate GTF file name')
    parser.add_argument('--nt_file', default='tmp.gtf', help='Input ntTraversalTime GTF file name')
    parser.add_argument('--experiment_time', type=int, default=15, help='Labeling time in minutes')
    parser.add_argument('--mRNAcount', type=int, default=5000, help='Number of mRNA molecules to simulate')  
    parser.add_argument('--seq_err', default='0.0001,0.0002', help='Illumina sequencing error rate range')
    parser.add_argument('--nt_inc_rate',type=str,default='0.09,0.1',help='Comma-separated range of nucleotide incorporation proportions')
    parser.add_argument('--subs_rate',type=str,default='0.95,1',help='Comma-separated range of nucleotide substitution proportions')
    parser.add_argument('--labeling_base', type=str, default='T', help='Specify the base analog that will get incorporated into the nascentRNA', required=False)
    parser.add_argument('--sub_base', type=str, help='Specify the identity of the base after conversion', required=False)
    parser.add_argument('--initi_rate',type=str,default='10,60',help='Comma-separated range of seconds for initiation rate')
    parser.add_argument('--bkg_molecules', type=float,default=0,help='Proportion of molecules that are background (0 or 1)')
    parser.add_argument('--path_to_tpm', type=str,help='Full path to the TPM per gene file')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument("--drb", action="store_true", help="DRB treatment experiment for transcription synchronization")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")
    parser.add_argument("--nosplicing", action="store_true", help="If set, skip all splicing simulation steps")

    args = parser.parse_args()

    df = pd.read_csv(args.region_file, sep='\t', comment='#')
    df.columns = ["chromosome", 'absolute_start', 'absolute_end', "region_start_coord", "region_end_coord", 'strand', "rate_initial", "rate_final", "rate_change_per_nt", "time_to_traverse", "sequence"]

    if args.seed is not None:
        gene_signature = f"{df.iloc[0]['chromosome']}_{df.iloc[0]['absolute_start']}_{df.iloc[0]['absolute_end']}"
        hasher = hashlib.sha256(gene_signature.encode('utf-8'))
        gene_content_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_content_hash) % (2**32)
        np.random.seed(unique_seed)
        rng = np.random.default_rng(unique_seed)
    else:
        rng = np.random.default_rng()

    nt_incorporation_rate = list(map(float, args.nt_inc_rate.split(',')))
    sub_rate_percent_range = list(map(float, args.subs_rate.split(',')))
    seqerr_range = list(map(float, args.seq_err.split(',')))
    
    df_nt = pd.read_csv(args.nt_file, sep='\t', comment='#')
    df_nt.columns = ["chromosome", "absolute_position", 'strand', "region_number", "nucleotide_coord", "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    
    dna_sequence_arr = getDNAseq_array(df.sequence)
    times_raw = df_nt['time_for_this_nt'].values
    cumsum_times_padded = np.concatenate(([0], np.cumsum(times_raw)))

    filename = os.path.splitext(os.path.basename(args.region_file))[0]
    base_filename = filename.split("_")[0]
    output_filename = os.path.join(args.o, base_filename + ".tsv.gz")

    data_list = []
    from_base = args.labeling_base
    to_base = args.sub_base
    max_transcript_length = df["region_end_coord"].max()

    low, high = map(int, args.initi_rate.split(","))
    initiation_interval_sec = rng.integers(low, high + 1)
    initiation_rate_per_sec = 1.0 / initiation_interval_sec
    
    input_path = args.region_file
    pausing_file_path = input_path.replace('_VariableElongationRateRegions.gtf', '_total_pausing_time.tsv')
    total_pausing_times = 0.0
    try:
        pausing_df = pd.read_csv(pausing_file_path, sep='\t')
        if not pausing_df.empty:
            total_pausing_times = pausing_df['total_pause_time_minutes'].iloc[0]
    except:
        pass 

    mean_elong_rate = df_nt['rate_for_this_nt'].mean()
    n_pol_per_unit = (max_transcript_length / mean_elong_rate + total_pausing_times) / (initiation_interval_sec / 60)
    
    bolus_per_unit = 0
    if args.drb:
        n_pol_per_unit = 0 
        bolus_per_unit = 1 

    target_sample_size = int(args.mRNAcount * args.experiment_time / 10)
    needed_molecules = int(target_sample_size * 1.2)
    if args.bkg_molecules == 1:
        needed_molecules = 0
    elif args.bkg_molecules > 0:
        needed_molecules = int(needed_molecules / (1 - args.bkg_molecules))
    
    initiations_per_unit = args.experiment_time / (initiation_interval_sec / 60)
    total_molecules_per_unit = n_pol_per_unit + bolus_per_unit + initiations_per_unit

    if total_molecules_per_unit > 0:
        scaling_factor = math.ceil(needed_molecules / total_molecules_per_unit)
    else:
        scaling_factor = 1 

    print(f"Optimizing run: Target={target_sample_size}, Scaler set to {scaling_factor} (was 50)")
    
    already_initiated = int(n_pol_per_unit * scaling_factor)
    n_mol_to_initiate = scaling_factor * args.experiment_time / (initiation_interval_sec / 60)

    single_nt_mask = df['region_start_coord'] == df['region_end_coord']
    inverse_weights = df['time_to_traverse']
    total_inverse_weight = inverse_weights.sum()
    expected_counts = (inverse_weights * already_initiated / total_inverse_weight).round().astype(int)
    
    single_nt_sample_counts = expected_counts[single_nt_mask].clip(upper=int(scaling_factor))
    single_nt_samples = df.loc[single_nt_sample_counts.index].loc[single_nt_sample_counts.index.repeat(single_nt_sample_counts)]

    n_single_nt = single_nt_sample_counts.sum()
    other_regions = df[~single_nt_mask]
    n_to_sample = already_initiated - n_single_nt

    if n_to_sample > 0:
        weights_other = inverse_weights[~single_nt_mask].values
        weights_other = weights_other / weights_other.sum()
        chosen_indices = rng.choice(len(other_regions), size=n_to_sample, p=weights_other, replace=True)
        other_samples = other_regions.iloc[chosen_indices]
    else:
        other_samples = pd.DataFrame()

    chosen_samples = pd.concat([single_nt_samples, other_samples])

    if not chosen_samples.empty:
        random_positions = rng.integers(low=chosen_samples['region_start_coord'].values, high=chosen_samples['region_end_coord'].values + 1)
        start_label_pos_list = random_positions.tolist()
    else:
        start_label_pos_list = []

    initiation_times = []
    experiment_time_sec = args.experiment_time * 60
    max_value = int(total_pausing_times * 60)
    
    if len(single_nt_samples) > 0:
        initiation_times.extend(experiment_time_sec + rng.integers(0, max_value + 1, size=len(single_nt_samples)))
    if len(other_samples) > 0:
        initiation_times.extend([experiment_time_sec] * len(other_samples))

    num_cells_needed = scaling_factor
    all_new_initiations = []

    if num_cells_needed > 0:
        for _ in range(int(num_cells_needed)):
            if args.drb:
                all_new_initiations.append(experiment_time_sec)
            if n_mol_to_initiate > 0:
                current_time = 0.0
                while current_time < experiment_time_sec:
                    wait = rng.exponential(1.0 / initiation_rate_per_sec)
                    current_time += wait
                    if current_time < experiment_time_sec:
                        labeling_duration_sec = experiment_time_sec - current_time
                        all_new_initiations.append(labeling_duration_sec)

    initiation_times.extend(all_new_initiations)
    start_label_pos_list.extend([1] * len(all_new_initiations))

    total_molecules = len(initiation_times)
    molecule_ids = generate_ids_bulk(total_molecules, size=12, rng=rng)

    # --- MAIN LOOP ---
    for i in range(total_molecules):
        initiation = initiation_times[i] / 60 
        start_label_pos = int(start_label_pos_list[i])
        stop_label_pos = findStopnt_optimized(times_raw, cumsum_times_padded, initiation, start_label_pos, rng)
        
        # Call the new refactored function
        pct_sub, conv_pos, inc_pos, pct_err, err_pos = generate_mutated_read(
            dna_sequence_arr,
            int(stop_label_pos),
            int(start_label_pos),
            nt_incorporation_rate,
            sub_rate_percent_range,
            seqerr_range,
            from_base,
            to_base,
            rng
        )
        
        if initiation > args.experiment_time:
            initiation = args.experiment_time

        row_data = [
            initiation,
            molecule_ids[i],
            df_nt.loc[0, 'strand'],
            sub_rate_percent_range,
            start_label_pos,
            stop_label_pos,
            conv_pos,
            inc_pos,
            pct_err,
            err_pos
        ] # full_molecule_sequence has been DROPPED
        data_list.append(row_data)

    # --- Background Molecules ---
    bkg_molecules = args.bkg_molecules
    if bkg_molecules > 0:
        bg_data_list = []
        output_filename_bg = os.path.join(args.o, base_filename + "_background.tsv.gz")
        num_to_sample = round(len(data_list) * (1 - bkg_molecules))
        
        indices = np.arange(len(data_list))
        rng.shuffle(indices)
        kept_indices = indices[:num_to_sample]
        data_list = [data_list[i] for i in kept_indices]
        
        bkg_molecules_to_simulate = 100

        # Retrieve lengths quickly without materializing background string loops
        idx = output_filename.find('pre-mRNA')
        up_to = output_filename[:idx]
        gtf_path_in = os.path.join(up_to, 'gtf', base_filename + '.tsv.gz')

        gtf_gene = pd.read_csv(gtf_path_in, sep='\t', comment='#')
        if args.nosplicing:
            bg_len = gtf_gene['sequence'].str.len().sum()
        else:
            exon_gtf = gtf_gene[gtf_gene['feature'] == 'exon']
            bg_len = exon_gtf['sequence'].str.len().sum()
        
        bg_ids = generate_ids_bulk(bkg_molecules_to_simulate, size=12, rng=rng)
        
        for i in range(bkg_molecules_to_simulate):
            seq_err_rate = rng.uniform(seqerr_range[0], seqerr_range[1])
            n_err = rng.binomial(bg_len, seq_err_rate)
            err_pos = []
            if n_err > 0:
                err_locs = rng.choice(bg_len, size=n_err, replace=False)
                err_pos = [int(x) for x in err_locs]
            
            bkg_row_data = [
                0,
                bg_ids[i],
                df_nt.loc[0, 'strand'],
                0,
                0,
                0,
                'NA',
                'NA',
                seq_err_rate,
                err_pos
            ] # full_molecule_sequence has been DROPPED
            bg_data_list.append(bkg_row_data)

        # Updated export columns
        bg_data_list_for_export = pd.DataFrame(bg_data_list, columns=[
            "initiation_time", 'molecule_id', 'strand', "sub_rate_percent_range", 
            'start_label_pos', 'stop_label_pos', 'converted_positions', 
            'incorporated_positions', 'percentage_seq_err', 'seq_err_positions'
        ])
        bg_data_list_for_export.to_csv(output_filename_bg, sep='\t', index=False, compression='gzip')

    df_for_export = pd.DataFrame(data_list, columns=[
        "initiation_time", 'molecule_id', 'strand', "sub_rate_percent_range", 
        'start_label_pos', 'stop_label_pos', 'converted_positions', 
        'incorporated_positions', 'percentage_seq_err', 'seq_err_positions'
    ])

    if len(df_for_export) > target_sample_size:
        df_for_export = df_for_export.sample(n=target_sample_size, random_state=42)

    df_for_export.to_csv(output_filename, sep='\t', index=False, compression='gzip')

    summary_path = os.path.join(args.o, "initiation_rates_all.tsv")
    df_summary = pd.DataFrame([{"gene_id": base_filename, "initiation_rate": initiation_interval_sec}])
    file_exists = os.path.exists(summary_path)
    df_summary.to_csv(summary_path, sep="\t", index=False, mode="a", header=not file_exists)
