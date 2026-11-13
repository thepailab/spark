import pandas as pd
import argparse
import math
import os
import string
import numpy as np
import hashlib
import pysam
from spark_labeling import parse_conversion_map


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

def getDNAseq_array(dna_sequence):
    return np.array(list(''.join(str(i) for i in dna_sequence)), dtype='U1')

def reverse_complement(seq):
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]

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

def generate_mutated_read(seq_array_ref, read_stop_idx, label_start_idx,
                          nt_inc_range, sub_rate_range, seq_err_range,
                          from_bases, to_base, rng, label=True):
    stop_safe = min(len(seq_array_ref), read_stop_idx)
    lbl_start_safe = max(0, label_start_idx)
    
    nt_inc_rate_val = 0.0
    incorporated_absolute = []
    converted_absolute = []
    
    if label and lbl_start_safe < stop_safe:
        label_view = seq_array_ref[lbl_start_safe:stop_safe]
        from_indices = np.where(np.isin(label_view, from_bases))[0]
        n_candidates = len(from_indices)
        
        nt_inc_rate_val = rng.uniform(nt_inc_range[0], nt_inc_range[1])
        n_inc = int(np.ceil(n_candidates * nt_inc_rate_val))
        
        if n_inc > 0:
            chosen_indices_rel = rng.choice(from_indices, size=n_inc, replace=False)
            incorporated_absolute = (lbl_start_safe + chosen_indices_rel).tolist()
            
            if to_base and to_base != "":
                conv_frac = rng.uniform(sub_rate_range[0], sub_rate_range[1])
                n_conv = int(np.floor(n_inc * conv_frac))
                
                if n_conv > 0:
                    conv_indices_rel = rng.choice(chosen_indices_rel, size=n_conv, replace=False)
                    converted_absolute = (lbl_start_safe + conv_indices_rel).tolist()

    read_len = stop_safe - lbl_start_safe
    if read_len <= 0:
        return nt_inc_rate_val, converted_absolute, incorporated_absolute, 0.0, []

    seq_err_rate = rng.uniform(seq_err_range[0], seq_err_range[1])
    n_err = rng.binomial(read_len, seq_err_rate)
    
    seq_err_absolute = []
    if n_err > 0:
        err_locs = rng.choice(read_len, size=n_err, replace=False)
        seq_err_absolute = [int(x + lbl_start_safe) for x in err_locs]
    
    return nt_inc_rate_val, converted_absolute, incorporated_absolute, seq_err_rate, seq_err_absolute

def generate_ids_bulk(n, size=12, rng=None):
    chars = list(string.ascii_uppercase + string.digits)
    rand_chars = rng.choice(chars, (n, size-1))
    return ["@" + "".join(row) for row in rand_chars]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--region_file', default='tmp.gtf')
    parser.add_argument('--nt_file', default='tmp.gtf')
    parser.add_argument('--experiment_time', type=int, default=15)
    parser.add_argument('--experiment_type', choices=['metaboliclabeling', 'ttseq', 'netseq', 'proseq','chromatin_associated'], default='metaboliclabeling')
    parser.add_argument('--mRNAcount', type=int, default=5000)  
    parser.add_argument('--seq_err', default='0.0001,0.0002')
    parser.add_argument('--nt_inc_rate',type=str,default='0.09,0.1')
    parser.add_argument('--subs_rate',type=str,default='0.95,1')
    parser.add_argument('--labeling_base', type=str, default='T')
    parser.add_argument('--sub_base', type=str)
    parser.add_argument('--initi_rate',type=str,default='10,60')
    parser.add_argument('--bkg_molecules', type=float,default=0)
    parser.add_argument('--path_to_tpm', type=str)
    parser.add_argument('--o', type=str, default='./')
    parser.add_argument("--drb", action="store_true")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--nosplicing", action="store_true")
    parser.add_argument("--cleavage_half_life", type=str, default="0.0,0.0")
    parser.add_argument("--pas_elong_rate", type=str, default="0.5,2.0")
    parser.add_argument("--readthrough_deg_rate", type=float, default=1.56)
    parser.add_argument("--nocleavage", action="store_true")
    parser.add_argument("--immediate_cleavage", action="store_true")
    parser.add_argument("--polyA_length", type=str, default="50,150")
    parser.add_argument("--genome_fasta", type=str)
    parser.add_argument("--reuse_parameters", action="store_true")

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
    polyA_len_min, polyA_len_max = map(int, args.polyA_length.split(','))
    
    df_nt = pd.read_csv(args.nt_file, sep='\t', comment='#')
    df_nt.columns = ["chromosome", "absolute_position", 'strand', "region_number", "nucleotide_coord", "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    
    dna_sequence_arr = getDNAseq_array(df.sequence)
    
    times_raw = df_nt['time_for_this_nt'].values * 60.0
    cumsum_times_padded = np.concatenate(([0], np.cumsum(times_raw)))

    filename = os.path.basename(args.region_file)
    base_filename = filename.replace("_VariableElongationRateRegions.gtf", "")
    output_filename = os.path.join(args.o, base_filename + ".tsv.gz")

    max_transcript_length = df["region_end_coord"].max()
    
    abs_chrom = str(df.iloc[0]['chromosome'])
    abs_strand = str(df.iloc[0]['strand'])
    abs_end_pos = int(df['absolute_end'].max())
    abs_start_pos = int(df['absolute_start'].min())

    idx = output_filename.find('pre-mRNA')
    up_to = output_filename[:idx] if idx != -1 else os.path.dirname(args.o)
    gtf_path_in = os.path.join(up_to, 'gtf', base_filename + '.tsv.gz')
    
    try:
        gtf_gene = pd.read_csv(gtf_path_in, sep='\t', comment='#')
        gtf_gene = gtf_gene[gtf_gene['feature'] != 'readthrough']
        if 'pas_coordinate' in gtf_gene.columns:
            PAS_coord = int(gtf_gene['pas_coordinate'].iloc[0])
        else:
            PAS_coord = max_transcript_length
    except:
        gtf_gene = pd.DataFrame()
        PAS_coord = max_transcript_length

    data_list = []
    from_bases, conv_map = parse_conversion_map(args.labeling_base, args.sub_base)
    to_base = args.sub_base if conv_map else None

    parameters_path = os.path.join(os.path.dirname(args.region_file), f"{base_filename}_parameters.tsv")
    if args.reuse_parameters and os.path.exists(parameters_path):
        stored_parameters = pd.read_csv(parameters_path, sep='\t').iloc[0]
        initiation_interval_sec = float(stored_parameters['initiation_interval_sec'])
        gene_pas_rate = float(stored_parameters['pas_elongation_rate_nt_per_sec'])
        gene_cleavage_half_life = float(stored_parameters['cleavage_half_life_sec'])
    else:
        low, high = map(float, args.initi_rate.split(","))
        initiation_interval_sec = float(np.exp(rng.uniform(np.log(low), np.log(high))))

        pas_rate_kb_min = list(map(float, args.pas_elong_rate.split(',')))
        gene_pas_rate = rng.uniform(pas_rate_kb_min[0] * 1000 / 60, pas_rate_kb_min[1] * 1000 / 60)

        cleavage_hl_range = list(map(float, args.cleavage_half_life.split(',')))
        if cleavage_hl_range[1] > 0:
            hl_min, hl_max = cleavage_hl_range[0], cleavage_hl_range[1]
            sampled_hl_min = rng.uniform(hl_min, hl_max)
            gene_cleavage_half_life = sampled_hl_min * 60.0
        else:
            gene_cleavage_half_life = 0.0

        atomic_to_csv(pd.DataFrame([{
            "gene_id": base_filename,
            "initiation_interval_sec": initiation_interval_sec,
            "pas_elongation_rate_nt_per_sec": gene_pas_rate,
            "cleavage_half_life_sec": gene_cleavage_half_life,
        }]), parameters_path, sep="\t", index=False)

    initiation_rate_per_sec = 1.0 / initiation_interval_sec
    base_readthrough_deg_rate_nt_sec = (args.readthrough_deg_rate * 1000) / 60.0

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
    n_pol_per_unit = cumsum_times_padded[-1] / initiation_interval_sec
    
    bolus_per_unit = 0
    if args.drb:
        n_pol_per_unit = 0 
        bolus_per_unit = 1 

    target_sample_size = int(args.mRNAcount * args.experiment_time / 10)

    if args.experiment_type in ['netseq', 'proseq']:
        target_sample_size = target_sample_size * 10

    bkg_molecules_to_simulate = int(target_sample_size * args.bkg_molecules)
    nascent_to_simulate = target_sample_size - bkg_molecules_to_simulate

    if args.bkg_molecules == 1.0:
        needed_molecules = 0
    else:
        needed_molecules = int(nascent_to_simulate * 1.2)

    initiations_per_unit = args.experiment_time / (initiation_interval_sec / 60)
    total_molecules_per_unit = n_pol_per_unit + bolus_per_unit + initiations_per_unit

    if args.experiment_type in ['netseq', 'proseq']:
        effective_molecules_per_unit = n_pol_per_unit + bolus_per_unit
    else:
        effective_molecules_per_unit = total_molecules_per_unit

    if effective_molecules_per_unit > 0:
        scaling_factor = math.ceil(needed_molecules / effective_molecules_per_unit)
    else:
        scaling_factor = math.ceil(needed_molecules / total_molecules_per_unit) if total_molecules_per_unit > 0 else 1

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
    
    n_single = len(single_nt_samples)
    
    if n_single > 0:
        single_positions = start_label_pos_list[:n_single]
        times_to_reach_single = cumsum_times_padded[single_positions]
        initiation_times.extend(experiment_time_sec + times_to_reach_single + rng.integers(0, max_value + 1, size=n_single))
        
    if len(other_samples) > 0:
        other_positions = start_label_pos_list[n_single:]
        times_to_reach_other = cumsum_times_padded[other_positions]
        initiation_times.extend(experiment_time_sec + times_to_reach_other)

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
                        all_new_initiations.append(labeling_duration_sec + cumsum_times_padded[1])

    initiation_times.extend(all_new_initiations)
    start_label_pos_list.extend([1] * len(all_new_initiations))

    total_molecules = len(initiation_times)
    molecule_ids = generate_ids_bulk(total_molecules, size=12, rng=rng)

    time_to_PAS = cumsum_times_padded[PAS_coord] if PAS_coord < len(cumsum_times_padded) else cumsum_times_padded[-1]
    time_to_end = cumsum_times_padded[-1]

    def emit_molecule(physical_start, physical_stop, label_start, mol_id_suffix, polyA_len, init_sec, m_id, time_since_cleaved_sec):
        if physical_start >= physical_stop: return
        eff_label_start = max(label_start, physical_start)
        
        pct_sub, conv_pos, inc_pos, pct_err, err_pos = generate_mutated_read(
            dna_sequence_arr, int(physical_stop), int(eff_label_start),
            nt_incorporation_rate, sub_rate_percent_range, seqerr_range,
            from_bases, to_base, rng, label=args.experiment_type != 'netseq'
        )
        
        spliced_introns_str = f"[[0, {int(physical_start)}]]" if physical_start > 0 else "[]"
        
        row = [
            init_sec / 60.0,
            f"{m_id}{mol_id_suffix}",
            df_nt.loc[0, 'strand'],
            sub_rate_percent_range,
            int(physical_start),  #Pol II position when teh exp starts
            int(physical_stop),
            conv_pos,
            inc_pos,
            pct_err,
            err_pos,
            spliced_introns_str,
            polyA_len,
            time_since_cleaved_sec / 60.0
        ]
        data_list.append(row)

    molecules_to_emit = []
    max_physical_stop = PAS_coord

    for i in range(total_molecules):
        initiation_sec = initiation_times[i]
        start_label_pos = int(start_label_pos_list[i])
        
        time_spent_at_start = cumsum_times_padded[start_label_pos] if start_label_pos < len(cumsum_times_padded) else cumsum_times_padded[-1]
        duration_transcribing = max(0.0, initiation_sec - time_spent_at_start)
        
        if initiation_sec <= time_to_PAS:
            stop_label_pos = findStopnt_optimized(times_raw, cumsum_times_padded, duration_transcribing, start_label_pos, rng)
            if args.immediate_cleavage:
                stop_label_pos = min(stop_label_pos, int(PAS_coord))
        else:
            if args.immediate_cleavage:
                stop_label_pos = int(PAS_coord)
            else:
                time_past_PAS = initiation_sec - time_to_PAS
                mol_pas_rate = gene_pas_rate * rng.uniform(0.9, 1.1)
                dist_past = int(time_past_PAS * mol_pas_rate)
                stop_label_pos = int(PAS_coord + dist_past) 
            
        cleaved = False
        deg_front = PAS_coord
        time_since_cleavage_sec = 0.0 #get the time after cleavage occured
        
        if not args.nocleavage:
            if args.immediate_cleavage:
                if stop_label_pos >= PAS_coord:
                    cleaved = True
                    deg_front = PAS_coord
                    time_since_cleavage_sec = initiation_sec - time_to_PAS 
            elif stop_label_pos > PAS_coord:
                time_past_PAS = initiation_sec - time_to_PAS 
                
                if time_past_PAS > 0:
                    if gene_cleavage_half_life > 0:
                        mean_life = gene_cleavage_half_life / np.log(2)
                        u = rng.random()
                        t_cleave = -mean_life * np.log(u) if u > 0 else float('inf')
                        
                        if t_cleave < time_past_PAS:
                            cleaved = True
                            time_since_cleavage_sec = time_past_PAS - t_cleave 
                            
                            t_deg = time_since_cleavage_sec 
                            mol_deg_rate = base_readthrough_deg_rate_nt_sec * rng.uniform(0.9, 1.1)
                            degraded_dist = int(t_deg * mol_deg_rate)
                            deg_front = min(PAS_coord + degraded_dist, stop_label_pos)

        if start_label_pos >= stop_label_pos:
            continue
            
        max_physical_stop = max(max_physical_stop, stop_label_pos)
        
        molecules_to_emit.append({
            'initiation_sec': initiation_sec,
            'start_label_pos': start_label_pos,
            'stop_label_pos': stop_label_pos,
            'cleaved': cleaved,
            'deg_front': deg_front,
            'molecule_id': molecule_ids[i],
            'time_since_cleavage_sec': time_since_cleavage_sec 
        })

    extra_needed = max_physical_stop - len(dna_sequence_arr)
    
    if extra_needed > 0:
        buffer_pad = 150
        extra_needed += buffer_pad
        appended_seq = ""
        fetch_start = abs_end_pos
        fetch_end = abs_end_pos + extra_needed

        if args.genome_fasta:
            try:
                fasta = pysam.FastaFile(args.genome_fasta)
                query_chrom = abs_chrom
                if query_chrom not in fasta.references:
                    alt_chrom = f"chr{query_chrom}" if not query_chrom.startswith("chr") else query_chrom.replace("chr", "")
                    if alt_chrom in fasta.references: query_chrom = alt_chrom

                if query_chrom in fasta.references:
                    chrom_len = fasta.get_reference_length(query_chrom)
                    if abs_strand == '+':
                        fetch_start = abs_end_pos
                        fetch_end = min(abs_end_pos + extra_needed, chrom_len)
                        true_extra = fetch_end - fetch_start
                        fetched_seq = fasta.fetch(query_chrom, fetch_start, fetch_end).upper()
                        appended_seq = fetched_seq
                        extra_needed = true_extra
                    else: 
                        fetch_end = abs_start_pos - 1
                        fetch_start = max(0, fetch_end - extra_needed)
                        true_extra = fetch_end - fetch_start
                        fetched_seq_raw = fasta.fetch(query_chrom, fetch_start, fetch_end).upper()
                        appended_seq = reverse_complement(fetched_seq_raw)
                        extra_needed = true_extra
                fasta.close()
            except Exception as e:
                print(f"Warning: Failed to fetch downstream sequence via pysam: {e}")
                appended_seq = 'N' * extra_needed
        else:
            appended_seq = 'N' * extra_needed
            
        if appended_seq:
            extra_seq_arr = getDNAseq_array(appended_seq)
            dna_sequence_arr = np.concatenate((dna_sequence_arr, extra_seq_arr))
            
            if not gtf_gene.empty:
                max_pos_index = gtf_gene['position'].max()
                new_row = pd.DataFrame([{
                    'chr': abs_chrom, 'start': int(fetch_start + 1), 'end': int(fetch_end),
                    'gene_id': base_filename, 'feature': 'readthrough',
                    'position': int(max_pos_index + 1), 'strand': abs_strand,
                    'sequence': appended_seq, 'pas_coordinate': int(PAS_coord)
                }])
                gtf_gene = pd.concat([gtf_gene, new_row], ignore_index=True)
                atomic_to_csv(gtf_gene, gtf_path_in, sep='\t', index=False)

    for mol in molecules_to_emit:
        if mol['cleaved']:
            base_polyA = rng.integers(polyA_len_min, polyA_len_max + 1) if polyA_len_max >= polyA_len_min else 0
            polyA_len = int(base_polyA * rng.uniform(0.9, 1.1))
            
            if args.experiment_type in ['metaboliclabeling', 'ttseq','chromatin_associated']:
                emit_molecule(0, PAS_coord, mol['start_label_pos'], "_cleaved_body", polyA_len, mol['initiation_sec'], mol['molecule_id'], mol['time_since_cleavage_sec'])
                emit_molecule(mol['deg_front'], mol['stop_label_pos'], mol['start_label_pos'], "_readthrough_tail", 0, mol['initiation_sec'], mol['molecule_id'], mol['time_since_cleavage_sec'])
            elif args.experiment_type in ['netseq', 'proseq']:
                emit_molecule(mol['deg_front'], mol['stop_label_pos'], mol['start_label_pos'], "_readthrough_tail", 0, mol['initiation_sec'], mol['molecule_id'], mol['time_since_cleavage_sec'])
        else:
            emit_molecule(0, mol['stop_label_pos'], mol['start_label_pos'], "_uncleaved", 0, mol['initiation_sec'], mol['molecule_id'], 0.0)

    #export molecules
    df_for_export = pd.DataFrame(data_list, columns=[
        "initiation_time", 'molecule_id', 'strand', "sub_rate_percent_range", 
        'start_label_pos', 'stop_label_pos', 'converted_positions', 
        'incorporated_positions', 'percentage_seq_err', 'seq_err_positions', 'spliced_introns', 'polyA_length', 'time_since_cleaved_min'
    ])

    n_emitted = len(df_for_export)
    if n_emitted > nascent_to_simulate:
        df_for_export = df_for_export.sample(n=nascent_to_simulate, random_state=42)
    simulated_cells = scaling_factor * len(df_for_export) / n_emitted if n_emitted > 0 else 0.0

    atomic_to_csv(df_for_export, output_filename, sep='\t', index=False)
    atomic_to_csv(pd.DataFrame([{
        "gene_id": base_filename,
        "simulated_cells": simulated_cells,
        "initiations_per_min": 60.0 / initiation_interval_sec,
    }]), os.path.join(os.path.dirname(args.region_file), f"{base_filename}_cells.tsv"), sep="\t", index=False)

    #capture mean seq_err across all nascent molecules for the ground truth
    _mean_seq_err = float(df_for_export['percentage_seq_err'].mean()) if not df_for_export.empty else 0.0

    if args.bkg_molecules > 0.0 and bkg_molecules_to_simulate > 0:
        bg_data_list = []
        output_filename_bg = os.path.join(args.o, base_filename + "_background.tsv.gz")
        
        try:
            if args.nosplicing:
                bg_len = gtf_gene[gtf_gene['feature'] != 'readthrough']['sequence'].str.len().sum()
            else:
                exon_gtf = gtf_gene[gtf_gene['feature'] == 'exon']
                bg_len = exon_gtf['sequence'].str.len().sum()
        except:
            bg_len = PAS_coord
        
        bg_ids = generate_ids_bulk(bkg_molecules_to_simulate, size=12, rng=rng)
        
        for i in range(bkg_molecules_to_simulate):
            seq_err_rate = rng.uniform(seqerr_range[0], seqerr_range[1])
            n_err = rng.binomial(bg_len, seq_err_rate)
            err_pos = []
            if n_err > 0:
                err_locs = rng.choice(bg_len, size=n_err, replace=False)
                err_pos = [int(x) for x in err_locs]
            
            base_bg_polyA = rng.integers(polyA_len_min, polyA_len_max + 1) if polyA_len_max >= polyA_len_min else 0
            bg_polyA_len = int(base_bg_polyA * rng.uniform(0.9, 1.1))
            
            bkg_row_data = [
                0, bg_ids[i], df_nt.loc[0, 'strand'], 0, 0, bg_len,
                '[]', '[]', seq_err_rate, err_pos, '[]', bg_polyA_len, 0.0
            ]
            bg_data_list.append(bkg_row_data)

        bg_data_list_for_export = pd.DataFrame(bg_data_list, columns=[
            "initiation_time", 'molecule_id', 'strand', "sub_rate_percent_range", 
            'start_label_pos', 'stop_label_pos', 'converted_positions', 
            'incorporated_positions', 'percentage_seq_err', 'seq_err_positions', 'spliced_introns', 'polyA_length', 'time_since_cleaved_min'
        ])
        atomic_to_csv(bg_data_list_for_export, output_filename_bg, sep='\t', index=False)

    parent_dir = os.path.dirname(args.o.rstrip('/'))
    temp_dir = os.path.join(parent_dir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)

    init_path = os.path.join(temp_dir, f"temp_{base_filename}_initiation.tsv")
    df_summary = pd.DataFrame([{"gene_id": base_filename, "initiation_rate": initiation_interval_sec}])
    atomic_to_csv(df_summary, init_path, sep="\t", index=False)

    elong_path = os.path.join(temp_dir, f"temp_{base_filename}_elongation.tsv")
    df_elongation = pd.DataFrame([{"gene_id": base_filename, "mean_elongation_rate": round(float(mean_elong_rate), 4)}])
    atomic_to_csv(df_elongation, elong_path, sep="\t", index=False)

    seq_err_path = os.path.join(temp_dir, f"temp_{base_filename}_seq_err.tsv")
    df_seq_err = pd.DataFrame([{"gene_id": base_filename, "mean_seq_err": round(_mean_seq_err, 8)}])
    atomic_to_csv(df_seq_err, seq_err_path, sep="\t", index=False)

    pas_elong_path = os.path.join(temp_dir, f"temp_{base_filename}_pas_elong.tsv")
    pas_elong_kb_min = round(gene_pas_rate * 60 / 1000, 4)  #convert nt/sec back to kb/min
    df_pas_elong = pd.DataFrame([{"gene_id": base_filename, "PAS_elongation_rate": pas_elong_kb_min}])
    atomic_to_csv(df_pas_elong, pas_elong_path, sep="\t", index=False)

    if not args.nocleavage and (gene_cleavage_half_life > 0 or args.immediate_cleavage):
        cleavage_path = os.path.join(temp_dir, f"temp_{base_filename}_cleavage.tsv")
        
        hl_out = 0.0 if args.immediate_cleavage else (gene_cleavage_half_life / 60.0)
        
        df_cleavage = pd.DataFrame([{
            "gene_id": base_filename,
            "pas_coordinates": int(PAS_coord),
            "half_life": hl_out
        }])
        atomic_to_csv(df_cleavage, cleavage_path, sep="\t", index=False)
