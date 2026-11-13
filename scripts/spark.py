import argparse
import subprocess
import random
import os
import glob
import sys
import numpy as np
import pandas as pd
import gzip
import logging
from tqdm import tqdm
from datetime import datetime
from multiprocessing import Pool
from spark_labeling import parse_conversion_map

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
POOL_MARGIN = 1.2
args = None 

PRESETS = {
    "default": {
    },
    "DRB4sUseq": {        
        "experiment_type": "metaboliclabeling",
        "drb": True,
        "labeling_base": "T",
        "nt_inc_prob": "0.005,0.023",
        "bkg_molecules": 0.15
    },
    "4sUseq": {
        "experiment_type": "metaboliclabeling",
        "labeling_base": "T",
        "nt_inc_prob": "0.005,0.023",
        "bkg_molecules": 0.15
    },
    "TTseq": {
        "experiment_type": "ttseq",
        "labeling_base": "T",
        "nt_inc_prob": "0.005,0.023",
        "bkg_molecules": 0.10
    },
    "BrUseq": {        
        "experiment_type": "metaboliclabeling",
        "labeling_base": "T",
        "nt_inc_prob": "0.001,0.005",
        "bkg_molecules": 0.15
    },
    "NETseq": {                
        "experiment_type": "netseq",
        "bkg_molecules": 0.10
    },
    "proseq": {        
        "experiment_type": "proseq",
        "labeling_base": "A,C,G,T",
        "bkg_molecules": 0.02
    },
    "chromatinassociated": {        
        "experiment_type": "chromatin_associated",
        "bkg_molecules": 0.20
    },
    "TT-TimelapseSeq": {
        "experiment_type": "ttseq",
        "labeling_base": "T",
        "sub_base": "C",
        "nt_inc_prob": "0.005,0.023",
        "bkg_molecules": 0.05
    },
    "eSLAM-seq": {
        "experiment_type": "metaboliclabeling",
        "labeling_base": "T",
        "sub_base": "C",
        "nt_inc_prob": "0.005,0.023",
        "bkg_molecules": 0.10
    },
    "polyA": {
        "experiment_type": "metaboliclabeling",
        "bkg_molecules": 1
    }
}

def init_worker(shared_args):
    global args
    args = shared_args

def worker_pool():
    return Pool(processes=args.threads, initializer=init_worker, initargs=(args,))

def c_log(msg):
    print(f"[SPARK] {msg}")
    logging.info(msg)

def run_cmd(cmd):
    logging.info(f"CMD: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        logging.error(f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        raise

def labeling_args():
    parts = [f"--labeling_base {args.labeling_base}"]
    if args.sub_base:
        parts.append(f"--sub_base {args.sub_base}")
    return parts

def remove_files(paths):
    for path in paths:
        try: os.remove(path)
        except FileNotFoundError: pass

def remove_gene_intermediates(gene_id):
    output_dir = args.o.rstrip('/')
    remove_files([
        os.path.join(output_dir, "gtf", f"{gene_id}.tsv.gz"),
        os.path.join(output_dir, "mRNA", f"{gene_id}.tsv.gz"),
        os.path.join(output_dir, "mRNA", f"{gene_id}_background.tsv.gz"),
        os.path.join(output_dir, "proseqmRNAs", f"{gene_id}_proseq.tsv.gz"),
        os.path.join(output_dir, "netseqmRNAs", f"{gene_id}_netseq.tsv.gz"),
    ])

def sweep_intermediates(output_dir, targets):
    for subdir, pattern in targets:
        remove_files(glob.glob(os.path.join(output_dir, subdir, pattern)))
    for subdir in dict.fromkeys(s for s, _ in targets):
        try: os.rmdir(os.path.join(output_dir, subdir))
        except OSError: pass

def count_reads(fastq_file):
    if not os.path.exists(fastq_file): return 0
    try:
        ps = subprocess.Popen(f"gzip -dc {fastq_file} | wc -l", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, _ = ps.communicate()
        lines = int(stdout.strip())
        return lines // 4
    except Exception as e:
        logging.error(f"Error counting reads in {fastq_file}: {e}")
        return 0

def downsample_library(r1_file, r2_file, target_depth, seed):
    current_reads = count_reads(r1_file)
    if current_reads <= target_depth:
        c_log(f"Total reads ({current_reads}) <= target depth ({target_depth}). No downsampling needed.")
        return

    c_log(f"Downsampling from {current_reads} to approx {target_depth} reads...")
    
    frac = target_depth / current_reads
    
    r1_out = r1_file.replace(".fastq.gz", ".tmp.fastq.gz")
    r2_out = r2_file.replace(".fastq.gz", ".tmp.fastq.gz") if r2_file else None
    
    rng = random.Random(seed) if seed is not None else random.Random()
    
    try:
        with gzip.open(r1_file, 'rt') as f1_in, gzip.open(r1_out, 'wt') as f1_out:
            if r2_file:
                with gzip.open(r2_file, 'rt') as f2_in, gzip.open(r2_out, 'wt') as f2_out:
                    while True:
                        r1_rec = [f1_in.readline() for _ in range(4)]
                        if not r1_rec[0]: break 
                        
                        r2_rec = [f2_in.readline() for _ in range(4)]
                        
                        if rng.random() < frac:
                            for line in r1_rec: f1_out.write(line)
                            for line in r2_rec: f2_out.write(line)
            else:
                while True:
                    rec = [f1_in.readline() for _ in range(4)]
                    if not rec[0]: break 
                    
                    if rng.random() < frac:
                        for line in rec: f1_out.write(line)

        os.replace(r1_out, r1_file)
        if r2_file:
            os.replace(r2_out, r2_file)
        
        c_log("Downsampling Complete. Files updated.")

    except Exception as e:
        logging.error(f"Downsampling failed: {e}")
        if os.path.exists(r1_out): os.remove(r1_out)
        if r2_file and os.path.exists(r2_out): os.remove(r2_out)

def consolidate_ground_truth(output_dir, args):
    temp_dir = os.path.join(output_dir, "temp")
    master_path = os.path.join(output_dir, "ground_truth_per_gene.tsv")
    
    #load temp files for ground truth
    def load_temp_files(pattern):
        files = glob.glob(os.path.join(temp_dir, pattern))
        if not files:
            return pd.DataFrame(), []
        df_list = [pd.read_csv(f, sep='\t') for f in files]
        df = pd.concat(df_list, ignore_index=True).drop_duplicates()
        if 'gene_id' in df.columns:
            df = df[~df['gene_id'].astype(str).str.endswith('_background')]
        return df, files

    init_df, init_files = load_temp_files("temp_*_initiation.tsv")
    cleavage_df, cleavage_files = load_temp_files("temp_*_cleavage.tsv")
    rpkm_df, rpkm_files = load_temp_files("temp_*_rpkm.tsv")
    if not rpkm_df.empty and 'rpkm' in rpkm_df.columns:
        rpkm_sum = rpkm_df['rpkm'].sum()
        rpkm_df['tpm'] = (1e6 * rpkm_df['rpkm'] / rpkm_sum) if rpkm_sum > 0 else 0.0
    elongation_df, elong_files = load_temp_files("temp_*_elongation.tsv")
    tss_pause_df, tss_pause_files = load_temp_files("temp_*_tss_pause.tsv")
    seq_err_df, seq_err_files = load_temp_files("temp_*_seq_err.tsv")
    pas_elong_df, pas_elong_files = load_temp_files("temp_*_pas_elong.tsv")
    chromatin_hl_df, chromatin_hl_files = load_temp_files("temp_*_chromatin_halflife.tsv")

    #rename half_life and drop extra columns from cleavage
    if not cleavage_df.empty and 'half_life' in cleavage_df.columns:
        cleavage_df.rename(columns={'half_life': 'cleavage_half_life'}, inplace=True)
    if not cleavage_df.empty and 'pas_coordinates' in cleavage_df.columns:
        cleavage_df.drop(columns=['pas_coordinates'], inplace=True)

    #gather all new data that exists
    new_dataframes = [df for df in [init_df, cleavage_df, rpkm_df, elongation_df, tss_pause_df, seq_err_df, pas_elong_df, chromatin_hl_df] if not df.empty]
    
    if not new_dataframes:
        return 
        
    c_log("Consolidating ground truth metrics into ground_truth_per_gene.tsv...")
    
    #outer merge all new metrics on gene_id
    new_master = new_dataframes[0]
    for df in new_dataframes[1:]:
        new_master = pd.merge(new_master, df, on="gene_id", how="outer")

    desired_cols = [
        'gene_id',
        'initiation_rate', 'tss_pause_position', 'tss_pause_duration_min',
        'min_elongation_rate', 'max_elongation_rate', 'mean_elongation_rate',
        'num_pauses', 'pause_duration_min', 'pause_positions',
        'cleavage_half_life', 'PAS_elongation_rate',
        'rpkm', 'tpm', 'mean_seq_err',
    ]
    if getattr(args, 'experiment_type', None) == 'chromatin_associated':
        desired_cols.append('chromatin_release_halflife_min')
    final_cols = [c for c in desired_cols if c in new_master.columns]
    new_master = new_master[final_cols]

    #ff an older master ground truth exists, update it
    if os.path.exists(master_path):
        old_master = pd.read_csv(master_path, sep='\t')
        old_master.set_index("gene_id", inplace=True)
        new_master.set_index("gene_id", inplace=True)
        combined = new_master.combine_first(old_master).reset_index()
    else:
        combined = new_master

    #format column order
    combined_cols = ['gene_id'] + [c for c in desired_cols if c in combined.columns and c != 'gene_id']
    combined = combined[combined_cols]
    combined = combined[~combined['gene_id'].astype(str).str.endswith('_background')]

    combined.to_csv(master_path, sep='\t', index=False)
    logging.info(f"Ground Truth saved to {master_path}")

    #cleanup the individual temp files
    for f in (init_files + cleavage_files + rpkm_files + elong_files + tss_pause_files + seq_err_files + pas_elong_files + chromatin_hl_files):
        try: os.remove(f)
        except OSError: pass

def process_rates_file(tsv):
    cmd_rates = [
        f"python {os.path.join(SCRIPT_DIR, 'rates_per_region.py')}",
        f"--tsv {tsv}", f"--elong_rate_range {args.elong_rate_range}",
        f"--pause_time {args.pause_time}", f"--o {args.o.rstrip('/')}",
        f"--promoter_pause_position {args.promoter_pause_position}",
        f"--promoter_pause_duration {args.promoter_pause_duration}"
    ]
    if args.flat_rates: cmd_rates.append("--flat_rates")
    if args.seed is not None: cmd_rates.extend(['--seed', str(args.seed)])
    if args.region_size_range: cmd_rates.extend(['--region_size_range', args.region_size_range])
    elif args.num_regions: cmd_rates.extend(['--num_regions', str(args.num_regions)])
    else: raise ValueError("You must provide either --region_size_range or --num_regions")
    if args.num_pauses is not None and args.num_pauses > 0:
        cmd_rates.append(f"--num_pauses {args.num_pauses}")
    run_cmd(" ".join(cmd_rates))

def process_mrna_file(region_file, mrna_count=None, reuse_parameters=False):
    output_dir = args.o.rstrip("/")
    dir_region_files = f"{output_dir}/rate_per_gene/"
    premRNA_out_dir = os.path.join(output_dir, "pre-mRNA")
    
    gene_id = os.path.basename(region_file).split("_")[0]
    nt_file = os.path.join(dir_region_files, f"{gene_id}_RatesandTraversalTimes.gtf")

    if not os.path.isfile(nt_file):
        logging.warning(f"nt file missing for {gene_id}, skipping.")
        return

    cmd_label = [
        f"python {os.path.join(SCRIPT_DIR, 'mRNA_generator.py')}",
        f"--region_file {region_file}", f"--nt_file {nt_file}",
        f"--experiment_time {args.experiment_time}", f"--mRNAcount {mrna_count or args.mRNAcount}",
        f"--seq_err {args.seq_err}", f"--nt_inc_rate {args.nt_inc_prob}",
        f"--subs_rate {args.subs_prob}", f"--labeling_base {args.labeling_base}",
        f"--bkg_molecules {args.bkg_molecules}", f"--initi_rate {args.initi_rate}",
        f"--experiment_type {args.experiment_type}",
        f"--cleavage_half_life {args.cleavage_half_life}",
        f"--pas_elong_rate {args.pas_elong_rate}",
        f"--readthrough_deg_rate {args.readthrough_deg_rate}",
        f"--polyA_length {args.polyA_length}", 
        f"--genome_fasta {args.genome_fasta}",
        f"--o {premRNA_out_dir}/"
    ]
    if args.seed is not None: cmd_label.extend(['--seed', str(args.seed)])
    if args.drb: cmd_label.append('--drb')
    if args.sub_base: cmd_label.extend(['--sub_base', args.sub_base])
    if args.nosplicing: cmd_label.append("--nosplicing") 
    if args.nocleavage: cmd_label.append("--nocleavage")
    if args.immediate_cleavage: cmd_label.append("--immediate_cleavage")
    if reuse_parameters: cmd_label.append("--reuse_parameters")
    run_cmd(" ".join(cmd_label))

def process_chromatin_isolation(mrna_file, reuse_parameters=False):
    script_path = os.path.join(SCRIPT_DIR, 'chromatin_isolation.py')
    carna_out_dir = os.path.join(args.o.rstrip('/'), 'caRNA_filtered')
    os.makedirs(carna_out_dir, exist_ok=True)
    
    cmd = [
        f"python {script_path}",
        f"--input_df {mrna_file}",
        f"--release_half_life {args.chromatin_release_half_life}",
        f"--o {carna_out_dir}"
    ]
    if args.seed is not None:
        cmd += ["--seed", str(args.seed)]
    if reuse_parameters:
        cmd.append("--reuse_parameters")
        
    run_cmd(" ".join(cmd))
    if args.remove_temp: remove_files([mrna_file])

def process_splicing(mrna_file):
    script_path = os.path.join(SCRIPT_DIR, 'spark_spliceosome.py')
    splice_out_dir = os.path.join(args.o.rstrip('/'), 'mRNA')
    halflife_out_dir = os.path.join(args.o.rstrip('/'), 'intron_half_lives')

    os.makedirs(splice_out_dir, exist_ok=True)
    os.makedirs(halflife_out_dir, exist_ok=True)

    base = os.path.splitext(os.path.basename(mrna_file))[0].split('.')[0]

    cmd = [
        f"python {script_path}",
        f"--mRNA_df {mrna_file}",
        f"--intron_half_life {args.intron_half_life}",
        f"--o {args.o.rstrip('/')}"
    ]
    if args.seed is not None:
        cmd += ["--seed", str(args.seed)]
    if args.nosplicing:
        cmd += ["--nosplicing"]
    if args.mRNA_coordinates:
        cmd += ["--mRNA_coordinates"]

    run_cmd(" ".join(cmd))
    if args.remove_temp: remove_files([mrna_file])

def mrna_to_proseq(tsv):
    proseq_mrnas_dirout = os.path.join(args.o.rstrip('/'), "proseqmRNAs")
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_proseq.py')}", f"--input_df {tsv}", f"--o {proseq_mrnas_dirout}"]
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))

def read_phase_args(reads):
    if reads is None:
        return ["--prepare"]
    return [f"--nascent_reads {reads[0]}", f"--bkg_reads {reads[1]}"]

def allocate_reads(gene_ids):
    output_dir = args.o.rstrip('/')
    rows = []
    for gene_id in gene_ids:
        units_path = os.path.join(output_dir, "temp", "pools", f"{gene_id}_units.tsv")
        cells_path = os.path.join(output_dir, "rate_per_gene", f"{gene_id}_cells.tsv")
        if not (os.path.exists(units_path) and os.path.exists(cells_path)):
            logging.warning(f"Read pool or simulated cell count missing for {gene_id}, skipping.")
            continue
        row = pd.read_csv(units_path, sep='\t').iloc[0].to_dict()
        row.update(pd.read_csv(cells_path, sep='\t').iloc[0].to_dict())
        gtf_df = pd.read_csv(os.path.join(output_dir, "gtf", f"{gene_id}.tsv.gz"), sep='\t')
        row['exon_kb'] = gtf_df.loc[gtf_df['feature'] == 'exon', 'sequence'].str.len().sum() / 1000
        row['gene_id'] = gene_id
        rows.append(row)
    if not rows:
        return {}
    genes = pd.DataFrame(rows)

    nascent_weight = (genes['nascent_units'] / genes['simulated_cells'].where(genes['simulated_cells'] > 0)).fillna(0.0)
    if args.bkg_molecules < 1.0 and nascent_weight.sum() > 0:
        nascent_reads = (args.seq_depth * (1 - args.bkg_molecules) * nascent_weight / nascent_weight.sum()).round()
        bkg_reads = (nascent_reads * args.bkg_molecules / (1 - args.bkg_molecules)).round()
    else:
        bkg_weight = (genes['initiations_per_min'] * genes['bg_units'] / genes['bg_molecules'].where(genes['bg_molecules'] > 0)).fillna(0.0)
        nascent_reads = pd.Series(0.0, index=genes.index)
        bkg_reads = (args.seq_depth * bkg_weight / bkg_weight.sum()).round() if bkg_weight.sum() > 0 else nascent_reads.copy()

    genes['nascent_reads'] = nascent_reads.astype(int)
    genes['bkg_reads'] = bkg_reads.astype(int)
    genes['rpkm'] = ((genes['nascent_reads'] + genes['bkg_reads']) / genes['exon_kb'].where(genes['exon_kb'] > 0) / (args.seq_depth / 1e6)).fillna(0.0)

    temp_dir = os.path.join(output_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    for gene in genes.itertuples():
        pd.DataFrame([{"gene_id": gene.gene_id, "rpkm": gene.rpkm}]).to_csv(
            os.path.join(temp_dir, f"temp_{gene.gene_id}_rpkm.tsv"), sep='\t', index=False)
    c_log(f"Allocated {int(genes['nascent_reads'].sum())} nascent and {int(genes['bkg_reads'].sum())} background reads across {len(genes)} genes")
    return {gene.gene_id: (int(gene.nascent_reads), int(gene.bkg_reads)) for gene in genes.itertuples()}

def enlargement_factors(allocation):
    pool_dir = os.path.join(args.o.rstrip('/'), "temp", "pools")
    factors = {}
    for gene_id, (nascent_reads, bkg_reads) in allocation.items():
        units = pd.read_csv(os.path.join(pool_dir, f"{gene_id}_units.tsv"), sep='\t').iloc[0]
        need = max(nascent_reads / units['nascent_units'] if units['nascent_units'] > 0 else 0.0,
                   bkg_reads / units['bg_units'] if units['bg_units'] > 0 else 0.0)
        if need > 1:
            factors[gene_id] = need * POOL_MARGIN
    return factors

def resimulate_gene(task):
    gene_id, factor = task
    output_dir = args.o.rstrip('/')
    premRNA_out_dir = os.path.join(output_dir, "pre-mRNA")
    os.makedirs(premRNA_out_dir, exist_ok=True)
    region_file = os.path.join(output_dir, "rate_per_gene", f"{gene_id}_VariableElongationRateRegions.gtf")
    process_mrna_file(region_file, mrna_count=int(np.ceil(args.mRNAcount * factor)), reuse_parameters=True)

    mrna_files = [path for path in (os.path.join(premRNA_out_dir, f"{gene_id}.tsv.gz"),
                                    os.path.join(premRNA_out_dir, f"{gene_id}_background.tsv.gz")) if os.path.exists(path)]
    if args.experiment_type == 'chromatin_associated':
        for path in mrna_files:
            process_chromatin_isolation(path, reuse_parameters=True)
        mrna_files = [os.path.join(output_dir, "caRNA_filtered", os.path.basename(path)) for path in mrna_files]
    for path in mrna_files:
        process_splicing(path)

def generate_reads(worker, tsv_files, suffix, desc, rebuild_pool=None):
    gene_ids = [os.path.basename(f).replace(suffix, "") for f in tsv_files]
    tsv_by_gene = dict(zip(gene_ids, tsv_files))
    with worker_pool() as pool:
        list(tqdm(pool.imap_unordered(worker, [(f, None) for f in tsv_files]), total=len(tsv_files), desc=f"{desc} (pools)"))
    allocation = allocate_reads(gene_ids)

    factors = enlargement_factors(allocation)
    if factors:
        c_log(f"{len(factors)} genes need more reads than their molecule pools hold; re-simulating them with up to {max(factors.values()):.1f}x more molecules")
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(resimulate_gene, factors.items()), total=len(factors), desc=f"{desc} (larger pools)"))
            if rebuild_pool is not None:
                mrna_dir = os.path.join(args.o.rstrip('/'), "mRNA")
                list(pool.imap_unordered(rebuild_pool, [os.path.join(mrna_dir, f"{g}.tsv.gz") for g in factors]))
            list(pool.imap_unordered(worker, [(tsv_by_gene[g], None) for g in factors]))

    tasks = [(f, allocation[g]) for f, g in zip(tsv_files, gene_ids) if g in allocation]
    with worker_pool() as pool:
        list(tqdm(pool.imap_unordered(worker, tasks), total=len(tasks), desc=desc))

def proseq_to_fastq(task):
    tsv, reads = task
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_proseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--s {args.s}", f"--o {args.o.rstrip('/')}/", f"--bkg_molecules {args.bkg_molecules}"]
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    cmd.extend(read_phase_args(reads))
    run_cmd(" ".join(cmd))
    if reads is not None and args.remove_temp: remove_gene_intermediates(os.path.basename(tsv).replace("_proseq.tsv.gz", ""))

def mrna_to_netseq(tsv):
    netseq_mrnas_dirout = os.path.join(args.o.rstrip('/'), "netseqmRNAs")
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_netseq.py')}", f"--input_df {tsv}", f"--o {netseq_mrnas_dirout}", f"--bkg_molecules {args.bkg_molecules}"]
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))
    
def netseq_to_fastq(task):
    tsv, reads = task
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_netseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--s {args.s}", f"--o {args.o.rstrip('/')}/"]
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    cmd.extend(read_phase_args(reads))
    run_cmd(" ".join(cmd))
    if reads is not None and args.remove_temp: remove_gene_intermediates(os.path.basename(tsv).replace("_netseq.tsv.gz", ""))

def nascent_to_fastq(task):
    tsv, reads = task
    output_dir = args.o.rstrip('/')
    if args.seq_tech == "longread":
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_long_read.py')}",
               f"--input_df {tsv}", f"--seq_type {args.seq_type}", f"-o {output_dir}",
               f"--bkg_molecules {args.bkg_molecules}", f"--experiment_type {args.experiment_type}"]
    else: 
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_short_read.py')}",
               f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
               f"--seq_type {args.seq_type}", f"--threads 1", f"--s {args.s}", f"--o {output_dir}/", f"--experiment_time {args.experiment_time}", f"--bkg_molecules {args.bkg_molecules}"]
        if args.fragments: cmd.append('--fragments')
        if args.sizeselectiontype: cmd.extend(["--sizeselectiontype", args.sizeselectiontype])
        if args.no_fragmentation: cmd.append("--no_fragmentation")
    cmd.extend(labeling_args())
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    cmd.extend(read_phase_args(reads))
    run_cmd(" ".join(cmd))
    if reads is not None and args.remove_temp: remove_gene_intermediates(os.path.basename(tsv).replace(".tsv.gz", ""))

def ttseq_to_fastq(task):
    tsv, reads = task
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_reads_ttseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--s {args.s}", f"--o {args.o.rstrip('/')}/", f"--bkg_molecules {args.bkg_molecules}"]
    if args.fragments: cmd.append('--fragments')
    if args.sizeselectiontype: cmd.extend(["--sizeselectiontype", args.sizeselectiontype])
    if args.no_fragmentation: cmd.append("--no_fragmentation")
    cmd.extend(labeling_args())
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    cmd.extend(read_phase_args(reads))
    run_cmd(" ".join(cmd))
    if reads is not None and args.remove_temp: remove_gene_intermediates(os.path.basename(tsv).replace(".tsv.gz", ""))

def carna_to_fastq(task):
    tsv, reads = task
    output_dir = args.o.rstrip('/')
    if args.seq_tech == "longread":
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_long_read.py')}",
               f"--input_df {tsv}", f"--seq_type {args.seq_type}", f"-o {output_dir}",
               f"--bkg_molecules {args.bkg_molecules}", f"--experiment_type chromatin_associated"]
    else: 
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_short_read.py')}",
               f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
               f"--seq_type {args.seq_type}", f"--threads 1", f"--s {args.s}", f"--o {output_dir}/", f"--experiment_time {args.experiment_time}", f"--bkg_molecules {args.bkg_molecules}", f"--experiment_type chromatin_associated"]
        if args.fragments: cmd.append('--fragments')
        if args.sizeselectiontype: cmd.extend(["--sizeselectiontype", args.sizeselectiontype])
        if args.no_fragmentation: cmd.append("--no_fragmentation")
    cmd.extend(labeling_args())
    if args.seed is not None: cmd.extend(['--seed', str(args.seed)])
    cmd.extend(read_phase_args(reads))
    run_cmd(" ".join(cmd))
    if reads is not None and args.remove_temp: remove_gene_intermediates(os.path.basename(tsv).replace(".tsv.gz", ""))

def main():
    global args
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group_general_arguments = parser.add_argument_group('general arguments')
    group_presets=parser.add_argument_group('presets')
    group_tsv_generation = parser.add_argument_group('gene selection')
    group_rates_per_region = parser.add_argument_group('rates per region')
    group_mRNAgeneration = parser.add_argument_group('pre-mRNA generation')
    group_mRNAprocessing = parser.add_argument_group('pre-mRNA processing')
    group_seq_tech = parser.add_argument_group('sequencing strategy')

    group_general_arguments.add_argument("--mode", choices=["fullpipeline", "gene_selection", "rates_per_region", "mRNAgeneration", "readgeneration"], help="pipeline step to run",default='fullpipeline')
    group_general_arguments.add_argument("-o", type=str, default="./", help="output directory")
    group_general_arguments.add_argument("--seed", type=int, help="random seed")
    group_general_arguments.add_argument("--experiment_type", default="metaboliclabeling", choices=["metaboliclabeling", "proseq", "ttseq", "netseq", "chromatin_associated"], help="type of experiment to simulate")
    group_general_arguments.add_argument("--threads", type=int, default=1, help="number of threads")
    group_general_arguments.add_argument("--remove_temp", action="store_true", help="delete intermediate files as soon as no later step needs them, keeping only the final libraries, ground truth files and the log. Files still needed by a later --mode step are kept until that step runs")
    
    group_presets.add_argument("--preset", choices=list(PRESETS.keys()), default="default", help="load a predefined set of parameters")
    group_presets.add_argument("--show_presets", action="store_true", help="display all available presets and their default values, then exit.")

    group_tsv_generation.add_argument("--gtf", type=str, help="genome annotation file")
    group_tsv_generation.add_argument("--genome_fasta", type=str, help="genome fasta file")
    group_tsv_generation.add_argument("--protein_coding_only", action="store_true", help="only use protein-coding genes")
    group_tsv_generation.add_argument("--n_gene_clusters", type=int, default=1, help="number of gene clusters to simulate")
    group_tsv_generation.add_argument("--n", type=int, default=0, help="number of genes to simulate (set to 0 to simulate all genes in the GTF)")

    group_rates_per_region.add_argument("--region_size_range", type=str, help="range of elongation rate region sizes to split the gene into (default: 500,5000)")
    group_rates_per_region.add_argument('--num_regions', type=int, help='number of elongation rate regions to split the gene into')
    group_rates_per_region.add_argument("--elong_rate_range", default="500,5000", help="range of elongation rates in nt/min")
    group_rates_per_region.add_argument('--num_pauses', type=int, default=0, help='number of pausing events across the gene')
    group_rates_per_region.add_argument('--pause_time',default="0.1,0.5", type=str, help="length of the pausing event in minutes")
    group_rates_per_region.add_argument("--flat_rates", action="store_true", help="all nucleotides have the same elongation rate")
    group_rates_per_region.add_argument('--promoter_pause_position', type=str, default='30,50', help='range of position for TSS pausing')
    group_rates_per_region.add_argument('--promoter_pause_duration', type=str, default='0,0', help='length of promoter pausing in minutes')
    group_rates_per_region.add_argument('--initi_rate',type=str,default='5,200',help='range of seconds between initiations, sampled log-uniformly per gene. Sets each gene\'s expression')

    group_mRNAgeneration.add_argument("--experiment_time", type=int, default=15, help="experiment duration in minutes")
    group_mRNAgeneration.add_argument("--seq_err", default="0.0001,0.0002", help="sequencing error rate")
    group_mRNAgeneration.add_argument('--nt_inc_prob',type=str,default='0.005,0.023',help='range of nucleotide incorporation probability')
    group_mRNAgeneration.add_argument('--subs_prob',type=str,default='0.95,1',help='range of nucleotide substitution probability if simulating nucleotide recoding')
    group_mRNAgeneration.add_argument('--labeling_base', type=str, default='T', help='base analog(s) incorporated into the nascent RNA. A single base or a comma-separated list, e.g. "T" or "A,T"')
    group_mRNAgeneration.add_argument('--sub_base', type=str, help='identity of the base(s) after conversion, matched positionally to --labeling_base, e.g. "C" or "G,C" (A->G and T->C). A single value is applied to every labeling base. Omit for incorporation without recoding')
    group_mRNAgeneration.add_argument('--bkg_molecules', type=float, default=0,help="fraction of each gene's reads that come from background (non-nascent) molecules")
    group_mRNAgeneration.add_argument("--drb", action="store_true", help="simulate a DRB treatment experiment for transcription synchronization")
    group_mRNAgeneration.add_argument("--mRNAcount", type=int, default=5000, help="size of the molecule pool simulated per gene, including background if specified. Sets simulation size, not expression")

    group_mRNAprocessing.add_argument("--cleavage_half_life", type=str, default="0.1,10.0", help="range of cleavage half-lives in min. Randomly chosen per gene")
    group_mRNAprocessing.add_argument("--chromatin_release_half_life", type=str, default="2.0,5.0", help="Min,Max chromatin release half-life in minutes")
    group_mRNAprocessing.add_argument("--pas_elong_rate", type=str, default="0.5,2.0", help="range of elongation rates in kb/min after passing the polyA site")
    group_mRNAprocessing.add_argument("--readthrough_deg_rate", type=float, default=1.56, help="degradation rate of read-through transcripts after cleavage in kb/min")
    group_mRNAprocessing.add_argument("--nocleavage", action="store_true", help="do not simulate cleavage")
    group_mRNAprocessing.add_argument("--immediate_cleavage", action="store_true", help="simulate immediate cleavage at the polyA site without a half-life delay")
    group_mRNAprocessing.add_argument('--intron_half_life', default='1,10',help="range of intron half life in minutes for the introns to simulate splicing")
    group_mRNAprocessing.add_argument("--nosplicing", action="store_true", help="do not simulate splicing")
    group_mRNAprocessing.add_argument("--polyA_length", type=str, default="50,150", help="range of polyA tail lengths to add to cleaved mRNAs (base range, subject to +/- 10%% jitter)")
    group_mRNAprocessing.add_argument("--mRNA_coordinates", action="store_true", help="output a TSV with the coordinates of the surviving features (exons + retained introns)")

    group_seq_tech.add_argument("--seq_tech", choices=["longread", "shortread"], default="shortread", help="sequencing technology")
    group_seq_tech.add_argument("--seq_type", choices=["RNA", "cDNA", "SE", "PE"], default="SE",help="sequencing type")
    group_seq_tech.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')
    group_seq_tech.add_argument("--insert_size", type=str, default="200,300", help="length of insert")
    group_seq_tech.add_argument("--read_length", type=int, default=100, help="length of each read")
    group_seq_tech.add_argument("--s", choices=["rf", "fr", "unstranded"], default="rf", help="library strandedness")
    group_seq_tech.add_argument("--fragments", action="store_true", help="export a ground truth for fragmentation & size selection")
    group_seq_tech.add_argument("--sizeselectiontype", choices=["none", "hardcut", "probabilistic"], default="probabilistic", help="type of size selection. 'hardcut' strictly filters fragments outside the --insert_size limits, while 'probabilistic' applies a double-sided sigmoid retention curve")
    group_seq_tech.add_argument("--no_fragmentation", action="store_true", help="if specified, do not fragment (only valid for metaboliclabeling and ttseq)")

    #TWO-STEP PARSING FOR PRESETS
    #check the known arguments
    temp_args, _ = parser.parse_known_args()

    #show presets
    if temp_args.show_presets:
        print("\n" + "="*45)
        print("           AVAILABLE SPARK PRESETS")
        print("="*45)
        for preset_name, preset_params in PRESETS.items():
            print(f"\n=> {preset_name}")
            if not preset_params:
                print("   (Uses standard CLI defaults)")
            else:
                for k, v in preset_params.items():
                    print(f"   --{k:<20} {v}")
        print("\n" + "="*45 + "\n")
        sys.exit(0)

    #override defaults if presets are chosen
    if temp_args.preset and temp_args.preset != "default":
        parser.set_defaults(**PRESETS[temp_args.preset])

    args = parser.parse_args()
    
    output_dir = args.o.rstrip("/")
    os.makedirs(output_dir, exist_ok=True)

    #start log file
    log_file = os.path.join(output_dir, "spark_pipeline.log")
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        filemode='a')

    c_log(f"Starting SPARK pipeline in mode: {args.mode}")
    if args.preset != "default":
        c_log(f"Loaded '{args.preset}' preset. Applying predefined parameters...")
    c_log(f"Global Settings: Experiment={args.experiment_type}, SeqTech={args.seq_tech}, Threads={args.threads}, RemoveTemp={args.remove_temp}")

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        logging.info(f"Global seed set to: {args.seed}")
    

    if args.region_size_range is None and args.num_regions is None:
        args.region_size_range = "500,5000"
    if args.region_size_range is not None and args.num_regions is not None:
        raise ValueError("Specify only one of --region_size_range or --num_regions, not both.")

    if args.n != 0 and args.n_gene_clusters > args.n:
        sys.exit(f"Error: --n_gene_clusters ({args.n_gene_clusters}) cannot be greater than -n ({args.n})")
    

    if args.experiment_type == 'netseq':
        logging.warning("nt_inc_prob, subs_prob, or sub_base are not used for NET-seq simulations")
    if args.seq_tech == "longread" and args.seq_type not in ["RNA", "cDNA"]:
        parser.error("When --seq_tech is longread, --seq_type must be either RNA or cDNA.")
    if args.seq_tech == "shortread" and args.seq_type not in ["SE", "PE"]:
        parser.error("When --seq_tech is shortread, --seq_type must be either SE or PE.")
    if not (0.0 <= args.bkg_molecules <= 1.0):
        sys.exit("Error: --bkg_molecules must be a float between 0 and 1")

    try:
        from_bases, _conv_map = parse_conversion_map(args.labeling_base, args.sub_base)
    except ValueError as e:
        parser.error(str(e))
    args.labeling_base = ",".join(from_bases)
    if args.sub_base:
        args.sub_base = ",".join(_conv_map[b] for b in from_bases)

    if args.mode in ["fullpipeline", "gene_selection"]:
        c_log("Running Gene Selection & Clustering...")
        logging.info(f"Settings: GTF={args.gtf}, n_genes={args.n}, clusters={args.n_gene_clusters}")
        cmd_seq = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'seq_and_clustering.R')}",
            f"-g {args.gtf}", f"-f {args.genome_fasta}", f"-o {output_dir}/",
            f"-t {args.threads}", f"-c {args.n_gene_clusters}", f"-n {args.n}"
        ]
        if args.protein_coding_only: cmd_seq.append("--protein_coding_only")
        if args.seed is not None: cmd_seq.extend(['--seed', str(args.seed)])
        run_cmd(" ".join(cmd_seq))
        if args.remove_temp:
            remove_files([os.path.join(output_dir, "metrics_per_transcript_df.tsv"),
                          os.path.join(output_dir, "temp", "longest_transcript_and_features.bed"),
                          os.path.join(output_dir, "temp", "nucleotide_metrics.tsv")])

    if args.mode in ["fullpipeline", "rates_per_region"]:
        c_log("Generating Elongation Rates per Region...")
        logging.info(f"Settings: elong_rate_range={args.elong_rate_range}, pause_time={args.pause_time}, flat_rates={args.flat_rates}")
        input_dir = f"{output_dir}/gtf"
        tsv_files = sorted(glob.glob(f"{input_dir}/*.tsv.gz"))
        if not tsv_files: sys.exit(f"No TSV files found in {input_dir}")
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(process_rates_file, tsv_files), total=len(tsv_files), desc="Rates per region"))

    if args.mode in ["fullpipeline", "mRNAgeneration"]:
        c_log("Generating Pre-mRNAs...")
        logging.info(f"Settings: experiment_time={args.experiment_time}, mRNAcount={args.mRNAcount}, bkg_molecules={args.bkg_molecules}")
        dir_region_files=f"{output_dir}/rate_per_gene/"
        region_files = sorted(glob.glob(os.path.join(dir_region_files, "*VariableElongationRateRegions.gtf")))
        premRNA_out_dir = os.path.join(output_dir, "pre-mRNA")
        os.makedirs(premRNA_out_dir, exist_ok=True)
        
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(process_mrna_file, region_files), total=len(region_files), desc="Pre-mRNA Gen"))

    if args.mode in ["fullpipeline", "mRNAgeneration"]:
        premRNA_out_dir = os.path.join(output_dir, "pre-mRNA")
        
        if args.experiment_type == 'chromatin_associated':
            c_log("Running Chromatin Isolation...")
            transcript_files_for_isolation = sorted(glob.glob(os.path.join(premRNA_out_dir, "*.tsv.gz")))
            with worker_pool() as pool:
                list(tqdm(pool.imap_unordered(process_chromatin_isolation, transcript_files_for_isolation), total=len(transcript_files_for_isolation), desc="Chromatin Iso"))

            splice_input_dir = os.path.join(output_dir, "caRNA_filtered")
        else:
            splice_input_dir = premRNA_out_dir
            
        c_log("Processing Splicing Events...")
        logging.info(f"Settings: intron_half_life={args.intron_half_life}, nosplicing={args.nosplicing}")
        transcript_files_for_splicing = sorted(glob.glob(os.path.join(splice_input_dir, "*.tsv.gz")))
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(process_splicing, transcript_files_for_splicing), total=len(transcript_files_for_splicing), desc="Splicing"))
        if args.remove_temp:
            sweep_intermediates(output_dir, [("pre-mRNA", "*.tsv.gz"), ("caRNA_filtered", "*.tsv.gz")])

    if args.mode in ["fullpipeline", "readgeneration"]:
        today_str = datetime.today().strftime("%Y%m%d")
        reads_out_dir = os.path.join(output_dir, "reads")
        libraries_out_dir = os.path.join(output_dir, "final_libraries")
        os.makedirs(reads_out_dir, exist_ok=True)
        os.makedirs(libraries_out_dir, exist_ok=True)
        
        c_log(f"Generating Sequences ({args.experiment_type})...")
        logging.info(f"Settings: seq_type={args.seq_type}, seq_depth={args.seq_depth}, insert_size={args.insert_size}")

    if args.experiment_type == 'proseq' and args.mode in ["fullpipeline", "readgeneration"]:
        proseq_mrnas_dirout = os.path.join(output_dir, "proseqmRNAs")
        os.makedirs(proseq_mrnas_dirout, exist_ok=True)
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(mrna_to_proseq, tsv_files), total=len(tsv_files), desc="mRNA -> PROseq"))
        
        proseq_files = [f for f in sorted(glob.glob(os.path.join(proseq_mrnas_dirout, "*_proseq.tsv.gz"))) if "background" not in os.path.basename(f)]
        generate_reads(proseq_to_fastq, proseq_files, "_proseq.tsv.gz", "PROseq -> FastQ", rebuild_pool=mrna_to_proseq)

    elif args.experiment_type == 'netseq' and args.mode in ["fullpipeline", "readgeneration"]:
        netseq_mrnas_dirout = os.path.join(output_dir, "netseqmRNAs")
        os.makedirs(netseq_mrnas_dirout, exist_ok=True)
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with worker_pool() as pool:
            list(tqdm(pool.imap_unordered(mrna_to_netseq, tsv_files), total=len(tsv_files), desc="mRNA -> NETseq"))

        netseq_files = [f for f in sorted(glob.glob(os.path.join(netseq_mrnas_dirout, "*_netseq.tsv.gz"))) if "background" not in os.path.basename(f)]
        generate_reads(netseq_to_fastq, netseq_files, "_netseq.tsv.gz", "NETseq -> FastQ", rebuild_pool=mrna_to_netseq)

    elif args.experiment_type == 'metaboliclabeling' and args.mode in ["fullpipeline", "readgeneration"]:
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "*.tsv.gz"))) if "background" not in os.path.basename(f)]
        generate_reads(nascent_to_fastq, tsv_files, ".tsv.gz", "NascentRNA -> FastQ")

    elif args.experiment_type == 'ttseq' and args.mode in ["fullpipeline", "readgeneration"]:
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "*.tsv.gz"))) if "background" not in os.path.basename(f)]
        generate_reads(ttseq_to_fastq, tsv_files, ".tsv.gz", "TTseq -> FastQ")
            
    elif args.experiment_type == 'chromatin_associated' and args.mode in ["fullpipeline", "readgeneration"]:
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "*.tsv.gz"))) if "background" not in os.path.basename(f)]
        generate_reads(carna_to_fastq, tsv_files, ".tsv.gz", "caRNA -> FastQ")
            
    if args.mode in ["fullpipeline", "readgeneration"]:
        if not glob.glob(f"{reads_out_dir}/*fastq*"):
            c_log("No read files found to concatenate. Skipping final step.")
            return

        c_log("Concatenating final libraries...")
        final_r1_path = None
        final_r2_path = None

        if args.seq_tech == "longread":
            output_filename = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}.fastq.gz"
            final_r1_path = os.path.join(libraries_out_dir, output_filename)
            run_cmd(f'cat "{reads_out_dir}"/*.fastq.gz > "{final_r1_path}"')
        
        elif args.seq_tech == "shortread" and args.seq_type == "PE":
            output_r1 = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}_R1.fastq.gz"
            output_r2 = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}_R2.fastq.gz"
            final_r1_path = os.path.join(libraries_out_dir, output_r1)
            final_r2_path = os.path.join(libraries_out_dir, output_r2)
            
            run_cmd(f'cat "{reads_out_dir}"/*_R1.fastq.gz > "{final_r1_path}"')
            run_cmd(f'cat "{reads_out_dir}"/*_R2.fastq.gz > "{final_r2_path}"')
        
        else: 
            output_se = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}.fastq.gz"
            final_r1_path = os.path.join(libraries_out_dir, output_se)
            run_cmd(f'cat "{reads_out_dir}"/*.fastq.gz > "{final_r1_path}"')
            
        if final_r1_path and os.path.exists(final_r1_path):
            downsample_library(final_r1_path, final_r2_path, args.seq_depth, args.seed)

    #consolidate ground truth data into a master file
    consolidate_ground_truth(output_dir, args)

    if args.mode in ["fullpipeline", "readgeneration"]:
        c_log("Cleaning up temporary read files...")
        reads_out_dir = os.path.join(output_dir, "reads")
        if os.path.exists(reads_out_dir):
            run_cmd(f'rm -f "{reads_out_dir}"/*fastq*')
        
        temp_dir = f"{output_dir}/temp"
        if os.path.isdir(temp_dir):
            run_cmd(f'rm -rf "{temp_dir}"')

        if args.remove_temp:
            c_log("Removing remaining intermediate files (--remove_temp)...")
            sweep_intermediates(output_dir, [
                ("gtf", "*.tsv.gz"), ("rate_per_gene", "*.tsv"),
                ("pre-mRNA", "*.tsv.gz"), ("caRNA_filtered", "*.tsv.gz"), ("mRNA", "*.tsv.gz"),
                ("proseqmRNAs", "*.tsv.gz"), ("netseqmRNAs", "*.tsv.gz"), ("reads", "*fastq*"),
            ])
            remove_files([os.path.join(output_dir, "metrics_per_transcript_df.tsv")])
            
    c_log("Pipeline run complete!")

if __name__ == "__main__":
    main()




