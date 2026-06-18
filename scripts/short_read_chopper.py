#!/usr/bin/env python3
"""
this was written to be able to chop acorss seq methods, such that:

tree output modes (--mode) hardcoded into the fastq generators that call thiss cript:
  normal  – used by fastq_generator_short_read.py
            Collapsed output (transcript_id → "start-end,start-end,..."), writes TPM,
            optionally writes coverage ground-truth files.
  proseq  – used by fastq_generator_proseq.py
            Per-row output (transcript_id, molecule_id, read_start, read_end).
  ttseq   – used by mRNA_to_reads_ttseq.py
            Collapsed output, no TPM, drops fragments shorter than 10 nt.

"""

import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import gamma as sp_gamma



def _apply_size_selection(lengths, ins_min, ins_max, sizeselectiontype, rng):
    """return a boolean to know which fragments pass size selection."""
    if sizeselectiontype == "none":
        return lengths >= 1
    elif sizeselectiontype == "probabilistic":
        p_left  = 1.0 / (1.0 + np.exp(-0.1 * (lengths - ins_min)))
        p_right = 1.0 / (1.0 + np.exp( 0.1 * (lengths - ins_max)))
        rolls   = rng.uniform(0.0, 1.0, len(lengths))
        return rolls <= (p_left * p_right)
    else:  # hardcut
        return (lengths >= ins_min) & (lengths <= ins_max)


def _fragment_molecule_raw(seq_len, ins_min, ins_max, rng):
    eta   = (ins_min + ins_max) / 2.0
    delta = np.log10(max(seq_len, 2))          # log10 must be > 0

    #number of internal cut-points
    n_cuts = max(0, round(seq_len / eta / sp_gamma(1.0 / delta + 1)) - 1)

    if n_cuts > 0:
        pts  = np.sort(rng.uniform(0.0, 1.0, n_cuts))
        xis  = np.diff(np.concatenate([[0.0], pts, [1.0]]))
    else:
        xis = np.array([1.0])

    xis_t = xis ** (1.0 / delta)
    s     = xis_t.sum()
    if s <= 0:
        xis_t = np.ones_like(xis)
        s     = float(len(xis))

    d       = np.round(seq_len * xis_t / s).astype(int)
    d       = np.maximum(d, 1)                #I don't want zero-size fragments
    n       = len(d)
    cumsums = np.cumsum(d)                    #length n


    if n > 1:
        s0_max     = max(1, min(ins_min, int(d[0])))
        random_s   = int(rng.integers(1, s0_max + 1))

        last_gap   = max(1, int(cumsums[-1]) - int(d[-1]))
        last_max   = max(1, min(ins_min, last_gap))
        random_off = int(rng.integers(1, last_max + 1))

        #starts: [random_s, cumsum[0], cumsum[1], ..., cumsum[n-2]]
        int_starts = np.concatenate([[random_s], cumsums[:-1]])
        #ends:   [cumsum[0], ..., cumsum[n-2], sum(d) - random_off]
        int_ends   = np.concatenate([cumsums[:-1], [int(cumsums[-1]) - random_off]])
    else:
        s0_max     = max(1, min(ins_min, int(d[0])))
        random_s   = int(rng.integers(1, s0_max + 1))
        int_starts = np.array([random_s])
        int_ends   = np.array([int(d[0])])

    all_starts = np.concatenate([[1],                int_starts, [int(int_ends[-1])]])
    all_ends   = np.concatenate([[int(int_starts[0])], int_ends,   [seq_len]])

    #filter out zero length entries
    keep       = all_ends > all_starts
    return all_starts[keep], all_ends[keep]


def fragment_molecule(seq_len, ins_min, ins_max, rng,
                      no_fragmentation=False, sizeselectiontype="probabilistic"):
    """
    fragment single molecules and size select
    """
    if seq_len < 1:
        empty = np.array([], dtype=int)
        return empty, empty, empty, empty

    if no_fragmentation:
        s = np.array([1],       dtype=int)
        e = np.array([seq_len], dtype=int)
        lengths = e - s
        passed  = _apply_size_selection(lengths, ins_min, ins_max, sizeselectiontype, rng)
        return s[passed], e[passed], s, e

    raw_s, raw_e = _fragment_molecule_raw(seq_len, ins_min, ins_max, rng)
    lengths      = raw_e - raw_s
    passed       = _apply_size_selection(lengths, ins_min, ins_max, sizeselectiontype, rng)
    return raw_s[passed], raw_e[passed], raw_s, raw_e

def compute_coverage(starts, ends):
    if len(starts) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    max_pos = int(ends.max())
    cov = np.zeros(max_pos + 2, dtype=np.int64)
    np.add.at(cov, starts.astype(int),     1)
    np.add.at(cov, ends.astype(int),      -1)
    cov = np.cumsum(cov)[:max_pos]
    pos = np.arange(1, max_pos + 1)
    mask = cov > 0
    return pos[mask], cov[mask]




def main():
    parser = argparse.ArgumentParser(
        description="Short-read chopper – Python port of the three R chopper scripts"
    )
    parser.add_argument("--tsv",             required=True,
                        help="Input TSV.gz with simulated mRNAs (must have sequence_length column)")
    parser.add_argument("--mode",            required=True,
                        choices=["normal", "proseq", "ttseq"],
                        help="Output mode: normal | proseq | ttseq")
    parser.add_argument("--insert_size",     default="200,300",
                        help="Insert-size range min,max  e.g. 200,300")
    parser.add_argument("--read_length",     type=int, default=100,
                        help="Read length (kept for CLI compatibility; not used in fragmentation)")
    parser.add_argument("--threads",         type=int, default=1,
                        help="Threads (kept for CLI compatibility)")
    parser.add_argument("--seq_depth",       type=int, default=20_000_000,
                        help="Total library sequencing depth (used for TPM in normal mode)")
    parser.add_argument("--tpm_lower_limit", type=int, default=5)
    parser.add_argument("--tpm_upper_limit", type=int, default=200)
    parser.add_argument("-o", "--dir_out",   default=".",
                        help="Output directory (same -o as fastq generators)")
    parser.add_argument("--seed",            type=int, default=None)
    parser.add_argument("--fragments",       type=str, default="no",
                        help="Pass 'with_ground_truth' to write coverage files (normal mode only)")
    parser.add_argument("--sizeselectiontype", default="probabilistic",
                        choices=["none", "hardcut", "probabilistic"])
    parser.add_argument("--no_fragmentation", action="store_true",
                        help="Skip fragmentation; treat each molecule as a single fragment")

    args = parser.parse_args()

    ins_min, ins_max = map(int, args.insert_size.split(","))


    gene_id_raw = os.path.basename(args.tsv)
    gene_id     = gene_id_raw
    for ext in (".tsv.gz", ".gz", ".tsv"):
        if gene_id.endswith(ext):
            gene_id = gene_id[: -len(ext)]
            break

    gene_id_for_hash = gene_id.split("_")[0]

    rng = np.random.default_rng()
    if args.seed is not None:
        h           = hashlib.sha256(gene_id_for_hash.encode("utf-8"))
        gene_hash   = int(h.hexdigest(), 16) % (2 ** 32)
        unique_seed = (args.seed + gene_hash) % (2 ** 32)
        rng         = np.random.default_rng(unique_seed)


    if args.mode == "proseq":
        base_part  = gene_id.split("_")[0]
        out_gene_id = (base_part + "_background") if "_background" in gene_id else base_part
    else:
        out_gene_id = gene_id

    try:
        df = pd.read_csv(args.tsv, sep='\t')
        
        if df.empty:
            print(f"Warning: {args.tsv} contains no molecules (likely 100% background). Exiting gracefully.")
            sys.exit(0)
            
        if 'sequence_length' not in df.columns:
            raise ValueError(f"'sequence_length' column not found in {args.tsv}")

    #this is here to account for cases were bkg 1 and nascent mRNA dfs are empty
    except pd.errors.EmptyDataError:
        print(f"Warning: {args.tsv} is completely empty. Ignore this warning if bkg 1 was set")
        sys.exit(0)

    if "full_molecule_sequence" in df.columns:
        df = df.drop(columns=["full_molecule_sequence"])
        
    df["transcript_id"] = range(1, len(df) + 1)
    
    if len(df) <= 5:
        return  


    sel_starts_list = []
    sel_ends_list   = []
    all_starts_list = []
    all_ends_list   = []
    tid_list        = []

    need_pre_cov = (args.mode == "normal" and args.fragments == "with_ground_truth")

    for row in df.itertuples(index=False):
        seq_len = int(row.sequence_length)
        tid     = int(row.transcript_id)

        s_sel, e_sel, s_all, e_all = fragment_molecule(
            seq_len, ins_min, ins_max, rng,
            no_fragmentation=args.no_fragmentation,
            sizeselectiontype=args.sizeselectiontype,
        )

        if len(s_sel) > 0:
            sel_starts_list.append(s_sel)
            sel_ends_list.append(e_sel)
            tid_list.append(np.full(len(s_sel), tid, dtype=int))

        if need_pre_cov and len(s_all) > 0:
            all_starts_list.append(s_all)
            all_ends_list.append(e_all)

    if not sel_starts_list:
        return  # nothing passed size selection

    sel_starts = np.concatenate(sel_starts_list)
    sel_ends   = np.concatenate(sel_ends_list)
    tids       = np.concatenate(tid_list)

    #TTseq: drop fragments shorter than 10nt
    if args.mode == "ttseq":
        keep = (sel_ends - sel_starts) > 10
        sel_starts, sel_ends, tids = sel_starts[keep], sel_ends[keep], tids[keep]
        if len(sel_starts) == 0:
            return

    if len(sel_starts) <= 10:
        return

    dir_out   = args.dir_out.rstrip("/")
    frags_dir = os.path.join(dir_out, "temp", "mRNAs_with_fragments")
    temp_dir  = os.path.join(dir_out, "temp")
    os.makedirs(frags_dir, exist_ok=True)
    os.makedirs(temp_dir,  exist_ok=True)
    out_path = os.path.join(frags_dir, f"{out_gene_id}_fragments.tsv")

    if args.mode == "normal":
        gene_length_kb = float(df["sequence_length"].mean()) / 1000.0
        gene_tpm       = float(rng.uniform(args.tpm_lower_limit, args.tpm_upper_limit))
        seq_depth_m    = args.seq_depth / 1e6
        reads_to_get   = gene_length_kb * gene_tpm * seq_depth_m

        #TPM temp file (may be overwritten by fastq_generator_short_read.py – that's ok)
        tpm_path = os.path.join(temp_dir, f"temp_{out_gene_id}_tpm.tsv")
        pd.DataFrame([{"gene_id": out_gene_id, "tpm": gene_tpm}]).to_csv(
            tpm_path, sep="\t", index=False
        )

        #ground truth if specified
        if args.fragments == "with_ground_truth":
            #oost-size-selection coverage
            pos_sel, cov_sel = compute_coverage(sel_starts, sel_ends)
            if len(pos_sel) > 0:
                gt_sel_dir = os.path.join(dir_out, "ground_truth_after_size_selection")
                os.makedirs(gt_sel_dir, exist_ok=True)
                pd.DataFrame({"position": pos_sel, "frequency": cov_sel}).to_csv(
                    os.path.join(gt_sel_dir, f"{out_gene_id}.tsv.gz"),
                    sep="\t", index=False
                )

            #pre-size-selection coverage (scaled to target depth)
            if all_starts_list:
                all_s = np.concatenate(all_starts_list)
                all_e = np.concatenate(all_ends_list)
                n_passed = len(sel_starts)
                scaling  = (reads_to_get / n_passed) if n_passed > 0 else 0.0
                pos_pre, cov_pre = compute_coverage(all_s, all_e)
                if len(pos_pre) > 0 and scaling > 0:
                    gt_pre_dir = os.path.join(dir_out, "ground_truth_pre_size_selection")
                    os.makedirs(gt_pre_dir, exist_ok=True)
                    pd.DataFrame({
                        "position":  pos_pre,
                        "frequency": (cov_pre * scaling).round().astype(int),
                    }).to_csv(
                        os.path.join(gt_pre_dir, f"{out_gene_id}.tsv.gz"),
                        sep="\t", index=False
                    )

        
        frag_df = pd.DataFrame({
            "transcript_id": tids.astype(int),
            "read_start":    sel_starts.astype(int),
            "read_end":      sel_ends.astype(int),
        })
        collapsed = (
            frag_df
            .assign(coord=frag_df["read_start"].astype(str) + "-" + frag_df["read_end"].astype(str))
            .groupby("transcript_id", sort=False)["coord"]
            .agg(",".join)
            .reset_index()
            .rename(columns={"coord": "read_coordinates"})
        )
        collapsed.to_csv(out_path, sep="\t", index=False)

    elif args.mode == "proseq":
        if "molecule_id" not in df.columns:
            raise ValueError("'molecule_id' column required for proseq mode")

        # PRO-seq captures the 3' end of the nascent RNA (where Pol II stalled).
        # Each fragment must be anchored at sequence_length, not randomly placed.
        seq_lengths = df["sequence_length"].values.astype(int)
        fragment_sizes = rng.integers(ins_min, ins_max + 1, size=len(df))
        read_ends = seq_lengths
        read_starts = np.maximum(1, seq_lengths - fragment_sizes)

        keep = (read_ends - read_starts) > 10
        frag_df = pd.DataFrame({
            "transcript_id": df["transcript_id"].values[keep].astype(int),
            "molecule_id":   df["molecule_id"].values[keep],
            "read_start":    read_starts[keep].astype(int),
            "read_end":      read_ends[keep].astype(int),
        })
        frag_df.to_csv(out_path, sep="\t", index=False)

    elif args.mode == "ttseq":
        frag_df = pd.DataFrame({
            "transcript_id": tids.astype(int),
            "read_start":    sel_starts.astype(int),
            "read_end":      sel_ends.astype(int),
        })
        collapsed = (
            frag_df
            .assign(coord=frag_df["read_start"].astype(str) + "-" + frag_df["read_end"].astype(str))
            .groupby("transcript_id", sort=False)["coord"]
            .agg(",".join)
            .reset_index()
            .rename(columns={"coord": "read_coordinates"})
        )
        collapsed.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
