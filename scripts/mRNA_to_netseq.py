import pandas as pd
import random
import argparse
import os
import numpy as np
import gzip
import hashlib

def chop_sequence_bg_coords(seq_len):
    """
    Determines coordinates for a random fragment (35-100nt) from a background sequence.
    Returns: (start_pos, end_pos)
    """
    fragment_length = random.randint(35, 100)
    
    
    if seq_len <= fragment_length:
        return 0, seq_len

    max_start_pos = seq_len - fragment_length
    start_pos = random.randint(0, max_start_pos)
    end_pos = start_pos + fragment_length
    return start_pos, end_pos

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert mRNAs to NET-seq fragments")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")
    parser.add_argument("--bkg_molecules", type=float, default=0.0, 
                        help="Proportion of total reads that should be background (0.0 to 1.0). Default 0.")

    args = parser.parse_args()

    
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_netseq.tsv.gz")

    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % 4294967295
        np.random.seed(unique_seed)
        random.seed(unique_seed)

    
    path_to_rates = os.path.join(os.path.dirname(args.o), "rate_per_gene", base_filename + "_RatesandTraversalTimes.gtf")
    
    
    if os.path.exists(path_to_rates):
        rates_df = pd.read_csv(path_to_rates, sep="\t", comment="#")
        rates_df["nucleotide_coord"] = rates_df["nucleotide_coord"].astype(int)
        max_length = rates_df["nucleotide_coord"].max()
    else:
        
        max_length = float('inf')

    
    output_columns = [
        "initiation_time", "molecule_id", "strand", "sub_rate_percent_range",
        "fragment_start_pos", "fragment_end_pos"
    ]

    with gzip.open(output_filename, 'wt') as f_out:
        f_out.write("\t".join(output_columns) + "\n")

        print(f"[NET-seq] Processing main gene coordinates: {base_filename}")
        
        df = pd.read_csv(args.input_df, sep="\t", comment="#")
        

        random_offsets = np.random.randint(35, 101, size=len(df))

        
        df["fragment_start_pos"] = (df["stop_label_pos"] - random_offsets).clip(lower=0)
        df["fragment_end_pos"] = df["stop_label_pos"]
        
        # size filter to avoid short inserts
        df["frag_len"] = df["fragment_end_pos"] - df["fragment_start_pos"]
        df_valid = df[df["frag_len"] > 10]
        
        for row in df_valid.itertuples():
            f_out.write(
                f"{row.initiation_time}\t"
                f"{row.molecule_id}\t"
                f"{row.strand}\t"
                f"{row.sub_rate_percent_range}\t"
                f"{row.fragment_start_pos}\t"
                f"{row.fragment_end_pos}\n"
            )
        
        del df
        del df_valid

        # bkg molecules

        if args.bkg_molecules > 0.0:
            path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
            
            if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
                df_bg = pd.read_csv(path_to_BGmRNAs, sep="\t", comment="#")

                if not df_bg.empty:
                    print(f"[NET-seq] Adding {len(df_bg)} background coordinate footprints...")

                    for row in df_bg.itertuples():
                        mol_id = str(row.molecule_id) + '_BG_'
                        
                        # Use the sequence_length column to determine bounds mathematically
                        seq_len = getattr(row, 'sequence_length', 1000) 
                        fstart, fend = chop_sequence_bg_coords(seq_len)
                        
                        f_out.write(
                            f"{row.initiation_time}\t"
                            f"{mol_id}\t"
                            f"{row.strand}\t"
                            f"{row.sub_rate_percent_range}\t"
                            f"{fstart}\t"
                            f"{fend}\n"
                        )
                    del df_bg
