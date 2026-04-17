import pandas as pd
import ast
import random
import argparse
import os
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to PRO-seq compatible coordinates")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()
    if args.seed is not None:
        np.random.seed(args.seed)

    df = pd.read_csv(args.input_df, sep="\t", comment="#")

    # If incorporated_positions is stored as string, convert to list safely
    df["incorporated_positions"] = df["incorporated_positions"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else (x if isinstance(x, list) else [])
    )

    # Keep only rows where an analog was actually incorporated (Polymerase stalled)
    df = df[df["incorporated_positions"].apply(lambda x: len(x) > 0)].copy()

    # The 3' end of the nascent RNA is precisely where the first analog was incorporated.
    # We set this as the 'sequence_length' so the downstream R chopper can fragment 
    # backwards mathematically from this exact physical coordinate.
    df["sequence_length"] = df["incorporated_positions"].apply(min)

    # Generate the output filename
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_proseq.tsv.gz")

    # Define columns to keep, dropping any legacy string variables
    columns_to_keep = [
        "initiation_time",
        "molecule_id",
        "strand",
        "sub_rate_percent_range",
        "start_label_pos",
        "stop_label_pos",
        "converted_positions",
        "incorporated_positions",
        "percentage_seq_err",
        "seq_err_positions",
        "sequence_length"
    ]
    
    # Ensure we only export columns that exist (in case of background vs main pipeline differences)
    final_cols = [c for c in columns_to_keep if c in df.columns]

    # Convert to DataFrame and export
    df_for_export = df[final_cols]
    df_for_export.to_csv(output_filename, sep="\t", index=False, compression="gzip")
