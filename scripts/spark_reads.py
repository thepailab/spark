import os
import numpy as np
import pandas as pd


def pool_paths(output_dir, gene_id):
    pool_dir = os.path.join(output_dir, "temp", "pools")
    os.makedirs(pool_dir, exist_ok=True)
    return (os.path.join(pool_dir, f"{gene_id}_nascent.tsv.gz"),
            os.path.join(pool_dir, f"{gene_id}_background.tsv.gz"),
            os.path.join(pool_dir, f"{gene_id}_units.tsv"))


def write_units(units_path, gene_id, nascent_units, bg_units, bg_molecules):
    pd.DataFrame([{
        "gene_id": gene_id,
        "nascent_units": int(nascent_units),
        "bg_units": int(bg_units),
        "bg_molecules": int(bg_molecules),
    }]).to_csv(units_path, sep="\t", index=False)


def read_pool(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            return pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            pass
    return pd.DataFrame()


def remove_paths(paths):
    for path in paths:
        try: os.remove(path)
        except OSError: pass


def sample_rows(df, n):
    if n <= 0 or df.empty:
        return df.iloc[0:0]
    if len(df) >= n:
        return df.sample(n=n, replace=False)
    return pd.concat([df, df.sample(n=n - len(df), replace=True)], ignore_index=True)


def count_rows(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            return len(pd.read_csv(path, sep="\t", usecols=[0]))
        except pd.errors.EmptyDataError:
            pass
    return 0


def explode_fragments(path):
    columns = ["transcript_id", "read_start", "read_end"]
    df = read_pool(path)
    if df.empty or "read_coordinates" not in df.columns:
        return pd.DataFrame(columns=columns)
    df["read_coordinates"] = df["read_coordinates"].astype(str).str.split(",")
    df = df.explode("read_coordinates").reset_index(drop=True)
    coords = df["read_coordinates"].str.split("-", expand=True)
    df["read_start"] = pd.to_numeric(coords[0], errors="coerce")
    df["read_end"] = pd.to_numeric(coords[1], errors="coerce")
    df["transcript_id"] = pd.to_numeric(df["transcript_id"], errors="coerce")
    df = df.dropna(subset=columns)
    return df[columns].astype(int)


def gene_references(gtf_path, bg_path):
    df_gtf = pd.read_csv(gtf_path, sep="\t").sort_values("position")
    ref_arr = np.array(list("".join(df_gtf["sequence"].astype(str).tolist())), dtype="U1")
    gene_rows = df_gtf[df_gtf["feature"] != "readthrough"]
    unspliced_length = int(gene_rows["sequence"].astype(str).str.len().sum())

    bg_unspliced = False
    if os.path.exists(bg_path) and os.path.getsize(bg_path) > 0:
        try:
            peek_df = pd.read_csv(bg_path, sep="\t", nrows=1)
            bg_unspliced = ("sequence_length" in peek_df.columns and not peek_df.empty
                            and int(peek_df.iloc[0]["sequence_length"]) == unspliced_length)
        except (pd.errors.EmptyDataError, ValueError):
            pass

    bg_rows = gene_rows if bg_unspliced else df_gtf[df_gtf["feature"] == "exon"]
    bg_ref_arr = np.array(list("".join(bg_rows["sequence"].astype(str).tolist())), dtype="U1")
    return ref_arr, bg_ref_arr


def coordinate_mapper(rates_df, extrapolate=False):
    if rates_df.empty:
        return lambda read_start, read_end: "NA"
    nucleotide_coords = rates_df["nucleotide_coord"].to_numpy()
    anchor = int(np.argmin(nucleotide_coords))
    first_nt, last_nt = int(nucleotide_coords[anchor]), int(nucleotide_coords.max())
    first_position = int(rates_df["absolute_position"].iloc[anchor])
    chrom = rates_df["chromosome"].iloc[0]
    strand = rates_df["strand"].iloc[0]
    direction = 1 if strand == "+" else -1

    def to_absolute(read_start, read_end):
        start, end = int(read_start), int(read_end)
        if start <= end and start <= last_nt and end >= first_nt:
            start, end = max(start, first_nt), min(end, last_nt)
        elif not extrapolate:
            return "NA"
        return f"{chrom}:{first_position + direction * (start - first_nt)}-{first_position + direction * (end - first_nt)}:{strand}"

    return to_absolute
