"""Flatten one ANNOVAR multianno table into the columns downstream filtering expects.

Port of varscan_multianno_to_df_no_filter_fast() from the notebook: keep the annotation
columns, then append FILTER plus Depth and VAF pulled out of the VCF sample field.

An empty ANNOVAR table is not an error -- a sample can legitimately have no variants -- so
in that case a header-only output is written and the pipeline carries on.
"""

import os
import sys

import pandas as pd
from pandas.errors import EmptyDataError

filename = snakemake.input.txt
vcf_file = snakemake.input.vcf
out_filename = snakemake.output.formatted

DEFAULT_VCF_HEADER = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE",
]

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    # Read as all strings for speed. header=None keeps pandas from misaligning columns:
    # ANNOVAR's header row has fewer fields than the data rows.
    try:
        df = pd.read_csv(filename, sep="\t", skiprows=1, header=None, dtype=str)
    except EmptyDataError:
        emit(f"Empty ANNOVAR table: {filename}. Writing header-only output.")
        pd.DataFrame().to_csv(out_filename, sep="\t", index=False)
        sys.exit(0)

    with open(filename) as f:
        txt_header = f.readline().strip("\n").split("\t")

    vcf_header = []
    if os.path.exists(vcf_file):
        with open(vcf_file) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_header = line.strip("\n").split("\t")
                    break
    if not vcf_header:
        emit(f"No #CHROM line in {vcf_file}. Using default VCF headers.")
        vcf_header = list(DEFAULT_VCF_HEADER)

    # Align columns from the left with the ANNOVAR header, pad if ANNOVAR inserted extra
    # columns, and pin the VCF headers to the end.
    num_data_cols = df.shape[1]
    num_vcf_cols = len(vcf_header)

    new_columns = []
    for i in range(num_data_cols - num_vcf_cols):
        new_columns.append(txt_header[i] if i < len(txt_header) else f"Padding_{i}")
    new_columns.extend(vcf_header)
    df.columns = new_columns

    if "Otherinfo1" in df.columns:
        # Keep the original behaviour: drop the column immediately before Otherinfo1.
        end_of_annotation_idx = max(0, list(df.columns).index("Otherinfo1") - 1)
    else:
        end_of_annotation_idx = num_data_cols - num_vcf_cols

    if "FORMAT" not in df.columns:
        sys.exit(f"FORMAT column not found in {filename}")

    format_col = df["FORMAT"].dropna()
    if format_col.empty:
        sys.exit(f"FORMAT column is entirely empty in {filename}")
    format_keys = str(format_col.iloc[0]).split(":")

    sample_expanded = df.iloc[:, -1].str.split(":", expand=True)

    def field(key, name):
        """Pull one FORMAT field out of the expanded sample column."""
        if key in format_keys:
            idx = format_keys.index(key)
            if idx in sample_expanded.columns:
                return sample_expanded[idx].rename(name)
        emit(f"FORMAT field {key} unavailable; filling {name} with '.'")
        return pd.Series(["."] * len(df), name=name, index=df.index)

    depth_col = field("DP", "Depth")
    vaf_col = field("FREQ", "VAF")
    filter_col = (
        df["FILTER"] if "FILTER" in df.columns
        else pd.Series(["."] * len(df), name="FILTER", index=df.index)
    )

    final_df = pd.concat(
        [df.iloc[:, :end_of_annotation_idx], filter_col, depth_col, vaf_col], axis=1
    )
    final_df.to_csv(out_filename, sep="\t", index=False)
    emit(f"Processed {os.path.basename(filename)} ({len(final_df)} rows) -> {out_filename}")
