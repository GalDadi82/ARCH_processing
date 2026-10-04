"""Merge every sample's formatted ANNOVAR table into one mutation table.

Port of the notebook's merge cell: insert Sample_Name as the 6th column and rename the
ANNOVAR coordinate columns to the names the Blood-2023 filtering expects.
"""

import os
import sys

import pandas as pd

RENAMES = {"Chr": "CHR", "Start": "POS", "Ref": "REF", "Alt": "ALT"}

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    formatted_files = list(snakemake.input.formatted)
    emit(f"Merging {len(formatted_files)} formatted files")

    merged_path = snakemake.output.merged
    header_written = False
    n_failed = 0

    for idx, file_path in enumerate(formatted_files):
        sample_name = os.path.basename(os.path.dirname(file_path))
        try:
            df_sample = pd.read_csv(file_path, sep="\t", dtype=str)
            if df_sample.empty:
                emit(f"Skipping {sample_name}: no variants")
                continue

            df_sample.insert(5, "Sample_Name", sample_name)
            df_sample = df_sample.rename(columns=RENAMES)

            df_sample.to_csv(
                merged_path,
                sep="\t",
                index=False,
                header=not header_written,
                mode="w" if not header_written else "a",
            )
            header_written = True
            emit(f"Finished {idx}: {sample_name} ({len(df_sample)} rows)")
        except Exception as exc:  # keep going; one bad sample should not sink the run
            n_failed += 1
            emit(f"ERROR in {idx}: {sample_name} ({file_path}): {exc}")

    if not header_written:
        sys.exit("No sample contributed any rows -- merged file would be empty")

    if n_failed:
        emit(f"WARNING: {n_failed} sample(s) failed to merge")
    emit(f"Merged file saved to: {merged_path}")
