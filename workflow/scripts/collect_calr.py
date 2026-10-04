"""Collect per-sample CALR Type-1 / Type-2 tables into the two run-level result files.

Port of the notebook's "Collect CALR results" cell.
"""

import os
import sys

import pandas as pd

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    def collect(files, out_path, label):
        frames = []
        for path in files:
            if not os.path.exists(path):
                emit(f"Alert: {path} does not exist")
                continue
            df = pd.read_csv(path, sep="\t")
            # The util scripts write Sample_Name themselves in some versions; set it from
            # the file name either way so the column is always present and correct.
            df["Sample_Name"] = os.path.basename(path).split("_CALR")[0]
            frames.append(df)

        if not frames:
            sys.exit(f"No {label} results found -- cannot write {out_path}")

        merged = pd.concat(frames, ignore_index=True)
        merged.to_csv(out_path, sep="\t", index=False)
        emit(f"{label}: collected {len(frames)} sample(s), {len(merged)} rows -> {out_path}")

    collect(snakemake.input.type1, snakemake.output.type1, "CALR type-1")
    collect(snakemake.input.type2, snakemake.output.type2, "CALR type-2")
