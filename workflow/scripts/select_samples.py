"""Select the samples to process from the QC classification table.

Port of the tail of the notebook's QC cell: keep samples whose dedup on-target coverage
fraction clears the (lower) processing cutoff, so poor samples are still processed to get
their info.
"""

import sys

import pandas as pd

min_base_cov = snakemake.params.min_base_cov
min_f_to_process = snakemake.params.min_f_dedup_cov_to_process
cov_col = f"base >={min_base_cov} dedup %on-target cov"

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    all_samples = pd.read_csv(snakemake.input.classified)

    if cov_col not in all_samples.columns:
        sys.exit(
            f"Column '{cov_col}' missing from {snakemake.input.classified}. "
            f"Available: {list(all_samples.columns)}"
        )

    selected = all_samples[all_samples[cov_col] >= min_f_to_process]
    selected.to_csv(snakemake.output[0], index=False)

    emit(f"{len(selected)}/{len(all_samples)} samples have {cov_col} >= {min_f_to_process}")
    if selected.empty:
        sys.exit("No samples passed the processing cutoff -- check min_f_dedup_cov_to_process")
    emit(f"Wrote {snakemake.output[0]}")
