"""Apply the final VAF threshold and drop known artefact / germline calls.

Port of the notebook's "VAF filter and remove excluded mutations" cell. Genes listed in
gene_min_vaf_overrides use their own threshold; everything else uses min_final_vaf.
"""

import sys

import pandas as pd

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    min_final_vaf = snakemake.params.min_final_vaf
    overrides = dict(snakemake.params.gene_min_vaf_overrides or {})
    exclusion_keys = list(snakemake.params.exclusion_keys)
    excluded = list(snakemake.params.excluded_mutations or [])

    final_df = pd.read_csv(snakemake.input.mutations, sep="\t")

    effective_min_vaf = final_df["Gene.refGene"].map(overrides).fillna(min_final_vaf)
    vaf_filtered = final_df[final_df["VAF"] > effective_min_vaf].copy()
    emit(
        f"After VAF filter (min_final_vaf={min_final_vaf}, overrides={overrides}): "
        f"{len(vaf_filtered)}"
    )

    if excluded:
        to_remove = pd.DataFrame(excluded)
        missing = [c for c in exclusion_keys if c not in to_remove.columns]
        if missing:
            sys.exit(f"excluded_mutations entries are missing exclusion_keys {missing}")

        # Match on the configured key columns, comparing as strings so that e.g. a POS
        # read as int64 still matches the YAML value.
        def key_index(df):
            return pd.MultiIndex.from_frame(df[exclusion_keys].astype(str))

        mask = ~key_index(vaf_filtered).isin(key_index(to_remove))
        n_before = len(vaf_filtered)
        vaf_filtered = vaf_filtered[mask]
        emit(
            f"After removing {len(to_remove)} excluded mutation(s): {len(vaf_filtered)} "
            f"({n_before - len(vaf_filtered)} row(s) dropped)"
        )

    vaf_filtered.to_csv(snakemake.output.master, index=False, sep="\t")

    emit(f"Removed {len(final_df) - len(vaf_filtered)} rows via VAF filter and manual exclusions")
    emit(
        f"Remaining: {len(vaf_filtered)} mutations across "
        f"{vaf_filtered['Sample_Name'].nunique()} samples"
    )
    emit(f"Saved to: {snakemake.output.master}")
