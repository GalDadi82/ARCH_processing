"""Remove recurrent mutations that lack haematopoietic COSMIC or ClinVar support.

Port of the notebook's recurrence cell. A variant seen in >=2 samples is kept only if it
has haematopoietic COSMIC occurrences or a real ClinVar disease annotation; single-sample
variants always pass. When the step is disabled the listed mutations are carried through
untouched.
"""

import re

import pandas as pd

VARIANT_COLS = ["CHR", "POS", "REF", "ALT", "Gene.refGene", "Selected.AA.change"]
NON_CLINVAR = [".", "not_provided", "not_specified"]

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    df = pd.read_csv(snakemake.input.listed, sep="\t")

    if not snakemake.params.enabled:
        df.to_csv(snakemake.output.filtered, index=False, sep="\t")
        emit("Recurrence filter skipped (filter_recurrent_without_support=false)")
        emit(f"Carrying through all {len(df)} mutations across {df['Sample_Name'].nunique()} samples")
    else:

        def extract_total_cosmic_count(cosmic_str):
            if pd.isna(cosmic_str) or cosmic_str == ".":
                return 0
            occ_match = re.search(r"OCCURENCE=([^;]+)", cosmic_str)
            if not occ_match:
                return 0
            return sum(int(c) for c in re.findall(r"(\d+)\(", occ_match.group(1)))

        # Count the samples each variant appears in.
        variant_sample_counts = (
            df.groupby(VARIANT_COLS)["Sample_Name"].nunique().reset_index(name="Sample_Count")
        )
        df_with_counts = pd.merge(df, variant_sample_counts, on=VARIANT_COLS, how="left")

        non_recurrent = df_with_counts[df_with_counts["Sample_Count"] == 1].copy()
        recurrent = df_with_counts[df_with_counts["Sample_Count"] >= 2].copy()

        recurrent_info = recurrent[VARIANT_COLS + ["cosmic103", "CLNDN"]].drop_duplicates(
            subset=VARIANT_COLS
        )
        recurrent_info["COSMIC_Total_Count"] = recurrent_info["cosmic103"].apply(
            extract_total_cosmic_count
        )
        recurrent_info["COSMIC_in_haematopoetic"] = recurrent_info["cosmic103"].str.contains(
            "haematopoietic_and_lymphoid_tissue", na=False
        )
        recurrent_info["in_clinvar"] = ~recurrent_info["CLNDN"].isin(NON_CLINVAR)
        recurrent_info["Keep"] = (
            (recurrent_info["COSMIC_Total_Count"] > 0)
            & recurrent_info["COSMIC_in_haematopoetic"]
        ) | recurrent_info["in_clinvar"]
        emit(f"Recurrent variant support:\n{recurrent_info['Keep'].value_counts().to_string()}")

        # Force-include specific mutations regardless of the recurrence criteria.
        for spec in snakemake.params.force_include or []:
            mask = pd.Series(True, index=recurrent_info.index)
            for col, val in spec.items():
                if col in recurrent_info.columns:
                    mask &= recurrent_info[col] == val
                else:
                    emit(f"WARNING: force_include column '{col}' not in the table; ignoring it")
            if mask.any():
                emit(f"Force-including {int(mask.sum())} variant(s) matching {spec}")
            else:
                emit(f"WARNING: force_include entry matched nothing: {spec}")
            recurrent_info.loc[mask, "Keep"] = True

        keep_variants = recurrent_info[recurrent_info["Keep"]][VARIANT_COLS]
        filtered_recurrent = pd.merge(recurrent, keep_variants, on=VARIANT_COLS, how="inner")

        final_df = pd.concat([non_recurrent, filtered_recurrent], ignore_index=True)
        final_df = final_df.drop(columns=["Sample_Count"])
        final_df.to_csv(snakemake.output.filtered, index=False, sep="\t")

        emit(
            f"Removed {len(df) - len(final_df)} recurrent variants without haematopoietic "
            "COSMIC or ClinVar support"
        )
        emit(
            f"Remaining: {len(final_df)} mutations across "
            f"{final_df['Sample_Name'].nunique()} samples"
        )
