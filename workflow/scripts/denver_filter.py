"""Cross-reference the master mutation table against the Denver normals cohort.

Port of the notebook's "Filter by Denver normals" cell. A mutation is valid if it was never
seen in Denver, or was seen and is on Denver's list of valid SNVs or indels; the final
is_valid also keeps anything with a narrow listed match.
"""

import pandas as pd

KEY_COLS = ["CHR", "POS", "REF", "ALT"]

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    vaf_filtered = pd.read_csv(snakemake.input.master, sep="\t")

    denver_df = pd.read_csv(snakemake.input.is_in_denver, sep="\t")
    denver_df = denver_df.rename(columns={"CHROM": "CHR"})

    def add_flag(df, path, flag):
        """Left-join a Denver pass-list, marking matched variants with `flag`."""
        passed = pd.read_csv(path, sep="\t")[KEY_COLS].drop_duplicates()
        passed[flag] = True
        out = df.merge(passed, how="left", on=KEY_COLS)
        # astype(bool) keeps the column a real boolean; a bare fillna leaves object dtype.
        out[flag] = out[flag].fillna(False).astype(bool)
        emit(f"{flag}: {int(out[flag].sum())} of {len(out)} matched {path}")
        return out

    denver_df = add_flag(denver_df, snakemake.input.filtered_snvs, "Valid SNV")
    denver_df = add_flag(denver_df, snakemake.input.filtered_indels, "Valid Indel")

    valid_by_denver = (
        (denver_df["match"] == "none") | denver_df["Valid SNV"] | denver_df["Valid Indel"]
    )
    denver_df["Valid by Denver"] = valid_by_denver
    emit(
        f"{int(valid_by_denver.sum())}/{len(denver_df)} mutations are either not in Denver "
        "or are valid SNVs/Indels in Denver"
    )

    vaf_filtered_with_denver = vaf_filtered.merge(
        denver_df, how="left", on=KEY_COLS, suffixes=("", "_denver")
    )
    vaf_filtered_with_denver["is_valid"] = (
        vaf_filtered_with_denver["Valid by Denver"]
        | (vaf_filtered_with_denver["Has.listed.match_type"] == "narrow")
    )
    vaf_filtered_with_denver.to_csv(snakemake.output.with_denver, index=False, sep="\t")

    denver_filtered = vaf_filtered_with_denver[vaf_filtered_with_denver["is_valid"]]
    denver_filtered.to_csv(snakemake.output.filtered, index=False, sep="\t")

    emit(
        f"Kept {len(denver_filtered)}/{len(vaf_filtered_with_denver)} mutations across "
        f"{denver_filtered['Sample_Name'].nunique()} samples"
    )
    emit(f"Saved to: {snakemake.output.filtered}")
