"""
Filter mutations from multiple samples (annovar output) by VAF, mutation type (e.g. stop
gain/loss, splicing), coordinate within specific genes and a curated list of mutations
(hotspots). Python port of Mutation_Curation/R/FilterMutations.R + the active
annotate_mutations_df() in Mutation_Curation/R/mut_utils.R.

Usage:
  python filter_mutations.py --input-file /path/to/merged_varscan_annotations.txt \
      --output-dir /path/to/5_vcfs --run-name MiseqR14_arch --sample-column Sample_Name
"""

import argparse
import os

import numpy as np
import pandas as pd

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_LISTED_MUTATIONS_FN = os.path.join(
    MODULE_DIR, "..", "Mutation_Curation", "data", "Vlasschaert_Blood_2023_data.xlsx"
)
DEFAULT_MANE_IFN = os.path.join(
    MODULE_DIR, "..", "Mutation_Curation", "data", "MANE.GRCh38.v1.4.summary.txt.gz"
)

# Matches mut_dict in mut_utils.R:300-305, including its underscore/space mismatch against the
# 'nonframeshift deletion'/'nonframeshift insertion' values actually present in the Include sheet
# (kept for parity with the R behavior — those rows silently get no mutation_func mapping there).
MUT_DICT = {
    "frameshift": "frameshift insertion;frameshift deletion",
    "nonsense": "stopgain;stoploss",
    "splicing": "splicing",
    "missense": "nonsynonymous SNV",
    "nonframeshift_deletion": "nonframeshift_deletion",
    "nonframeshift_insertion": "nonframeshift_insertion",
}


def _explode_sep(df, col, sep):
    df = df.copy()
    df[col] = df[col].apply(lambda v: str(v).split(sep) if pd.notna(v) else [None])
    return df.explode(col, ignore_index=True)


def _fmt(v):
    if pd.isna(v):
        return "NA"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def annotate_mutations_df(
    input_df,
    listed_mutations_fn=DEFAULT_LISTED_MUTATIONS_FN,
    mane_ifn=DEFAULT_MANE_IFN,
    mane_status_priorities=("MANE Plus Clinical", "MANE Select", "None"),
    ignore_transcript_id_version=True,
    mut_info_col="AAChange.refGene",
    non_exonic_mut_info_col="GeneDetail.refGene",
    merged_mut_info_col=None,
    func_col="Func.refGene",
    exonic_func_col="ExonicFunc.refGene",
    mut_info_sep=",",
    non_exonic_mut_info_sep=";",
    include_non_syn_by_loci=True,
):
    """Port of the active annotate_mutations_df() in Mutation_Curation/R/mut_utils.R:214-396.

    Adds Has.protein.coding, Top.MANE.status, Has.listed, Has.listed.match_type, Has.excluded,
    Selected.transcript, Selected.AA.change, Last.AA.Change and f"{mut_info_col}.Annotated" to
    input_df.

    Has.listed.match_type distinguishes *why* Has.listed is True: "narrow" means the match came
    from an Include-sheet rule with at least one restriction (exon/AA-range/DNA-range/missense-AA
    -change), i.e. a specific hotspot/domain; "broad" means it only matched an unrestricted rule
    (any occurrence of that mutation type anywhere in the gene qualifies) — e.g. CEBPA's
    frameshift/nonsense/splicing rule, or DNMT3A's LOF rule. Broad matches carry materially lower
    confidence than narrow ones and shouldn't be treated as equivalent evidence. NaN when
    Has.listed is False.
    """
    if merged_mut_info_col is None:
        merged_mut_info_col = f"{mut_info_col}.Annotated"

    if mut_info_col not in input_df.columns:
        raise ValueError(f"Column {mut_info_col} not found in input data frame.")
    if non_exonic_mut_info_col not in input_df.columns:
        raise ValueError(f"Column {non_exonic_mut_info_col} not found in input data frame.")

    input_df = input_df.reset_index(drop=True).copy()
    input_df["Last.AA.Change"] = input_df[mut_info_col].astype(str).str.replace(
        r"^.*:p\.", "", regex=True
    )
    input_df["__row_id__"] = input_df.index

    orig_columns = [c for c in input_df.columns if c not in (mut_info_col, "__row_id__")]

    non_exonic_fixed = (
        input_df[non_exonic_mut_info_col]
        .astype(str)
        .str.replace(non_exonic_mut_info_sep, mut_info_sep, regex=False)
    )
    transcript_info_col = pd.Series(
        np.where(
            input_df[func_col] == "splicing",
            non_exonic_fixed,
            np.where(input_df[func_col] == "exonic", input_df[mut_info_col].astype(str), "."),
        ),
        index=input_df.index,
    )

    df = input_df.copy()
    df["transcript_info_col"] = transcript_info_col
    df = _explode_sep(df, "transcript_info_col", mut_info_sep)

    is_splicing = df[func_col] == "splicing"
    df.loc[is_splicing, "transcript_info_col"] = (
        "None:" + df.loc[is_splicing, "transcript_info_col"].astype(str) + ":p.A-1A"
    )

    sub_fields = ["sub__gene", "sub__transcript", "sub__exon", "sub__cDNA_change", "sub__protein_change"]
    split_cols = df["transcript_info_col"].astype(str).str.split(":", n=4, expand=True)
    split_cols.columns = sub_fields[: split_cols.shape[1]]
    for col in sub_fields:
        if col not in split_cols.columns:
            split_cols[col] = np.nan
    df = pd.concat([df.reset_index(drop=True), split_cols[sub_fields].reset_index(drop=True)], axis=1)

    df["sub__protein_change"] = df["sub__protein_change"].str.replace(
        r"^p\.", "", regex=True
    )
    # Capture the digit run immediately after the leading AA letter(s), e.g. "W288Cfs*12" -> "288".
    # (Stripping a trailing letter run instead would miss this: frameshift suffixes like "fs*12"
    # end in digits, not letters, so the AA position would never be extracted.)
    aa_pos_str = df["sub__protein_change"].str.extract(r"^[A-Za-z]+(\d+)", expand=False)
    df["sub__aa_position"] = pd.to_numeric(aa_pos_str, errors="coerce")

    df["sub__cDNA_change"] = df["sub__cDNA_change"].str.replace(r"^c\.", "", regex=True)
    cdna_letters_stripped = df["sub__cDNA_change"].str.replace(r"[A-Za-z]", "", regex=True)
    cdna_start_str = cdna_letters_stripped.str.replace(r"_[0-9]+", "", regex=True)
    cdna_end_str = cdna_letters_stripped.str.replace(r"[0-9]+_", "", regex=True)
    df["sub__cDNA_start"] = pd.to_numeric(cdna_start_str, errors="coerce")
    df["sub__cDNA_end"] = pd.to_numeric(cdna_end_str, errors="coerce")

    df["Has.protein.coding"] = df["sub__transcript"].str.contains("NM_", na=False)
    df["mutation_func"] = np.where(
        df[func_col] == "splicing", "splicing", df[exonic_func_col].astype(str)
    )
    df = _explode_sep(df, "mutation_func", ";")
    df["row_num"] = np.arange(len(df))

    if ignore_transcript_id_version:
        df["transcript_for_join"] = df["sub__transcript"].str.replace(
            r"\.[0-9]+$", "", regex=True
        )
    else:
        df["transcript_for_join"] = df["sub__transcript"]

    if mane_ifn and os.path.exists(mane_ifn):
        mane = pd.read_csv(mane_ifn, sep="\t")
        if ignore_transcript_id_version:
            mane["transcript_for_join"] = mane["RefSeq_nuc"].astype(str).str.replace(
                r"\.[0-9]+$", "", regex=True
            )
        else:
            mane["transcript_for_join"] = mane["RefSeq_nuc"]
        df = df.merge(mane[["transcript_for_join", "MANE_status"]], on="transcript_for_join", how="left")
        df["MANE_status"] = df["MANE_status"].fillna("None")
        df["MANE_status"] = pd.Categorical(
            df["MANE_status"], categories=list(mane_status_priorities), ordered=True
        )
    else:
        df["MANE_status"] = pd.Categorical(
            [np.nan] * len(df), categories=list(mane_status_priorities), ordered=True
        )

    # Listed mutations (Include sheet)
    listed_mutations = pd.read_excel(listed_mutations_fn, sheet_name="Include")
    listed_mutations = listed_mutations[listed_mutations["Include"] == "Yes"][
        [
            "Gene",
            "Accession",
            "Mutation_types",
            "Restrict_to_exons",
            "Restrict_to_AA_range",
            "Restrict_to_DNA_range",
            "Restrict_missense_to_AA_change",
        ]
    ].copy()
    listed_mutations = _explode_sep(listed_mutations, "Mutation_types", ";")

    valid_types = {
        "frameshift",
        "nonsense",
        "splicing",
        "missense",
        "nonframeshift deletion",
        "nonframeshift insertion",
    }
    unknown_types = set(listed_mutations["Mutation_types"].dropna().unique()) - valid_types
    if unknown_types:
        raise ValueError(f"Unexpected Mutation_types values in Include sheet: {unknown_types}")

    listed_mutations["mutation_func"] = listed_mutations["Mutation_types"].map(MUT_DICT)
    listed_mutations = _explode_sep(listed_mutations, "mutation_func", ";")

    # is_broad must be computed per exploded mutation_func, not per raw Include row: a single row
    # can combine e.g. "frameshift;nonsense;splicing;missense" with Restrict_missense_to_AA_change
    # set, but that restriction only ever applies to the missense (nonsynonymous SNV) branch (see
    # the "mutation_func != 'nonsynonymous SNV'" bypass below) — the frameshift/nonsense/splicing
    # branches of that same row are still fully unrestricted (broad) despite the column being
    # non-null. Affects genes like PHF6, TP53, GATA2, RAD21, RUNX1, SMC3, EZH2 that combine LOF +
    # position-restricted missense on one row.
    missense_restricted = listed_mutations["Restrict_missense_to_AA_change"].notna() & (
        listed_mutations["mutation_func"] == "nonsynonymous SNV"
    )
    other_restrict_cols = ["Restrict_to_exons", "Restrict_to_AA_range", "Restrict_to_DNA_range"]
    listed_mutations["is_broad"] = listed_mutations[other_restrict_cols].isna().all(axis=1) & ~missense_restricted

    for col in ("Restrict_to_exons", "Restrict_to_AA_range", "Restrict_to_DNA_range", "Restrict_missense_to_AA_change"):
        listed_mutations = _explode_sep(listed_mutations, col, ";")

    listed_mutations["AA_range_start"] = pd.to_numeric(
        listed_mutations["Restrict_to_AA_range"].astype(str).str.replace(r"-.*", "", regex=True),
        errors="coerce",
    )
    listed_mutations["AA_range_end"] = pd.to_numeric(
        listed_mutations["Restrict_to_AA_range"].astype(str).str.replace(r".*-", "", regex=True),
        errors="coerce",
    )
    listed_mutations["DNA_range_start"] = pd.to_numeric(
        listed_mutations["Restrict_to_DNA_range"].astype(str).str.replace(r"-.*", "", regex=True),
        errors="coerce",
    )
    listed_mutations["DNA_range_end"] = pd.to_numeric(
        listed_mutations["Restrict_to_DNA_range"].astype(str).str.replace(r".*-", "", regex=True),
        errors="coerce",
    )
    listed_mutations["AA_change_start"] = pd.to_numeric(
        listed_mutations["Restrict_missense_to_AA_change"].astype(str).str.replace(
            r"[A-Z]", "", regex=True
        ),
        errors="coerce",
    )

    listed_df = df.merge(
        listed_mutations,
        left_on=["transcript_for_join", "mutation_func"],
        right_on=["Accession", "mutation_func"],
        how="inner",
    )
    listed_df = listed_df[
        listed_df["Restrict_to_exons"].isna()
        | (listed_df["sub__exon"] == ("exon" + listed_df["Restrict_to_exons"].astype(str)))
    ]
    listed_df = listed_df[
        listed_df["Restrict_to_AA_range"].isna()
        | (
            (listed_df["sub__aa_position"] >= listed_df["AA_range_start"])
            & (listed_df["sub__aa_position"] <= listed_df["AA_range_end"])
        )
    ]
    listed_df = listed_df[
        listed_df["Restrict_to_DNA_range"].isna()
        | (
            (listed_df["sub__cDNA_start"] >= listed_df["DNA_range_start"])
            & (listed_df["sub__cDNA_end"] <= listed_df["DNA_range_end"])
        )
    ]
    listed_df = listed_df[
        listed_df["Restrict_missense_to_AA_change"].isna()
        | (listed_df["mutation_func"] != "nonsynonymous SNV")
        | (
            (listed_df["mutation_func"] == "nonsynonymous SNV")
            & (
                (listed_df["sub__protein_change"] == listed_df["Restrict_missense_to_AA_change"])
                | (include_non_syn_by_loci & (listed_df["sub__aa_position"] == listed_df["AA_change_start"]))
            )
        )
    ]

    df_cols = list(df.columns)
    listed_row_nums = set(listed_df["row_num"])
    rest_df = df[~df["row_num"].isin(listed_row_nums)].copy()
    listed_df = listed_df[df_cols + ["is_broad"]].copy()
    listed_df["is_listed"] = "Listed"
    rest_df["is_listed"] = "Not listed"
    rest_df["is_broad"] = pd.NA

    df = pd.concat([listed_df, rest_df], ignore_index=True).sort_values("row_num").reset_index(drop=True)

    is_splicing = df["mutation_func"] == "splicing"
    base_aachange = np.where(is_splicing, ".", df["transcript_info_col"].astype(str))
    df["AAChange.refGene.for_join"] = pd.Series(base_aachange, index=df.index).str.replace(
        r":p\..+", "", regex=True
    )
    df["GeneDetail.refGene.for_join"] = np.where(
        is_splicing,
        df["sub__transcript"].astype(str) + ":" + df["sub__exon"].astype(str) + ":c." + df["sub__cDNA_change"].astype(str),
        ".",
    )

    excluded_mutations = pd.read_excel(listed_mutations_fn, sheet_name="Exclude")[
        ["Gene.refGene", "AAChange.refGene", "GeneDetail.refGene"]
    ].copy()
    excluded_mutations["AAChange.refGene"] = excluded_mutations["AAChange.refGene"].astype(str).str.replace(
        r":p\..+", "", regex=True
    )
    excluded_mutations.columns = ["Gene.refGene", "AAChange.refGene.for_join", "GeneDetail.refGene.for_join"]
    excluded_mutations = _explode_sep(excluded_mutations, "AAChange.refGene.for_join", mut_info_sep)
    excluded_mutations = _explode_sep(excluded_mutations, "GeneDetail.refGene.for_join", non_exonic_mut_info_sep)
    excluded_mutations["is_excluded"] = True
    excluded_mutations = excluded_mutations.drop_duplicates()

    df = df.merge(
        excluded_mutations,
        on=["Gene.refGene", "AAChange.refGene.for_join", "GeneDetail.refGene.for_join"],
        how="left",
    )
    df["is_excluded"] = df["is_excluded"].map({True: "Excluded"}).fillna("Not excluded")

    df["merged__sub"] = df.apply(
        lambda r: ":".join(
            _fmt(r[c])
            for c in (
                "sub__gene",
                "sub__transcript",
                "sub__exon",
                "sub__cDNA_change",
                "sub__protein_change",
                "sub__aa_position",
                "MANE_status",
                "is_listed",
                "is_excluded",
            )
        ),
        axis=1,
    )

    mane_priorities = list(mane_status_priorities)

    def _aggregate_group(g):
        has_protein_coding = bool(g["Has.protein.coding"].any())
        mane_codes = g["MANE_status"].cat.codes
        valid_codes = mane_codes[mane_codes >= 0]
        top_mane_status = mane_priorities[valid_codes.min()] if len(valid_codes) > 0 else None
        has_listed = bool((g["is_listed"] == "Listed").any())
        has_excluded = bool((g["is_excluded"] == "Excluded").any())

        if has_listed:
            # is_broad is only ever NA for "Not listed" rows, never for "Listed" ones, so no NA
            # handling is needed here.
            listed_is_broad = g.loc[g["is_listed"] == "Listed", "is_broad"].astype(bool)
            listed_match_type = "broad" if bool(listed_is_broad.all()) else "narrow"
        else:
            listed_match_type = np.nan

        if has_listed:
            sel_idx = g.index[g["is_listed"] == "Listed"][0]
        elif top_mane_status is None:
            sel_idx = None
        elif top_mane_status == "None":
            sel_idx = g.index[-1]
        else:
            matches = g.index[g["MANE_status"] == top_mane_status]
            sel_idx = matches[0] if len(matches) else g.index[-1]

        if sel_idx is None:
            selected_transcript = np.nan
            selected_aa_change = np.nan
        else:
            selected_transcript = g.loc[sel_idx, "sub__transcript"]
            selected_aa_change = g.loc[sel_idx, "sub__protein_change"]

        return pd.Series(
            {
                "Has.protein.coding": has_protein_coding,
                "Top.MANE.status": top_mane_status if top_mane_status is not None else np.nan,
                "Has.listed": has_listed,
                "Has.listed.match_type": listed_match_type,
                "Has.excluded": has_excluded,
                "Selected.transcript": selected_transcript,
                "Selected.AA.change": selected_aa_change,
                merged_mut_info_col: mut_info_sep.join(g["merged__sub"].tolist()),
                mut_info_col: mut_info_sep.join(g[mut_info_col].astype(str).tolist()),
            }
        )

    agg_df = df.groupby("__row_id__", sort=False).apply(_aggregate_group).reset_index()

    orig_static = input_df[["__row_id__"] + orig_columns]
    output_df = orig_static.merge(agg_df, on="__row_id__", how="left").drop(columns=["__row_id__"])
    return output_df


def _ensure_numeric_vaf(df, col, min_vaf, log):
    if not pd.api.types.is_numeric_dtype(df[col]):
        log(f"WARNING: {col} column is not numeric, checking if it contains percentage strings...")
        str_vals = df[col].astype(str)
        if str_vals.str.endswith("%").all():
            df[col] = str_vals.str.rstrip("%").astype(float) / 100
            log(f"Converted {col} percentage strings to numeric values.")
        else:
            raise ValueError(f"{col} column is not numeric and does not contain percentage strings")
    df = df[df[col] >= min_vaf]
    log(f"After {col} >= {min_vaf} filtering: {len(df)} rows")
    return df


def filter_mutations(
    input_file,
    output_dir,
    run_name,
    gene_column="Gene.refGene",
    chr_column="CHR",
    pos_column="POS",
    sample_column="SAMPLE_ID",
    listed_mutations_file_name=DEFAULT_LISTED_MUTATIONS_FN,
    mane_ifn=DEFAULT_MANE_IFN,
    overwrite=True,
    min_VAF=0.005,
    min_depth=20,
    max_gnomad41_AF=0.001,
    remove_intronic=True,
    remove_synonymous=True,
    include_non_syn_by_loci=True,
    dups_input_file=False,
    sample_specific_columns="VAF;Depth;Sample_Name",
    min_supporting_samples=1,
    min_avg_VAF=-1.0,
):
    """Port of Mutation_Curation/R/FilterMutations.R's main filtering flow.

    Returns a dict with 'listed_df', 'rest_df', 'listed_fn', 'rest_fn', 'log_fn'.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_fn = os.path.join(output_dir, f"FilterMutations_{run_name}_listed.tsv")
    rest_fn = os.path.join(output_dir, f"FilterMutations_{run_name}_filtered_yet_not_listed.tsv")
    log_fn = os.path.join(output_dir, f"FilterMutations_{run_name}.log")

    if not overwrite:
        for out_fn in (output_fn, rest_fn, log_fn):
            if os.path.exists(out_fn):
                raise FileExistsError(f"Output file {out_fn} exists and overwrite is set to False!")

    log_fh = open(log_fn, "w")

    def log(message):
        print(message)
        log_fh.write(message + "\n")

    def check_df(df):
        if len(df) == 0:
            log("No mutations left after filtering!")
            log_fh.close()
            raise ValueError("No mutations left after filtering!")

    try:
        log(f"Starting FilterMutations for run {run_name}")
        log("=" * 80)

        df = pd.read_csv(input_file, sep="\t")
        required_columns = [gene_column, sample_column, chr_column, pos_column, "REF", "ALT", "VAF", "gnomad41_genome_AF", "Func.refGene", "AAChange.refGene"]
        if dups_input_file:
            required_columns += ["VAF_DUP", "DEPTH_DUP"]
        else:
            required_columns += [c for c in sample_specific_columns.split(";") if c]
        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise ValueError(f"Column(s) {missing} not found in input file!")
        log(f"Input file {input_file} read: {len(df)} rows, {len(df.columns)} columns")
        check_df(df)

        df_f = df.copy()

        vaf_cols = ["VAF", "VAF_DUP"] if dups_input_file else ["VAF"]
        for c_vaf_col in vaf_cols:
            df_f = _ensure_numeric_vaf(df_f, c_vaf_col, min_VAF, log)
            check_df(df_f)

        depth_cols = ["DEPTH", "DEPTH_DUP"] if dups_input_file else ["Depth"]
        for c_depth_col in depth_cols:
            df_f = df_f[df_f[c_depth_col] >= min_depth]
            log(f"After {c_depth_col} >= {min_depth} filtering: {len(df_f)} rows")
            check_df(df_f)

        if remove_intronic:
            df_f = df_f[df_f["Func.refGene"].isin(["exonic", "splicing"])]
            log(f"After removing intronic mutations, left with exonic/splicing: {len(df_f)} rows")
            check_df(df_f)

        if remove_synonymous:
            df_f = df_f[df_f["ExonicFunc.refGene"] != "synonymous SNV"]
            log(f"After removing synonymous mutations: {len(df_f)} rows")
            check_df(df_f)

        df_f = annotate_mutations_df(
            input_df=df_f,
            listed_mutations_fn=listed_mutations_file_name,
            mane_ifn=mane_ifn,
            include_non_syn_by_loci=include_non_syn_by_loci,
        )
        log("Annotated mutations with listed/excluded mutations and MANE information")

        n_excluded = int(df_f["Has.excluded"].sum())
        df_f = df_f[~df_f["Has.excluded"]]
        log(f"After removing {n_excluded} mutations in the exclusion list left with: {len(df_f)} rows")
        check_df(df_f)

        if dups_input_file:
            if min_supporting_samples > 2:
                raise ValueError("min_supporting_samples must be <= 2 when dups_input_file is True")
            df_f["avg_VAF"] = df_f[["VAF", "VAF_DUP"]].mean(axis=1, skipna=True)
            df_f["mean_depth"] = df_f[["DEPTH", "DEPTH_DUP"]].mean(axis=1)
            log("Calculated average VAF and DEPTH across duplicate samples")
        else:
            df_f["avg_VAF"] = df_f["VAF"]
            df_f["mean_depth"] = df_f["Depth"]
            log("Calculated average VAF and Depth across duplicate samples")

        if min_avg_VAF == -1:
            f_min_avg_VAF = df_f["avg_VAF"] >= (min_depth / 10) / (min_depth + df_f["mean_depth"])
        else:
            f_min_avg_VAF = df_f["avg_VAF"] >= min_avg_VAF

        numeric_gnomad41_AF = pd.to_numeric(df_f["gnomad41_genome_AF"], errors="coerce")

        suspected_common = df_f[df_f["Has.listed"] & f_min_avg_VAF & numeric_gnomad41_AF.notna() & (numeric_gnomad41_AF > max_gnomad41_AF)]
        if len(suspected_common) > 0:
            sus_info = suspected_common[["Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "Selected.AA.change", "avg_VAF"]].astype(str).agg(":".join, axis=1)
            log(
                f"Found {len(suspected_common)} listed mutations with population AF above "
                f"{max_gnomad41_AF} cutoff, will be filtered out: {','.join(sus_info)}"
            )

        df_f = df_f[numeric_gnomad41_AF.isna() | (numeric_gnomad41_AF <= max_gnomad41_AF)]
        f_min_avg_VAF = f_min_avg_VAF.loc[df_f.index]
        log(f"After gnomad41_AF <= {max_gnomad41_AF} filtering: {len(df_f)} rows")
        check_df(df_f)

        df_ff = df_f[df_f["Has.listed"] & f_min_avg_VAF].copy()
        log(f"Listed mutations found: {len(df_ff)} rows")

        df_f_rest = df_f[~(df_f["Has.listed"] & f_min_avg_VAF)].copy()
        log(f"Non-listed mutations yet survived filters so far found: {len(df_f_rest)} rows")

        df_ff.to_csv(output_fn, sep="\t", na_rep="NA", index=False)
        df_f_rest.to_csv(rest_fn, sep="\t", na_rep="NA", index=False)
    finally:
        log_fh.close()

    return {"listed_df": df_ff, "rest_df": df_f_rest, "listed_fn": output_fn, "rest_fn": rest_fn, "log_fn": log_fn}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter mutations from multiple samples (annovar output) by VAF, mutation type, "
        "coordinate within specific genes and a curated list of mutations (hotspots).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-file", required=True, help="Input mutation file (tab-delimited, annovar output columns).")
    parser.add_argument("--gene-column", default="Gene.refGene")
    parser.add_argument("--chr-column", default="CHR")
    parser.add_argument("--pos-column", default="POS")
    parser.add_argument("--sample-column", default="SAMPLE_ID")
    parser.add_argument("--listed-mutations-file-name", default=DEFAULT_LISTED_MUTATIONS_FN)
    parser.add_argument("--mane-ifn", default=DEFAULT_MANE_IFN)
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--run-name", required=True, help="Name of run, will be part of output files.")
    parser.add_argument("--overwrite", type=lambda v: v.lower() != "false", default=True)
    parser.add_argument("--min-VAF", type=float, default=0.005)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--max-gnomad41-AF", type=float, default=0.001)
    parser.add_argument("--remove-intronic", type=lambda v: v.lower() != "false", default=True)
    parser.add_argument("--remove-synonymous", type=lambda v: v.lower() != "false", default=True)
    parser.add_argument("--include-non-syn-by-loci", type=lambda v: v.lower() != "false", default=True)
    parser.add_argument("--dups-input-file", type=lambda v: v.lower() != "false", default=False)
    parser.add_argument("--sample-specific-columns", default="VAF;Depth;Sample_Name")
    parser.add_argument("--min-supporting-samples", type=int, default=1)
    parser.add_argument("--min-avg-VAF", type=float, default=-1.0)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    filter_mutations(
        input_file=args.input_file,
        output_dir=args.output_dir,
        run_name=args.run_name,
        gene_column=args.gene_column,
        chr_column=args.chr_column,
        pos_column=args.pos_column,
        sample_column=args.sample_column,
        listed_mutations_file_name=args.listed_mutations_file_name,
        mane_ifn=args.mane_ifn,
        overwrite=args.overwrite,
        min_VAF=args.min_VAF,
        min_depth=args.min_depth,
        max_gnomad41_AF=args.max_gnomad41_AF,
        remove_intronic=args.remove_intronic,
        remove_synonymous=args.remove_synonymous,
        include_non_syn_by_loci=args.include_non_syn_by_loci,
        dups_input_file=args.dups_input_file,
        sample_specific_columns=args.sample_specific_columns,
        min_supporting_samples=args.min_supporting_samples,
        min_avg_VAF=args.min_avg_VAF,
    )
