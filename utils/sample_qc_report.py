#!/usr/bin/env python3
"""
Sample pass/fail QC report for a sequencing run.

A library "works" if >= MIN_FRACTION (default 80%) of on-target panel bases
have deduplicated coverage >= MIN_BASE_COV (default 20x) -- i.e. the
'base >=20 dedup %on-target cov' column produced by utils/mip_run_statistics.py.

Input: either {panel}_{seq_run}_statistics.csv or {seq_run}_dedup_coverage.csv
(both contain the needed column).

Output:
  - console summary (counts + list of failed libraries)
  - qc_pass_fail_{run}.xlsx with 3 sheets: Libraries, Donors, Summary
  - qc_pass_fail_{run}_libraries.csv (same as the Libraries sheet, for grep/awk)

Usage on the cluster:
  python sample_qc_report.py AML24_NovaseqR132_statistics.csv --run-name NovaseqR132
  python sample_qc_report.py stats.csv --min-base-cov 20 --min-fraction 0.8
"""

import argparse
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Pass/fail QC report (coverage-based)")
    ap.add_argument("stats_csv", help="statistics CSV from mip_run_statistics.py")
    ap.add_argument("--min-base-cov", type=int, default=20,
                    help="per-base dedup coverage threshold (default 20)")
    ap.add_argument("--min-fraction", type=float, default=0.8,
                    help="minimum fraction of panel bases at that coverage (default 0.8)")
    ap.add_argument("--run-name", default=None, help="run name for output file names")
    ap.add_argument("--out-dir", default=".", help="where to write the report files")
    args = ap.parse_args()

    cov_col = f"base >={args.min_base_cov} dedup %on-target cov"
    run = args.run_name or os.path.basename(args.stats_csv).split("_statistics")[0]

    df = pd.read_csv(args.stats_csv)
    if "Sample" not in df.columns:
        sys.exit(f"No 'Sample' column in {args.stats_csv}. Columns: {list(df.columns)}")
    if cov_col not in df.columns:
        sys.exit(
            f"Column '{cov_col}' not found in {args.stats_csv}.\n"
            f"Available: {list(df.columns)}\n"
            f"Re-run mip_run_statistics.py with --min-base-cov {args.min_base_cov}, "
            f"or pass the matching --min-base-cov here."
        )

    lib = df[["Sample", cov_col]].copy()
    lib[cov_col] = pd.to_numeric(lib[cov_col], errors="coerce")
    lib["Status"] = (lib[cov_col] >= args.min_fraction).map({True: "Pass", False: "Fail"})
    lib.loc[lib[cov_col].isna(), "Status"] = "Fail (no coverage data)"

    # Sample names look like EXP_PANEL_DONOR_INDEX, e.g. AEX35_ARCH_PBH13_S148.
    toks = lib["Sample"].astype(str).str.split("_")
    lib["donor"] = toks.str[2]
    lib["is_control"] = lib["donor"].fillna("").str.startswith(("PC", "NC"))
    lib = lib.sort_values(cov_col, ascending=False)

    # Donor level: a donor works if at least one of its technical libraries passes.
    donors = (
        lib[~lib["is_control"]]
        .groupby("donor", dropna=False)
        .agg(
            n_libraries=("Sample", "size"),
            n_pass=("Status", lambda s: int((s == "Pass").sum())),
            best_coverage_fraction=(cov_col, "max"),
            libraries=("Sample", lambda s: ";".join(s)),
        )
        .reset_index()
    )
    donors["Donor_Status"] = (donors["n_pass"] > 0).map({True: "Pass", False: "Fail"})
    donors = donors.sort_values(["Donor_Status", "best_coverage_fraction"],
                                ascending=[True, False])

    real = lib[~lib["is_control"]]
    summary = pd.DataFrame(
        [
            ("Criterion", f">= {args.min_fraction:.0%} of panel bases at >= {args.min_base_cov}x (dedup)"),
            ("Libraries total (excl. controls)", len(real)),
            ("Libraries pass", int((real["Status"] == "Pass").sum())),
            ("Libraries fail", int((real["Status"] != "Pass").sum())),
            ("Library pass rate", f"{(real['Status'] == 'Pass').mean():.1%}" if len(real) else "n/a"),
            ("Donors total", len(donors)),
            ("Donors pass (>=1 library)", int((donors["Donor_Status"] == "Pass").sum())),
            ("Donors fail", int((donors["Donor_Status"] == "Fail").sum())),
            ("Control libraries (PC/NC)", int(lib["is_control"].sum())),
        ],
        columns=["Metric", "Value"],
    )

    os.makedirs(args.out_dir, exist_ok=True)
    xlsx_path = os.path.join(args.out_dir, f"qc_pass_fail_{run}.xlsx")
    csv_path = os.path.join(args.out_dir, f"qc_pass_fail_{run}_libraries.csv")
    with pd.ExcelWriter(xlsx_path) as xl:
        lib.to_excel(xl, sheet_name="Libraries", index=False)
        donors.to_excel(xl, sheet_name="Donors", index=False)
        summary.to_excel(xl, sheet_name="Summary", index=False)
    lib.to_csv(csv_path, index=False)

    print(summary.to_string(index=False))
    failed = real[real["Status"] != "Pass"]
    if len(failed):
        print(f"\nFailed libraries ({len(failed)}):")
        print(failed[["Sample", cov_col, "Status"]].to_string(index=False))
    else:
        print("\nAll libraries passed.")
    print(f"\nWrote: {xlsx_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
