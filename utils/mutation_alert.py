"""
Mutation alert — scan a curated mutation table for a configured list of target
mutations, plus optional CALR type-1/type-2 result tables for high-VAF calls,
and email a configured list of recipients if any are found.

Usage:
  python mutation_alert.py --mutations-file /path/to/FilterMutations_..._VAF_0.05.tsv \
      --config mutation_alert_config.yaml --run-name PERIBLOOD_merged_arch \
      --calr-type1-file /path/to/CALR_type1_results_{seq_run}.tsv \
      --calr-type2-file /path/to/CALR_type2_results_{seq_run}.tsv

Config file (see mutation_alert_config.yaml) defines:
  - recipients: list of email addresses to notify
  - target_mutations: list of dicts, each a subset of
    {CHR, POS, REF, ALT, Gene.refGene, Selected.AA.change} — a mutation row
    matches if ALL fields specified in a dict are equal to that row's values
  - min_CALR_del52_VAF / min_CALR_ins5_VAF: VAF thresholds for the CALR
    type-1/type-2 result files; leave unset to skip that check
  - healthy_donor_prefixes / min_healthy_donor_VAF: alert on ANY mutation
    (regardless of gene) with avg_VAF >= min_healthy_donor_VAF found in a
    donor whose identifier (3rd underscore token of Sample_Name) starts with
    one of healthy_donor_prefixes; leave the prefix list empty or the
    threshold unset to skip this check
  - smtp: host/port/use_tls/username/password/from_address
"""

import argparse
import io
import os
import smtplib
import sys
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import pandas as pd
import yaml

REPORT_COLUMNS = ["Sample_Name", "Gene.refGene", "Selected.AA.change", "CHR", "POS", "REF", "ALT", "avg_VAF", "mean_depth"]

def load_config(config_path: str) -> dict:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f) or {}
    config.setdefault("recipients", [])
    config.setdefault("target_mutations", [])
    config.setdefault("min_CALR_del52_VAF", None)
    config.setdefault("min_CALR_ins5_VAF", None)
    config.setdefault("healthy_donor_prefixes", [])
    config.setdefault("min_healthy_donor_VAF", None)
    config.setdefault("smtp", {})
    smtp = config["smtp"]
    smtp.setdefault("host", "localhost")
    smtp.setdefault("port", 25)
    smtp.setdefault("use_tls", False)
    smtp.setdefault("username", None)
    smtp.setdefault("password", None)
    smtp.setdefault("from_address", "arch-alerts@localhost")
    return config


def find_matches(df: pd.DataFrame, target_mutations: list) -> pd.DataFrame:
    if not target_mutations:
        return df.iloc[0:0]

    combined_mask = pd.Series(False, index=df.index)
    for spec in target_mutations:
        mask = pd.Series(True, index=df.index)
        for col, val in spec.items():
            if col not in df.columns:
                mask &= False
                continue
            mask &= df[col].astype(str) == str(val)
        combined_mask |= mask

    return df[combined_mask].drop_duplicates()


def find_healthy_donor_alerts(df: pd.DataFrame, healthy_donor_prefixes: list, min_vaf, vaf_column: str = "avg_VAF") -> pd.DataFrame:
    if not healthy_donor_prefixes or min_vaf is None:
        return df.iloc[0:0]
    if "Sample_Name" not in df.columns or vaf_column not in df.columns:
        return df.iloc[0:0]

    donor = df["Sample_Name"].str.split("_").str[2]
    is_healthy = donor.str.startswith(tuple(healthy_donor_prefixes)).fillna(False)
    high_vaf = df[vaf_column].astype(float) >= float(min_vaf)
    return df[is_healthy & high_vaf]


def find_calr_alerts(path: str, vaf_column: str, min_vaf) -> pd.DataFrame:
    if not path or min_vaf is None:
        return pd.DataFrame()
    if not os.path.exists(path):
        print(f"CALR results file not found, skipping: {path}", file=sys.stderr)
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", dtype=str)
    return df[df[vaf_column].astype(float) > float(min_vaf)]


def send_alert_email(sections: list, recipients: list, run_name: str, smtp_cfg: dict):
    if not recipients:
        print("No recipients configured — skipping email send.", file=sys.stderr)
        return

    total = sum(len(section["df"]) for section in sections)
    body_lines = [f"Alert(s) found in run: {run_name or '(unknown)'}", ""]
    for section in sections:
        body_lines += [
            f"{section['label']} ({len(section['df'])} row(s)):",
            section["df"].to_string(index=False),
            "",
        ]
    body_lines.append("Full details for each section are attached as separate CSV files.")

    msg = MIMEMultipart()
    msg["Subject"] = f"[ARCH mutation alert] {total} alert(s) found in {run_name or 'run'}"
    msg["From"] = smtp_cfg["from_address"]
    msg["To"] = ", ".join(recipients)
    msg.attach(MIMEText("\n".join(body_lines)))

    for section in sections:
        csv_buffer = io.StringIO()
        section["df"].to_csv(csv_buffer, index=False)
        attachment = MIMEApplication(csv_buffer.getvalue().encode("utf-8"), _subtype="csv")
        attachment.add_header("Content-Disposition", "attachment", filename=section["filename"])
        msg.attach(attachment)

    try:
        with smtplib.SMTP(smtp_cfg["host"], smtp_cfg["port"]) as server:
            if smtp_cfg.get("use_tls"):
                server.starttls()
            if smtp_cfg.get("username"):
                server.login(smtp_cfg["username"], smtp_cfg["password"])
            server.sendmail(smtp_cfg["from_address"], recipients, msg.as_string())
        print(f"Alert email sent to: {', '.join(recipients)}")
    except Exception as e:
        print(f"ERROR: failed to send alert email: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Scan a mutation table (and optional CALR result tables) for alerts and email recipients if found.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mutations-file", required=True, help="Path to the mutation TSV to scan")
    parser.add_argument(
        "--config",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "mutation_alert_config.yaml"),
        help="Path to the mutation_alert YAML config file",
    )
    parser.add_argument("--run-name", default=None, help="Run name, included in the email subject/body")
    parser.add_argument("--calr-type1-file", default=None, help="Path to CALR_type1_results_{seq_run}.tsv")
    parser.add_argument("--calr-type2-file", default=None, help="Path to CALR_type2_results_{seq_run}.tsv")

    args = parser.parse_args()

    config = load_config(args.config)

    df = pd.read_csv(args.mutations_file, sep="\t", dtype=str)
    target_matches = find_matches(df, config["target_mutations"])
    calr1_matches = find_calr_alerts(args.calr_type1_file, "VAF_DEL52", config["min_CALR_del52_VAF"])
    calr2_matches = find_calr_alerts(args.calr_type2_file, "VAF_INS5", config["min_CALR_ins5_VAF"])
    healthy_donor_matches = find_healthy_donor_alerts(df, config["healthy_donor_prefixes"], config["min_healthy_donor_VAF"])

    if target_matches.empty and calr1_matches.empty and calr2_matches.empty and healthy_donor_matches.empty:
        print("No alerts found — no email sent.")
        return

    sections = []
    if not target_matches.empty:
        report_cols = [c for c in REPORT_COLUMNS if c in target_matches.columns]
        sort_spec = [(col, asc) for col, asc in [("Gene.refGene", True), ("avg_VAF", False)] if col in report_cols]
        target_df = target_matches[report_cols]
        if sort_spec:
            target_df = target_df.sort_values(by=[c for c, _ in sort_spec], ascending=[a for _, a in sort_spec])
        sections.append({
            "label": "Target mutations",
            "df": target_df,
            "filename": f"mutation_alert_{args.run_name or 'run'}.csv",
        })
    if not calr1_matches.empty:
        sections.append({
            "label": "CALR type-1 (del52) high-VAF samples",
            "df": calr1_matches.sort_values(by="VAF_DEL52", ascending=False),
            "filename": f"CALR_type1_alert_{args.run_name or 'run'}.csv",
        })
    if not calr2_matches.empty:
        sections.append({
            "label": "CALR type-2 (ins5) high-VAF samples",
            "df": calr2_matches.sort_values(by="VAF_INS5", ascending=False),
            "filename": f"CALR_type2_alert_{args.run_name or 'run'}.csv",
        })
    if not healthy_donor_matches.empty:
        report_cols = [c for c in REPORT_COLUMNS if c in healthy_donor_matches.columns]
        healthy_df = healthy_donor_matches[report_cols]
        if "avg_VAF" in report_cols:
            healthy_df = healthy_df.sort_values(by="avg_VAF", ascending=False)
        sections.append({
            "label": "Healthy-donor high-VAF mutations",
            "df": healthy_df,
            "filename": f"healthy_donor_alert_{args.run_name or 'run'}.csv",
        })

    for section in sections:
        print(f"{section['label']}: {len(section['df'])} row(s)")
        print(section["df"].to_string(index=False))

    send_alert_email(sections, config["recipients"], args.run_name, config["smtp"])


if __name__ == "__main__":
    main()
