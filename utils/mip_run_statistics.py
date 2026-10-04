import os
import re
import subprocess
import time
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

# ==========================================
# Utility Functions
# ==========================================
def get_mip_run_statistics_df(run_name, runs_base_dir, run_type='arch', stat_base_fn='total_Statistics.csv', rebuild=False):
    analysis_dir = os.path.join(runs_base_dir, run_name, f"analysis_{run_type}")
    stat_ifn = os.path.join(analysis_dir, stat_base_fn)

    if not os.path.exists(stat_ifn) or rebuild:
        if os.path.exists(stat_ifn):
            os.remove(stat_ifn)

        print('Collecting statistics from jobs...')
        # List files/directories in the output path that start with "job"
        job_dirs = [d for d in os.listdir(analysis_dir) if d.startswith("job")]
        job_ids = sorted([int(d.replace("job", "")) for d in job_dirs])

        if run_type == 'geno':
            run_type = 'genotype'

        for job_id in job_ids:
            job_stat_ifn = os.path.join(analysis_dir, f"job{job_id}", "paired", "4_report", f"{run_type}_{run_name}_{job_id}_Statistics.xlsx")
            
            if not os.path.exists(job_stat_ifn):
                continue
                
            # Read Excel
            df = pd.read_excel(job_stat_ifn, sheet_name="Summary")
            
            # Write to CSV (append mode, writes header only if file doesn't exist yet)
            df.to_csv(stat_ifn, index=False, mode='a', header=not os.path.exists(stat_ifn))

    # Read the combined dataframe back
    print(f'Reading combined statistics from {stat_ifn}')
    stat_df = pd.read_csv(stat_ifn)
    return stat_df

def plot_scatter_with_repel(data, x_col, y_col, title, hue_by='exp', style_by='control', pt_size=100, pdf=None, h_lines=None, v_lines=None, label_mask=None):
    fig = plt.figure(figsize=(10, 6))

    # Seaborn natively handles coloring by 'exp' and shaping by 'control'
    sns.scatterplot(data=data, x=x_col, y=y_col, hue=hue_by, style=style_by, s=pt_size)

    # Add text labels only for samples flagged via label_mask (None = label everyone)
    texts = []
    for _, row in data.iterrows():
        if pd.notna(row[x_col]) and pd.notna(row[y_col]):
            if label_mask is None or label_mask.loc[row.name]:
                texts.append(plt.text(row[x_col], row[y_col], str(row['donor']), fontsize=9))

    # Repel texts so they don't overlap
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5),
                    prevent_crossings=False, time_lim=0.5)
    
    if h_lines:
        for y in h_lines:
            plt.axhline(y, color='black', linestyle='--')
    if v_lines:
        for x in v_lines:
            plt.axvline(x, color='black', linestyle='--')

    plt.title(title)
    # Move legend outside the plot area to match ggplot defaults
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    if pdf is not None:
        pdf.savefig(fig)
        plt.close(fig)
    else:
        plt.show()

# ==========================================
# Library Complexity & Dedup Coverage
# ==========================================
def _wait_for_lsf_job(job_id, poll_interval=30):
    print(f"Waiting for LSF job array {job_id} to finish (polling every {poll_interval}s)...")
    while True:
        result = subprocess.run(
            ["bjobs", "-noheader", str(job_id)],
            capture_output=True, text=True
        )
        # All tasks purged from active list — check for any that exited with failure
        if result.returncode != 0 or "not found" in result.stderr.lower():
            exit_check = subprocess.run(
                ["bjobs", "-noheader", "-x", str(job_id)],
                capture_output=True, text=True
            )
            if exit_check.stdout.strip():
                raise RuntimeError(f"LSF job {job_id} has failed tasks:\n{exit_check.stdout.strip()}")
            break

        all_lines = result.stdout.strip().splitlines()
        lines = [l for l in all_lines if len(l.split()) > 2]
        skipped = [l for l in all_lines if l.strip() and len(l.split()) <= 2]
        if skipped:
            print(f"  [bjobs] skipped {len(skipped)} short line(s): {skipped}")
        statuses = [line.split()[2] for line in lines]

        n_exit = sum(1 for s in statuses if s == 'EXIT')
        if n_exit:
            raise RuntimeError(f"LSF job {job_id} has {n_exit} failed task(s).")

        n_active = sum(1 for s in statuses if s in {'PEND', 'RUN', 'SSUSP', 'USUSP', 'PSUSP'})
        n_done = sum(1 for s in statuses if s == 'DONE')
        if n_active == 0:
            break

        print(f"  Job {job_id}: {n_active} active, {n_done} done...")
        time.sleep(poll_interval)
    print(f"LSF job {job_id} completed.")


def estimate_library_complexity(run_dir, rebuild):
    comp_estimate_script = os.path.join(run_dir, "submit_lib_size_estimate_array.bsub")
    if not os.path.exists(comp_estimate_script):
        raise FileNotFoundError(
            f"File not found: {comp_estimate_script}, it should have been generated by utils/gen_mip_processing_cmds.py"
        )
    lib_comp_odir = os.path.join(run_dir, 'results')
    if rebuild or not os.path.exists(lib_comp_odir):
        with open(comp_estimate_script, 'r') as f:
            result = subprocess.run(["bsub"], stdin=f, check=True, capture_output=True, text=True, cwd=run_dir)
        bsub_output = result.stdout + result.stderr
        print(bsub_output.strip())
        match = re.search(r'Job <(\d+)>', bsub_output)
        if not match:
            raise RuntimeError(f"Could not parse job ID from bsub output: {bsub_output!r}")
        _wait_for_lsf_job(match.group(1))
    else:
        print(f"Library complexity results already exist at {lib_comp_odir}, skipping.")

    # Collect Picard EstimateLibraryComplexity results (row 7 = header, row 8 = data)
    metric_suffix = ".filtered.sorted.rg.realigned.rn_complexity_metrics.txt"
    metric_files = list(Path(lib_comp_odir).rglob(f"*{metric_suffix}"))
    if not metric_files:
        print(f"Warning: no complexity metrics files found under {lib_comp_odir}")
        return None

    dfs = []
    for f in metric_files:
        sample_name = f.name.replace(metric_suffix, '')
        df = pd.read_csv(f, sep='\t', skiprows=6, nrows=1)
        df.insert(0, 'Sample', sample_name)
        dfs.append(df)

    complexity_df = pd.concat(dfs, ignore_index=True)
    print(f"Collected library complexity metrics for {len(complexity_df)} samples.")
    return complexity_df


def check_dedup_coverage(run_dir, seq_run_name, min_base_cov, output_dir):
    file_pattern = f"job*/paired/3_stats/*_{seq_run_name}_*_per_base_dedup.tsv"
    print(file_pattern)
    cov_files = list(Path(run_dir).glob(file_pattern))
    print(f"Found {len(cov_files)} coverage files")
    samples_counts_dict = {}

    for cov_file in cov_files:
        cov_df = pd.read_csv(cov_file, sep='\t')
        cov_df.set_index([cov_df.columns[0], cov_df.columns[1]], inplace=True)
        for col in cov_df.columns:
            samples_counts_dict[col] = (cov_df[col] >= min_base_cov).mean()

    samples_counts_df = pd.DataFrame(
        list(samples_counts_dict.items()),
        columns=['Sample', f'base >={min_base_cov} dedup %on-target cov']
    )
    samples_counts_df = samples_counts_df.sort_values(
        by=f'base >={min_base_cov} dedup %on-target cov', ascending=False
    )
    print(samples_counts_df.to_string())

    dedup_out = os.path.join(output_dir, f"{seq_run_name}_dedup_coverage.csv")
    samples_counts_df.to_csv(dedup_out, index=False)
    print(f"Exported dedup coverage to {dedup_out}")

    return samples_counts_df


# ==========================================
# FastQC Summary
# ==========================================
def get_fastqc_summary_df(fastqc_dir):
    tsv_files = sorted(Path(fastqc_dir).glob("fastqc_summary_*.tsv"))
    if not tsv_files:
        print(f"Warning: no FastQC summary files found under {fastqc_dir}")
        return None

    read_dfs = []
    for tsv_file in tsv_files:
        read = re.match(r"fastqc_summary_(.+)\.tsv$", tsv_file.name).group(1)
        df = pd.read_csv(tsv_file, sep='\t')
        content_col = df.columns[0]
        sample_names = df[content_col].str.split(f'_{read}').str[0]
        df = df.drop(columns=[content_col])
        df = df.rename(columns={c: f'{read}_{c}' for c in df.columns})
        df.insert(0, 'Sample', sample_names)
        read_dfs.append(df)

    fastqc_df = read_dfs[0]
    for df in read_dfs[1:]:
        fastqc_df = pd.merge(fastqc_df, df, on='Sample', how='outer')

    print(f"Collected FastQC summary metrics for {len(fastqc_df)} samples across {len(tsv_files)} read file(s).")
    return fastqc_df

# ==========================================
# Main Execution Flow
# ==========================================
def main(output_dir, seq_run_name, panel, rebuild, runs_base_dir, min_f_dedup_cov, sample_fixed_column, min_base_cov, min_lib_size, fastqc_dir=None):

    run_dir = os.path.join(runs_base_dir, seq_run_name, f"analysis_{panel}")

    if fastqc_dir is None:
        fastqc_dir = os.path.join(runs_base_dir, seq_run_name, 'fastq_PERIBLOOD', 'Reports', 'FastQC')


    # ------------------------------------------
    # 3. Load Data
    # ------------------------------------------
    print(f"Fetching statistics for {seq_run_name} ({panel})...")
    stat_base_fn = f"{panel}_{seq_run_name}_statistics.csv"
    stat_df = get_mip_run_statistics_df(
        run_name=seq_run_name, 
        runs_base_dir=runs_base_dir, 
        run_type=panel, 
        stat_base_fn=stat_base_fn,
        rebuild=rebuild
    )
    print(f'stat_df.columns: {stat_df.columns}')

    # ------------------------------------------
    # 2. Estimate Library Complexity
    # ------------------------------------------
    lib_comp_df = estimate_library_complexity(run_dir, rebuild)
    if lib_comp_df is not None:
        print(f"Library complexity metrics found.")
    else:
        print("No library complexity metrics found")

    # ------------------------------------------
    # 3. Check Dedup Coverage
    # ------------------------------------------
    dedup_f_cov_df = check_dedup_coverage(run_dir, seq_run_name, min_base_cov, output_dir)
    if dedup_f_cov_df is not None:
        print(f"Dedup coverage metrics found.")
    else:
        print("No dedup coverage metrics found")

    # ------------------------------------------
    # 4. Collect FastQC Summary
    # ------------------------------------------
    fastqc_df = get_fastqc_summary_df(fastqc_dir)
    if fastqc_df is not None:
        print(f"FastQC summary metrics found.")
    else:
        print("No FastQC summary metrics found")

    # Data Manipulation
    x = stat_df.copy()

    # Merge library complexity metrics
    if lib_comp_df is None:
        raise RuntimeError(f"No library complexity metrics found under {os.path.join(run_dir, 'results')} — run estimate_library_complexity first.")
    x = pd.merge(x, lib_comp_df, on='Sample', how='left')

    # Merge dedup coverage metrics
    x = pd.merge(x, dedup_f_cov_df, on='Sample', how='left')

    # Merge FastQC summary metrics
    if fastqc_df is not None:
        x = pd.merge(x, fastqc_df, on='Sample', how='left')

    # Split Sample into exp, panel_name, donor, index_num
    x[['exp', 'panel_name', 'donor', 'index_num']] = x[sample_fixed_column].str.split('_', expand=True)

    # Create control column by checking if 'donor' contains 'PC' or 'NC'
    x['control'] = x['donor'].str.contains('[PN]C', regex=True)

    # Pre-calculate log columns for later plotting
    x['log2_Total'] = np.log2(x['Total'])
    x['log2_Mapped'] = np.log2(x['Mapped'])

    # ------------------------------------------
    # Classify Samples (moved up so plots below can label only samples of interest)
    # ------------------------------------------
    base_cov_col = f'base >={min_base_cov} dedup %on-target cov'
    choices = ['Pass', 'Repeat Infiniseq', 'Repeat DNA extraction']

    class_conditions = [
        x[base_cov_col] >= min_f_dedup_cov,
        (x[base_cov_col] <= min_f_dedup_cov) & (x['ESTIMATED_LIBRARY_SIZE'] >= min_lib_size),
        (x[base_cov_col] <= min_f_dedup_cov) & (x['ESTIMATED_LIBRARY_SIZE'] < min_lib_size),
    ]

    x['Status'] = np.select(class_conditions, choices, default='Unknown')

    # Only label samples that need attention (+ controls) on the QC scatter plots below —
    # labeling every sample doesn't scale (adjust_text's label-overlap detection blows up
    # well past quadratic cost as the number of mutually-overlapping labels grows).
    x['needs_label'] = (x['Status'] != 'Pass') | x['control']
    MAX_LABELED_SAMPLES = 150
    n_flagged = x['needs_label'].sum()
    if n_flagged > MAX_LABELED_SAMPLES:
        print(f"Warning: {n_flagged} samples flagged for labeling, exceeds cap of {MAX_LABELED_SAMPLES}; "
              f"labeling only the first {MAX_LABELED_SAMPLES}.")
        flagged_idx = x[x['needs_label']].index[:MAX_LABELED_SAMPLES]
        x['needs_label'] = x.index.isin(flagged_idx)

    # ------------------------------------------
    # PDF Report Setup
    # ------------------------------------------
    pdf_path = os.path.join(output_dir, f"{panel}_{seq_run_name}_report.pdf")
    print(f"Generating PDF report for {len(x)} samples at: {pdf_path}")
    
    with PdfPages(pdf_path) as pdf:
        # ------------------------------------------
        # 2. Mapping Success
        # ------------------------------------------
        plot_scatter_with_repel(x, 'Mapped%', 'On Target%', f"{seq_run_name}: {panel}", pdf=pdf, label_mask=x['needs_label'])
        plot_scatter_with_repel(x, 'Mapped%', '>x100(%).fragments', f"{seq_run_name}: {panel}", pdf=pdf, label_mask=x['needs_label'])
        plot_scatter_with_repel(x, 'Uniformity.per.base(%)', 'On Target%', f"{seq_run_name}: {panel}", pdf=pdf, label_mask=x['needs_label'])
        
        # ------------------------------------------
        # 3. Sequencing Depth
        # ------
        # ------------------------------------
        plot_scatter_with_repel(x, 'log2_Total', '>x100(%).fragments', f"{seq_run_name}: {panel}", pdf=pdf, label_mask=x['needs_label'])

        # ------------------------------------------
        # 4. % On-Target vs # Reads
        # ------------------------------------------
        plot_scatter_with_repel(x, 'log2_Mapped', 'On Target%', f"{seq_run_name}: {panel}", pdf=pdf, label_mask=x['needs_label'])

        # ------------------------------------------
        # 5. ESTIMATED_LIBRARY_SIZE vs PERCENT_DUPLICATION
        # ------------------------------------------
        plot_scatter_with_repel(x, 'ESTIMATED_LIBRARY_SIZE', 'PERCENT_DUPLICATION', f"{seq_run_name}: {panel}", pdf=pdf, v_lines=[min_lib_size], label_mask=x['needs_label'])

        # ------------------------------------------
        # 6. base > min_base_cov dedup %on-target cov vs ESTIMATED_LIBRARY_SIZE
        # ------------------------------------------
        plot_scatter_with_repel(x, 'ESTIMATED_LIBRARY_SIZE', base_cov_col, f"{seq_run_name}: {panel}", pdf=pdf, v_lines=[min_lib_size], h_lines=[min_f_dedup_cov], label_mask=x['needs_label'])

        # ------------------------------------------
        # 7. base > min_base_cov dedup %on-target cov vs PERCENT_DUPLICATION
        # ------------------------------------------
        plot_scatter_with_repel(x, 'PERCENT_DUPLICATION', base_cov_col, f"{seq_run_name}: {panel}", pdf=pdf, h_lines=[min_f_dedup_cov], label_mask=x['needs_label'])

        # ------------------------------------------
        # 6. Histograms & Filtering
        # ------------------------------------------
        # Histogram
        fig_hist = plt.figure(figsize=(8, 5))
        sns.histplot(data=x, x=base_cov_col, binwidth=0.05)
        plt.axvline(min_f_dedup_cov, color='black', linestyle='--')
        plt.title(f"{seq_run_name}: {panel}")
        plt.tight_layout()
        pdf.savefig(fig_hist)
        plt.close(fig_hist)

        # Jitter plot (Strip plot)
        fig_strip = plt.figure(figsize=(8, 5))
        sns.stripplot(data=x, x='exp', y=base_cov_col, hue='exp', dodge=True, jitter=True, alpha=0.8)
        plt.axhline(min_f_dedup_cov, color='black', linestyle='-')
        plt.title(f"{seq_run_name}: {panel}")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig(fig_strip)
        plt.close(fig_strip)

    print("PDF generation complete.")

    # ------------------------------------------
    # Export Data
    # ------------------------------------------
    all_out = os.path.join(output_dir, f"{panel}_{seq_run_name}_all_classified.csv")
    x.to_csv(all_out, index=False)
    print(f"Exported all samples to {all_out}")

    valid_samples = x[x['Status'] == 'Pass']
    valid_out = os.path.join(output_dir, f"{panel}_{seq_run_name}_valid_samples.csv")
    valid_samples.to_csv(valid_out, index=False)
    print(f"Exported {len(valid_samples)} valid samples to {valid_out}")

    for_reseq = x[x['Status'] != 'Pass']
    reseq_out = os.path.join(output_dir, f"{panel}_{seq_run_name}_to_reseq.csv")
    for_reseq.to_csv(reseq_out, index=False)
    print(f"Exported {len(for_reseq)} samples for re-sequencing to {reseq_out}")


# ==========================================
# Entry Point
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute and plot MIP run statistics for a sequencing run."
    )
    parser.add_argument(
        '-o', '--output-dir',
        default='./MIP_run_statistics',
        help="Directory where output statistics and plots are written.",
    )
    parser.add_argument(
        '-r', '--seq-run-name',
        default='NovaseqR117',
        help="Name of the sequencing run to process.",
    )
    parser.add_argument(
        '-p', '--panel',
        default='arch',
        help="Panel / run type to analyze ('arch' or 'geno').",
    )
    parser.add_argument(
        '--rebuild',
        action='store_true',
        help="Rebuild the combined statistics CSV from per-job Excel files.",
    )
    parser.add_argument(
        '-b', '--runs-base-dir',
        default='/home/projects/shlush/shared/runs_analysis/',
        help="Base directory containing the run analysis folders.",
    )
    parser.add_argument(
        '-m', '--min-f-dedup-cov',
        type=float,
        default=0.95,
        help="Minimum dedup coverage > min_base_cov percentage threshold for valid samples.",
    )
    parser.add_argument(
        '-s', '--sample-fixed-column',
        default='FixedSample',
        help="Column name in the statistics file that contains the fixed sample names.",
    )
    parser.add_argument(
        '--min-base-cov',
        type=int,
        default=20,
        help="Minimum per-base dedup coverage depth for the dedup coverage check.",
    )
    parser.add_argument(
        '--min-lib-size',
        type=int,
        default=70000,
        help="Minimum estimated library size for valid samples.",
    )
    parser.add_argument(
        '--fastqc-dir',
        default=None,
        help="Directory with FastQC summary tsv files (fastqc_summary_READ.tsv). "
             "Defaults to <runs-base-dir>/<seq-run-name>/fastq_PERIBLOOD/Reports/FastQC.",
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    main(
        output_dir=args.output_dir,
        seq_run_name=args.seq_run_name,
        panel=args.panel,
        rebuild=args.rebuild,
        runs_base_dir=args.runs_base_dir,
        min_f_dedup_cov=args.min_f_dedup_cov,
        sample_fixed_column=args.sample_fixed_column,
        min_base_cov=args.min_base_cov,
        min_lib_size=args.min_lib_size,
        fastqc_dir=args.fastqc_dir,
    )