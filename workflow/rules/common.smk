# Shared paths, helpers and the sample-list plumbing for the ARCH pipeline.
#
# Every path the notebook defined in its parameters cell is derived here, so the rest of
# the workflow refers to names rather than to f-strings.

import glob
import os
import re

import pandas as pd


# ----------------------------------------------------------------------- basics
SEQ_RUN = config["seq_run"]
PANEL = config["panel"]
RUNS_BASE_DIR = config["runs_base_dir"]

# utils/mip_run_statistics.py derives the analysis dir as analysis_{panel}; keep in step.
ANALYSIS_DIR = os.path.join(RUNS_BASE_DIR, SEQ_RUN, f"analysis_{PANEL}")
OUTPUT_BASE_DIR = os.path.join(ANALYSIS_DIR, "All_Samples")
VCFS_DIR = os.path.join(OUTPUT_BASE_DIR, "5_vcfs")
MAPPED_DIR = os.path.join(OUTPUT_BASE_DIR, "2_mapped")
FASTQ_DIR = os.path.join(RUNS_BASE_DIR, SEQ_RUN, config["fastq_subdir"])

# Repo root: workflow/ lives directly under it.
REPO_ROOT = os.path.abspath(os.path.join(workflow.basedir, ".."))
UTILS_DIR = os.path.join(REPO_ROOT, "utils")


def from_repo(path):
    """Resolve a config path relative to the repo root unless it is already absolute."""
    return path if os.path.isabs(path) else os.path.join(REPO_ROOT, path)


# ------------------------------------------------------------------- QC outputs
LINK_MANIFEST = os.path.join(OUTPUT_BASE_DIR, "sample_links.tsv")
ALL_CLASSIFIED = os.path.join(VCFS_DIR, f"{PANEL}_{SEQ_RUN}_all_classified.csv")
SAMPLES_FOR_PROCESSING = os.path.join(VCFS_DIR, f"{PANEL}_{SEQ_RUN}_samples_for_processing.csv")

# ------------------------------------------------------------- merged / filtered
MERGED_VARSCAN = os.path.join(VCFS_DIR, "merged_varscan_annotations.txt")
MULTIANNO_REPORT = os.path.join(VCFS_DIR, "multianno_file_report.xlsx")

FILTER_RUN_NAME = f"{SEQ_RUN}_{config['filter_mutations']['run_suffix']}"
LISTED = os.path.join(VCFS_DIR, f"FilterMutations_{FILTER_RUN_NAME}_listed.tsv")
NOT_LISTED = os.path.join(VCFS_DIR, f"FilterMutations_{FILTER_RUN_NAME}_filtered_yet_not_listed.tsv")

# "_FilteredRecurrent" drops out of the names when the recurrence step is skipped, so a
# run's filtered and unfiltered outputs can coexist on disk.
RECURRENCE_TAG = "_FilteredRecurrent" if config["filter_recurrent_without_support"] else ""
OUT_PREFIX = os.path.join(VCFS_DIR, f"FilterMutations_{SEQ_RUN}_{config['output_suffix']}")

NO_RECURRENT = f"{OUT_PREFIX}_listed{RECURRENCE_TAG}.tsv"
MASTER = f"{OUT_PREFIX}_listed{RECURRENCE_TAG}_VAF_{config['min_final_vaf']}.tsv"

IS_IN_DENVER = MASTER.replace(".tsv", "_is_in_denver.tsv")
MASTER_WITH_DENVER = MASTER.replace(".tsv", "_with_denver.tsv")
DENVER_FILTERED = MASTER.replace(".tsv", "_Denver_filtered.tsv")

# The alert reads the Denver-filtered table when Denver is enabled, else the master table.
FINAL_MUTATIONS = DENVER_FILTERED if config["denver"]["enabled"] else MASTER

# ------------------------------------------------------------------------- CALR
CALR_DIR = os.path.join(VCFS_DIR, "calr_results")
CALR_TYPE1_RESULTS = os.path.join(VCFS_DIR, f"CALR_type1_results_{SEQ_RUN}.tsv")
CALR_TYPE2_RESULTS = os.path.join(VCFS_DIR, f"CALR_type2_results_{SEQ_RUN}.tsv")

QC_PLOTS = os.path.join(VCFS_DIR, f"{PANEL}_{SEQ_RUN}_mutation_plots.pdf")


# --------------------------------------------------------------- sample listing
def read_processing_samples():
    """Sample names selected for processing. Only valid after the checkpoint has run."""
    df = pd.read_csv(checkpoints.select_samples.get().output[0])
    return sorted(df["Sample"].astype(str).unique())


def calr_samples():
    """Samples that get CALR calls -- the util scripts filter by name pattern too."""
    pattern = config["calr"]["sample_pattern"]
    samples = read_processing_samples()
    if pattern:
        samples = [s for s in samples if re.search(pattern, s)]
    return samples


# ------------------------------------------------------- per-sample input paths
# Mirrors PAIR_RE in utils/CALR_type1_from_fastq.py.
_PAIR_RE = re.compile(r"^(?P<base>.+)_R(?P<mate>[12]).*\.fastq\.gz$")


def find_fastq_pair(sample):
    """Return (R1, R2) for a sample, mirroring find_fastq_pairs() in the CALR type-1 util."""
    r1 = r2 = None
    for fn in sorted(os.listdir(FASTQ_DIR)):
        m = _PAIR_RE.match(fn)
        if not m or m.group("base") != sample:
            continue
        if m.group("mate") == "1":
            r1 = os.path.join(FASTQ_DIR, fn)
        else:
            r2 = os.path.join(FASTQ_DIR, fn)
    if not (r1 and r2):
        raise FileNotFoundError(f"No complete FASTQ pair for {sample} under {FASTQ_DIR}")
    return r1, r2


def pick_bam(sample):
    """Preferred BAM for a sample, else the newest .bam -- mirrors the CALR type-2 util."""
    sample_dir = os.path.join(MAPPED_DIR, sample)
    preferred = sorted(
        p for p in glob.glob(os.path.join(sample_dir, "*.sorted.rg*.bam"))
        if not p.endswith(".bai")
    )
    if preferred:
        return preferred[0]
    all_bams = [p for p in glob.glob(os.path.join(sample_dir, "*.bam")) if not p.endswith(".bai")]
    if not all_bams:
        raise FileNotFoundError(f"No BAM found for {sample} under {sample_dir}")
    return max(all_bams, key=os.path.getmtime)


def sample_vcf(wildcards):
    template = config["annovar"]["vcf_template"]
    return os.path.join(VCFS_DIR, wildcards.sample, template.format(sample=wildcards.sample))


# ANNOVAR names its table {out}.{buildver}_multianno.txt.
MULTIANNO_STEM = f"{config['annovar']['out_basename']}.{config['annovar']['buildver']}_multianno"
MULTIANNO_TXT = os.path.join(VCFS_DIR, "{sample}", MULTIANNO_STEM + ".txt")
MULTIANNO_VCF = os.path.join(VCFS_DIR, "{sample}", MULTIANNO_STEM + ".vcf")
MULTIANNO_FORMATTED = os.path.join(VCFS_DIR, "{sample}", MULTIANNO_STEM + ".formatted.txt")


# ------------------------------------------------------- checkpoint aggregation
def multianno_formatted(wildcards):
    """All per-sample formatted ANNOVAR tables -- forces the checkpoint to resolve first."""
    return expand(MULTIANNO_FORMATTED, sample=read_processing_samples())


def multianno_raw(wildcards):
    return expand(MULTIANNO_TXT, sample=read_processing_samples())


def calr_type1_files(wildcards):
    return expand(os.path.join(CALR_DIR, "{sample}_CALR.tsv"), sample=calr_samples())


def calr_type2_files(wildcards):
    return expand(os.path.join(CALR_DIR, "{sample}_CALR_type2.tsv"), sample=calr_samples())
