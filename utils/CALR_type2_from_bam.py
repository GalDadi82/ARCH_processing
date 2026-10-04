#!/usr/bin/env python3
"""
CALR Type-2 (ins5) detector — BAM-based, pileup window approach.
- Scans each sample BAM for 5 bp insertions matching a known motif.
- Default: submits one LSF job per sample (parallel).
- Worker mode: processes a single sample on the compute node.

Submit mode (default — queues one LSF job per sample):
  python CALR_type2_from_bam.py --bam_dir /path/to/samples

Worker mode (called automatically by each LSF job):
  python CALR_type2_from_bam.py --sample_name SAMPLE --bam_path /path/to/sample.bam --out out.tsv

Merge results once all jobs finish:
  awk '(NR==1)||(FNR>1)' calr_results/*_CALR_type2.tsv > CALR_type2_ins5_fromBAM_merged.tsv
"""

import os, re, sys, glob, argparse, subprocess
from collections import Counter
from typing import Optional
import pandas as pd
from tqdm import tqdm

try:
    import pysam
except ImportError:
    sys.stderr.write("pysam is required. Install with:\npython -m pip install --user pysam\n")
    sys.exit(1)

# ------------------ DEFAULTS ------------------
DEFAULT_BAM_DIR   = "/home/projects/shlush/shared/runs_analysis/MiseqR3/analysis/3_alignment_to_annotation_bwa"
DEFAULT_OUT_TSV   = "CALR_type2_ins5_fromBAM.tsv"

# IMPORTANT: BAM shows chr19 reads around 12943711 (not 13054550)
DEFAULT_CHROM     = "chr19"
DEFAULT_WIN_START = 12943000   # 1-based
DEFAULT_WIN_END   = 12945000   # 1-based

DEFAULT_INS_LEN   = 5
DEFAULT_MOTIF_FWD = "TTGTC"
DEFAULT_MOTIF_REV = "GACAA"    # reverse complement

DEFAULT_MIN_MAPQ  = 0
DEFAULT_SKIP_DUP  = False

DEFAULT_QUEUE   = "short"
DEFAULT_MEM_MB  = 32768
# ----------------------------------------------

try:
    sys.stderr.reconfigure(line_buffering=True)
except Exception:
    pass


def read_overlaps_window(read, start0, end0):
    if read.reference_start is None or read.reference_end is None:
        return False
    return not (read.reference_end <= start0 or read.reference_start >= end0)


def has_type2_insertion(
    read,
    ins_len: int = DEFAULT_INS_LEN,
    motif_fwd: str = DEFAULT_MOTIF_FWD,
    motif_rev: str = DEFAULT_MOTIF_REV,
) -> Optional[int]:
    """
    Detect ins_len bp insertion matching motif_fwd or motif_rev.
    read.query_sequence is stored as-sequenced; both strands are checked.
    Returns the 1-based reference position of the base immediately preceding
    the insertion (VCF-style anchor), or None if no match.
    """
    if read.cigartuples is None or read.query_sequence is None:
        return None

    qpos = 0
    rpos = read.reference_start  # 0-based
    seq = read.query_sequence.upper()

    for op, ln in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X: consume query + ref
            qpos += ln
            rpos += ln
        elif op == 1:  # I: consume query only
            ins_seq = seq[qpos:qpos + ln]
            if ln == ins_len and ins_seq in (motif_fwd, motif_rev):
                return rpos
            qpos += ln
        elif op in (2, 3):  # D, N: consume ref only
            rpos += ln
        elif op == 4:  # S consumes query
            qpos += ln
        # H, P: do not consume query or ref

    return None


def pick_bam(sample_dir):
    """Pick preferred BAM, else newest .bam (exclude .bai)."""
    preferred = [p for p in glob.glob(os.path.join(sample_dir, "*.sorted.rg*.bam")) if not p.endswith(".bai")]
    if preferred:
        preferred.sort()
        return preferred[0]

    all_bams = [p for p in glob.glob(os.path.join(sample_dir, "*.bam")) if not p.endswith(".bai")]
    if not all_bams:
        return None

    all_bams.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return all_bams[0]


def resolve_chrom(bam, chrom_pref):
    """Return the actual contig name present in BAM (chr19 vs 19)."""
    refs = set(bam.references)
    if chrom_pref in refs:
        return chrom_pref
    if chrom_pref.startswith("chr"):
        alt = chrom_pref[3:]
        if alt in refs:
            return alt
    else:
        alt = "chr" + chrom_pref
        if alt in refs:
            return alt
    return None


def count_type2_in_bam(
    bam_path: str,
    chrom: str = DEFAULT_CHROM,
    win_start: int = DEFAULT_WIN_START,
    win_end: int = DEFAULT_WIN_END,
    ins_len: int = DEFAULT_INS_LEN,
    motif_fwd: str = DEFAULT_MOTIF_FWD,
    motif_rev: str = DEFAULT_MOTIF_REV,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    skip_dup: bool = DEFAULT_SKIP_DUP,
):
    # convert to 0-based half-open for pysam.fetch
    s0 = win_start - 1
    e0 = win_end

    depth = 0
    ins5 = 0
    ins_positions = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_used = resolve_chrom(bam, chrom)
        if chrom_used is None:
            raise ValueError(f"Contig not found in BAM: {chrom} (available: {bam.references[:5]})")

        ncheck = sum(1 for _ in zip(range(5), bam.fetch(chrom_used, s0, e0)))
        if ncheck == 0:
            sys.stderr.write(f"WARNING: no reads in window for {os.path.basename(bam_path)} at {chrom_used}:{win_start}-{win_end}\n")

        for read in bam.fetch(chrom_used, s0, e0):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if skip_dup and read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if not read_overlaps_window(read, s0, e0):
                continue

            depth += 1
            pos = has_type2_insertion(read, ins_len, motif_fwd, motif_rev)
            if pos is not None:
                ins5 += 1
                ins_positions.append(pos)

    wt = max(0, depth - ins5)
    vaf = ins5 / depth if depth > 0 else 0.0

    if ins_positions:
        pos_counts = Counter(ins_positions)
        ins_pos, pos_support = pos_counts.most_common(1)[0]
        pos_frac = pos_support / len(ins_positions)
    else:
        ins_pos, pos_frac = None, 0.0

    return chrom_used, depth, ins5, wt, round(vaf, 6), ins_pos, round(pos_frac, 6)


def submit_jobs(
    bam_dir: str,
    results_dir: str = "calr_results",
    chrom: str = DEFAULT_CHROM,
    win_start: int = DEFAULT_WIN_START,
    win_end: int = DEFAULT_WIN_END,
    ins_len: int = DEFAULT_INS_LEN,
    motif_fwd: str = DEFAULT_MOTIF_FWD,
    motif_rev: str = DEFAULT_MOTIF_REV,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    skip_dup: bool = DEFAULT_SKIP_DUP,
    queue: str = DEFAULT_QUEUE,
    mem_mb: int = DEFAULT_MEM_MB,
    sample_pattern: str = "ARCH",
):
    """Submit mode: discover sample BAMs and queue one LSF job per sample."""
    if not os.path.isdir(bam_dir):
        sys.stderr.write(f"ERROR: bam_dir not found: {bam_dir}\n")
        sys.exit(1)

    samples = sorted(
        s for s in os.listdir(bam_dir)
        if os.path.isdir(os.path.join(bam_dir, s))
    )
    if not samples:
        sys.stderr.write(f"ERROR: no sample sub-directories in {bam_dir}\n")
        sys.exit(1)

    if sample_pattern:
        before = len(samples)
        samples = [s for s in samples if re.search(sample_pattern, s)]
        print(f"Filtered to {len(samples)}/{before} samples matching pattern '{sample_pattern}'")
        if not samples:
            sys.stderr.write(f"ERROR: no samples match pattern '{sample_pattern}'\n")
            sys.exit(1)

    os.makedirs("calr_logs", exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    script_path = os.path.abspath(__file__)
    python_exe  = sys.executable

    submitted = 0
    for sample in tqdm(samples, desc="Submitting samples"):
        bam = pick_bam(os.path.join(bam_dir, sample))
        if not bam:
            sys.stderr.write(f"WARNING: no BAM found for {sample}, skipping\n")
            continue

        out_tsv = os.path.join(results_dir, f"{sample}_CALR_type2.tsv")
        bsub_cmd = [
            "bsub",
            "-q", queue,
            "-M", str(mem_mb),
            "-R", f"rusage[mem={mem_mb}]",
            "-o", f"calr_logs/{sample}_type2.log",
            "-e", f"calr_logs/{sample}_type2.err",
            python_exe, script_path,
            "--sample_name", sample,
            "--bam_path",    bam,
            "--out",         out_tsv,
            "--chrom",       chrom,
            "--win_start",   str(win_start),
            "--win_end",     str(win_end),
            "--ins_len",     str(ins_len),
            "--motif_fwd",   motif_fwd,
            "--motif_rev",   motif_rev,
            "--min_mapq",    str(min_mapq),
        ]
        if skip_dup:
            bsub_cmd.append("--skip_dup")

        subprocess.run(bsub_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        submitted += 1

    print(f"Submitted {submitted} jobs to '{queue}' queue.")
    print("To merge results once done:")
    print(f"  awk '(NR==1)||(FNR>1)' {results_dir}/*_CALR_type2.tsv > CALR_type2_ins5_fromBAM_merged.tsv")


def run_worker(
    sample: str,
    bam_path: str,
    out_tsv: str,
    chrom: str = DEFAULT_CHROM,
    win_start: int = DEFAULT_WIN_START,
    win_end: int = DEFAULT_WIN_END,
    ins_len: int = DEFAULT_INS_LEN,
    motif_fwd: str = DEFAULT_MOTIF_FWD,
    motif_rev: str = DEFAULT_MOTIF_REV,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    skip_dup: bool = DEFAULT_SKIP_DUP,
):
    """Worker mode: process a single sample BAM on the compute node."""
    try:
        chrom_used, depth, ins5, wt, vaf, ins_pos, pos_frac = count_type2_in_bam(
            bam_path,
            chrom=chrom,
            win_start=win_start,
            win_end=win_end,
            ins_len=ins_len,
            motif_fwd=motif_fwd,
            motif_rev=motif_rev,
            min_mapq=min_mapq,
            skip_dup=skip_dup,
        )
        row = {
            "Sample_Name":        sample,
            "bam_used":           os.path.basename(bam_path),
            "CHR":         chrom_used,
            "POS":         ins_pos if ins_pos is not None else "",
            "POS_frac":    pos_frac,
            "Gene.refGene": "CALR",
            "Selected.AA.change": "Type-2 (ins5)",
            "window_start_1based": win_start,
            "window_end_1based":   win_end,
            "depth_window":        depth,
            "INS5_count":          ins5,
            "WT_count":            wt,
            "VAF_INS5":            vaf,
            "avg_VAF":             vaf,  # same as VAF_INS5 for single window
            "status":              "OK" if depth > 0 else "ZERO_DEPTH_WINDOW_WRONG",
        }
    except Exception as e:
        sys.stderr.write(f"ERROR {sample}: {e}\n")
        row = {
            "Sample_Name":        sample,
            "bam_used":           os.path.basename(bam_path),
            "CHR":         "",
            "POS":         "",
            "POS_frac":    0.0,
            "Gene.refGene": "",
            "Selected.AA.change": "",
            "window_start_1based": win_start,
            "window_end_1based":   win_end,
            "depth_window":        0,
            "INS5_count":          0,
            "WT_count":            0,
            "VAF_INS5":            0.0,
            "avg_VAF":             0.0,
            "status":              f"ERROR: {e}",
        }

    df = pd.DataFrame([row], columns=[
        "Sample_Name", "bam_used", "CHR", "POS", "POS_frac", "Gene.refGene", "Selected.AA.change",
        "window_start_1based", "window_end_1based",
        "depth_window", "INS5_count", "WT_count", "VAF_INS5", "avg_VAF", "status",
    ])
    df.to_csv(out_tsv, sep="\t", index=False)
    sys.stderr.write(f"Finished {sample} -> {out_tsv}\n")


def main():
    parser = argparse.ArgumentParser(
        description="CALR Type-2 (ins5) detector — BAM-based, LSF parallel mode.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- mode (mutually exclusive; default is submit) ---
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--sample_name", metavar="NAME",
                      help="Worker mode: sample name (supplied by LSF job)")
    mode.add_argument("--bam_dir", default=DEFAULT_BAM_DIR,
                      help="Submit mode: root directory of per-sample BAM sub-folders")

    # --- submit-mode LSF options ---
    parser.add_argument("--results_dir", default="calr_results",
                        help="Directory for per-sample output TSVs")
    parser.add_argument("--queue", default=DEFAULT_QUEUE,
                        help="LSF queue name")
    parser.add_argument("--mem_mb", type=int, default=DEFAULT_MEM_MB,
                        help="Memory limit per LSF job in MB")
    parser.add_argument("--sample_pattern", default="ARCH",
                        help="Regex (grep-like) filter applied to sample folder names discovered in bam_dir; only matching samples are submitted")

    # --- worker-mode inputs ---
    parser.add_argument("--bam_path", metavar="FILE.bam",
                        help="Worker mode: path to the sample BAM")
    parser.add_argument("--out", metavar="OUT.tsv",
                        help="Worker mode: output TSV path")

    # --- genomic window (shared) ---
    parser.add_argument("--chrom", default=DEFAULT_CHROM,
                        help="Chromosome / contig name (chr-prefix auto-resolved)")
    parser.add_argument("--win_start", type=int, default=DEFAULT_WIN_START,
                        help="Window start position (1-based)")
    parser.add_argument("--win_end", type=int, default=DEFAULT_WIN_END,
                        help="Window end position (1-based)")

    # --- insertion params (shared) ---
    parser.add_argument("--ins_len", type=int, default=DEFAULT_INS_LEN,
                        help="Expected insertion length in bp")
    parser.add_argument("--motif_fwd", default=DEFAULT_MOTIF_FWD,
                        help="Forward-strand insertion motif")
    parser.add_argument("--motif_rev", default=DEFAULT_MOTIF_REV,
                        help="Reverse-strand insertion motif (RC of motif_fwd)")

    # --- read filters (shared) ---
    parser.add_argument("--min_mapq", type=int, default=DEFAULT_MIN_MAPQ,
                        help="Minimum mapping quality to include a read")
    parser.add_argument("--skip_dup", action="store_true", default=DEFAULT_SKIP_DUP,
                        help="Skip reads flagged as PCR duplicates")

    args = parser.parse_args()

    algo = dict(
        chrom=args.chrom, win_start=args.win_start, win_end=args.win_end,
        ins_len=args.ins_len, motif_fwd=args.motif_fwd, motif_rev=args.motif_rev,
        min_mapq=args.min_mapq, skip_dup=args.skip_dup,
    )

    if args.sample_name:
        if not (args.bam_path and args.out):
            parser.error("Worker mode requires --bam_path and --out")
        run_worker(args.sample_name, args.bam_path, args.out, **algo)
    else:
        submit_jobs(
            bam_dir=args.bam_dir,
            results_dir=args.results_dir,
            queue=args.queue,
            mem_mb=args.mem_mb,
            sample_pattern=args.sample_pattern,
            **algo,
        )


if __name__ == "__main__":
    main()
