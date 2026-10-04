"""
CALR Type-1 (del52) detector — FASTQ-only, breakpoint (junction) based.
- Learns a panel of del52 junction motifs from a seed WT R1 sequence.
- Submits an LSF job per sample to process FASTQ pairs in parallel.

Submit mode (queues one LSF job per sample):
  python CALR_type1_from_fastq.py --submit --fastq_dir /path/to/fastqs

Worker mode (called automatically by each LSF job):
  python CALR_type1_from_fastq.py --sample_name SAMPLE --fq1 R1.fastq.gz --fq2 R2.fastq.gz --out out.tsv

Results land in calr_results/<sample>_CALR.tsv; merge with:
  awk '(NR==1)||(FNR>1)' calr_results/*_CALR.tsv > CALR_type1_fromFASTQ_merged.tsv
"""

import os, re, sys, gzip, subprocess, argparse
from typing import List, Tuple, Set
import pandas as pd
from tqdm import tqdm

# ------------------------- DEFAULTS -------------------------
DEFAULT_FASTQ_DIR = "/home/projects/shlush/shared/runs_analysis/PERIBLOOD_merged_arch/fastq_PERIBLOOD"
DEFAULT_SEED_WT_R1 = (
    "GCAGCAGAGAAACAAATGAAGGACAAACAGGACGAGGAGCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAGGCAGAGGACAAGGAGGAAGATGAGGACAAAGATGAGGATGAGGAGGATGAGGAGGACAAGG"
)
DEFAULT_DEL_LEN   = 52
DEFAULT_LEFT_LEN  = 12
DEFAULT_RIGHT_LEN = 12
DEFAULT_WT_TOL    = 6
DEFAULT_MIN_READLEN = 50

# ----------------------------------------------------------
def openq(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

PAIR_RE = re.compile(r"^(?P<base>.+)_R(?P<mate>[12]).*\.fastq\.gz$")

def find_fastq_pairs(folder: str):
    files = [f for f in os.listdir(folder) if f.endswith(".fastq.gz")]
    pairs = {}
    for f in files:
        m = PAIR_RE.match(f)
        if not m:
            continue
        base = m.group("base")
        mate = m.group("mate")
        entry = pairs.setdefault(base, {"R1": None, "R2": None})
        full = os.path.join(folder, f)
        if mate == "1":
            entry["R1"] = full
        elif mate == "2":
            entry["R2"] = full
    out = []
    for base, d in sorted(pairs.items()):
        if d["R1"] and d["R2"]:
            out.append((base, d["R1"], d["R2"]))
    return out

def build_junction_panel(seed: str, left_len: int, right_len: int, del_len: int) -> Tuple[Set[str], List[Tuple[str,str]]]:
    seed = seed.upper()
    J: Set[str] = set()
    PAIRS: List[Tuple[str,str]] = []
    max_i = len(seed) - (left_len + del_len + right_len)
    for i in range(0, max_i + 1):
        left = seed[i:i+left_len]
        right_start = i + left_len + del_len
        right = seed[right_start:right_start + right_len]
        if len(left) != left_len or len(right) != right_len:
            continue
        if not (set(left) <= set("ACGT") and set(right) <= set("ACGT")):
            continue
        PAIRS.append((left, right))
        J.add(left + right)
    return J, PAIRS

def count_sample(
    fq1: str,
    fq2: str,
    junctions: Set[str],
    lr_pairs: List[Tuple[str, str]],
    del_len: int = DEFAULT_DEL_LEN,
    wt_tol: int = DEFAULT_WT_TOL,
    min_readlen: int = DEFAULT_MIN_READLEN,
) -> Tuple[int, int]:
    WT = 0
    MUT = 0
    low = del_len - wt_tol
    high = del_len + wt_tol
    with openq(fq1) as f1, openq(fq2) as f2:
        while True:
            h1 = f1.readline(); h2 = f2.readline()
            if not h1 or not h2:
                break
            r1 = f1.readline().strip().upper(); r2 = f2.readline().strip().upper()
            f1.readline(); f1.readline()
            f2.readline(); f2.readline()

            classified = False
            for seq in (r1, r2):
                if len(seq) < min_readlen:
                    continue

                hit_mut = False
                for j in junctions:
                    if j in seq:
                        MUT += 1
                        hit_mut = True
                        classified = True
                        break
                if hit_mut:
                    break

                step = max(1, len(lr_pairs) // 100)
                for idx in range(0, len(lr_pairs), step):
                    L, R = lr_pairs[idx]
                    pL = seq.find(L)
                    if pL == -1:
                        continue
                    start = pL + len(L) + low
                    end   = pL + len(L) + high + len(R)
                    if start < 0: start = 0
                    if end > len(seq): end = len(seq)
                    segment = seq[start:end]
                    if R in segment:
                        WT += 1
                        classified = True
                        break
                if classified:
                    break
    return WT, MUT

def submit_jobs(
    fastq_dir: str,
    seed_wt: str = DEFAULT_SEED_WT_R1,
    del_len: int = DEFAULT_DEL_LEN,
    left_len: int = DEFAULT_LEFT_LEN,
    right_len: int = DEFAULT_RIGHT_LEN,
    wt_tol: int = DEFAULT_WT_TOL,
    min_readlen: int = DEFAULT_MIN_READLEN,
    queue: str = "short",
    mem_mb: int = 32768,
    sample_pattern: str = "ARCH",
):
    """Master mode: find all FASTQ pairs and submit one LSF job per sample."""
    if not os.path.isdir(fastq_dir):
        print(f"ERROR: fastq_dir not found: {fastq_dir}", file=sys.stderr)
        sys.exit(1)

    pairs = find_fastq_pairs(fastq_dir)
    if not pairs:
        print(f"ERROR: no paired FASTQs found in {fastq_dir}", file=sys.stderr)
        sys.exit(1)

    if sample_pattern:
        before = len(pairs)
        pairs = [p for p in pairs if re.search(sample_pattern, p[0])]
        print(f"Filtered to {len(pairs)}/{before} samples matching pattern '{sample_pattern}'")
        if not pairs:
            print(f"ERROR: no samples match pattern '{sample_pattern}'", file=sys.stderr)
            sys.exit(1)

    os.makedirs("calr_logs", exist_ok=True)
    os.makedirs("calr_results", exist_ok=True)

    script_path = os.path.abspath(__file__)
    python_exe = sys.executable

    print(f"Submitting {len(pairs)} jobs to the '{queue}' queue...")

    for base, fq1, fq2 in tqdm(pairs, desc="Submitting samples"):
        out_tsv = f"calr_results/{base}_CALR.tsv"
        bsub_cmd = [
            "bsub",
            "-q", queue,
            "-M", str(mem_mb),
            "-R", f"rusage[mem={mem_mb}]",
            "-o", f"calr_logs/{base}_type1.log",
            "-e", f"calr_logs/{base}_type1.err",
            python_exe, script_path,
            "--sample_name", base,
            "--fq1", fq1,
            "--fq2", fq2,
            "--out", out_tsv,
            "--seed_wt", seed_wt,
            "--del_len", str(del_len),
            "--left_len", str(left_len),
            "--right_len", str(right_len),
            "--wt_tol", str(wt_tol),
            "--min_readlen", str(min_readlen),
        ]
        subprocess.run(bsub_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"All {len(pairs)} jobs submitted.")
    print("To merge results once done:")
    print("  awk '(NR==1)||(FNR>1)' calr_results/*_CALR.tsv > CALR_type1_fromFASTQ_merged.tsv")


def run_worker(
    base: str,
    fq1: str,
    fq2: str,
    out_tsv: str,
    seed_wt: str = DEFAULT_SEED_WT_R1,
    del_len: int = DEFAULT_DEL_LEN,
    left_len: int = DEFAULT_LEFT_LEN,
    right_len: int = DEFAULT_RIGHT_LEN,
    wt_tol: int = DEFAULT_WT_TOL,
    min_readlen: int = DEFAULT_MIN_READLEN,
):
    """Worker mode: process a single sample on the compute node."""
    junctions, lr_pairs = build_junction_panel(seed_wt, left_len, right_len, del_len)
    if not junctions:
        print("ERROR: failed to derive junction panel from seed WT sequence.", file=sys.stderr)
        sys.exit(1)

    wt, mut = count_sample(fq1, fq2, junctions, lr_pairs, del_len, wt_tol, min_readlen)
    vaf = mut / (mut + wt) if (mut + wt) > 0 else 0.0

    df = pd.DataFrame([{
        "Sample_Name": base,
        "Gene.refGene": "CALR",
        "Selected.AA.change": "Type-1 (del52)",
        "CHR": "19",
        "POS": "",
        "WT_count": wt,
        "DEL52_count": mut,
        "VAF_DEL52": round(vaf, 6),
        "avg_VAF": round(vaf, 6),
    }])
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Finished {base} -> {out_tsv}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="CALR Type-1 (del52) detector — FASTQ-only, LSF parallel mode.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- mode ---
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--submit", action="store_true",
                      help="Submit mode: find all FASTQ pairs and queue LSF jobs")
    mode.add_argument("--sample_name", metavar="NAME",
                      help="Worker mode: sample name (set by LSF job)")

    # --- submit-mode inputs ---
    parser.add_argument("--fastq_dir", default=DEFAULT_FASTQ_DIR,
                        help="Directory containing *_R{1,2}.fastq.gz files")
    parser.add_argument("--queue", default="short",
                        help="LSF queue name")
    parser.add_argument("--mem_mb", type=int, default=32768,
                        help="Memory limit per job in MB")
    parser.add_argument("--sample_pattern", default="ARCH",
                        help="Regex (grep-like) filter applied to sample/base names discovered in fastq_dir; only matching samples are submitted")

    # --- worker-mode inputs ---
    parser.add_argument("--fq1", metavar="R1.fastq.gz",
                        help="Worker mode: path to R1 FASTQ")
    parser.add_argument("--fq2", metavar="R2.fastq.gz",
                        help="Worker mode: path to R2 FASTQ")
    parser.add_argument("--out", metavar="OUT.tsv",
                        help="Worker mode: output TSV path")

    # --- shared algorithm params ---
    parser.add_argument("--seed_wt", default=DEFAULT_SEED_WT_R1,
                        help="Seed WT R1 sequence used to build the junction panel")
    parser.add_argument("--del_len", type=int, default=DEFAULT_DEL_LEN,
                        help="Expected deletion length in bp")
    parser.add_argument("--left_len", type=int, default=DEFAULT_LEFT_LEN,
                        help="Left flank k-mer length")
    parser.add_argument("--right_len", type=int, default=DEFAULT_RIGHT_LEN,
                        help="Right flank k-mer length")
    parser.add_argument("--wt_tol", type=int, default=DEFAULT_WT_TOL,
                        help="Tolerance around del_len for WT classification (± bp)")
    parser.add_argument("--min_readlen", type=int, default=DEFAULT_MIN_READLEN,
                        help="Minimum read length to classify")

    args = parser.parse_args()

    if args.submit:
        submit_jobs(
            fastq_dir=args.fastq_dir,
            seed_wt=args.seed_wt,
            del_len=args.del_len,
            left_len=args.left_len,
            right_len=args.right_len,
            wt_tol=args.wt_tol,
            min_readlen=args.min_readlen,
            queue=args.queue,
            mem_mb=args.mem_mb,
            sample_pattern=args.sample_pattern,
        )
    else:
        if not (args.fq1 and args.fq2 and args.out):
            parser.error("Worker mode requires --fq1, --fq2, and --out")
        run_worker(
            base=args.sample_name,
            fq1=args.fq1,
            fq2=args.fq2,
            out_tsv=args.out,
            seed_wt=args.seed_wt,
            del_len=args.del_len,
            left_len=args.left_len,
            right_len=args.right_len,
            wt_tol=args.wt_tol,
            min_readlen=args.min_readlen,
        )


if __name__ == "__main__":
    main()
