# Stage 4: CALR Type-1 (del52, from FASTQ) and Type-2 (ins5, from BAM) detection.
#
# The util scripts have a --submit mode that bsubs one job per sample; here we call their
# worker mode directly and let Snakemake do the scatter.


rule calr_type1:
    """Detect the CALR Type-1 52bp deletion for one sample, from its FASTQ pair."""
    input:
        samples=SAMPLES_FOR_PROCESSING,
    output:
        tsv=os.path.join(CALR_DIR, "{sample}_CALR.tsv"),
    params:
        script=os.path.join(UTILS_DIR, "CALR_type1_from_fastq.py"),
        fq=lambda w: find_fastq_pair(w.sample),
        seed_wt=config["calr"]["type1"]["seed_wt"],
        del_len=config["calr"]["type1"]["del_len"],
        left_len=config["calr"]["type1"]["left_len"],
        right_len=config["calr"]["type1"]["right_len"],
        wt_tol=config["calr"]["type1"]["wt_tol"],
        min_readlen=config["calr"]["type1"]["min_readlen"],
    log:
        os.path.join(VCFS_DIR, "logs", "calr_type1", "{sample}.log"),
    resources:
        mem_mb=32768,
        runtime=240,
    shell:
        r"""
        python {params.script:q} \
            --sample_name {wildcards.sample:q} \
            --fq1 {params.fq[0]:q} \
            --fq2 {params.fq[1]:q} \
            --out {output.tsv:q} \
            --seed_wt {params.seed_wt:q} \
            --del_len {params.del_len} \
            --left_len {params.left_len} \
            --right_len {params.right_len} \
            --wt_tol {params.wt_tol} \
            --min_readlen {params.min_readlen} > {log:q} 2>&1
        """


rule calr_type2:
    """Detect the CALR Type-2 5bp insertion for one sample, from its BAM."""
    input:
        samples=SAMPLES_FOR_PROCESSING,
    output:
        tsv=os.path.join(CALR_DIR, "{sample}_CALR_type2.tsv"),
    params:
        script=os.path.join(UTILS_DIR, "CALR_type2_from_bam.py"),
        bam=lambda w: pick_bam(w.sample),
        chrom=config["calr"]["type2"]["chrom"],
        win_start=config["calr"]["type2"]["win_start"],
        win_end=config["calr"]["type2"]["win_end"],
        ins_len=config["calr"]["type2"]["ins_len"],
        motif_fwd=config["calr"]["type2"]["motif_fwd"],
        motif_rev=config["calr"]["type2"]["motif_rev"],
        min_mapq=config["calr"]["type2"]["min_mapq"],
        skip_dup="--skip_dup" if config["calr"]["type2"]["skip_dup"] else "",
    log:
        os.path.join(VCFS_DIR, "logs", "calr_type2", "{sample}.log"),
    resources:
        mem_mb=32768,
        runtime=240,
    shell:
        r"""
        python {params.script:q} \
            --sample_name {wildcards.sample:q} \
            --bam_path {params.bam:q} \
            --out {output.tsv:q} \
            --chrom {params.chrom:q} \
            --win_start {params.win_start} \
            --win_end {params.win_end} \
            --ins_len {params.ins_len} \
            --motif_fwd {params.motif_fwd:q} \
            --motif_rev {params.motif_rev:q} \
            --min_mapq {params.min_mapq} \
            {params.skip_dup} > {log:q} 2>&1
        """


rule collect_calr:
    """Gather per-sample CALR tables into the two run-level result files."""
    input:
        type1=calr_type1_files,
        type2=calr_type2_files,
        samples=SAMPLES_FOR_PROCESSING,
    output:
        type1=CALR_TYPE1_RESULTS,
        type2=CALR_TYPE2_RESULTS,
    log:
        os.path.join(VCFS_DIR, "logs", "collect_calr.log"),
    script:
        "../scripts/collect_calr.py"
