# Stage 1: link per-sample job output into All_Samples, then build the QC report and
# decide which samples get processed.


rule link_sample_dirs:
    """Populate All_Samples/{2_mapped,5_vcfs} with symlinks to each job's sample folders."""
    output:
        manifest=LINK_MANIFEST,
    params:
        analysis_dir=ANALYSIS_DIR,
        output_base_dir=OUTPUT_BASE_DIR,
    log:
        os.path.join(VCFS_DIR, "logs", "link_sample_dirs.log"),
    script:
        "../scripts/link_sample_dirs.py"


rule qc_statistics:
    """Run utils/mip_run_statistics.py to produce the run QC report and sample classification."""
    input:
        manifest=LINK_MANIFEST,
    output:
        classified=ALL_CLASSIFIED,
    params:
        script=os.path.join(UTILS_DIR, "mip_run_statistics.py"),
        output_dir=VCFS_DIR,
        seq_run=SEQ_RUN,
        panel=PANEL,
        runs_base_dir=RUNS_BASE_DIR,
        sample_fixed_column=config["qc"]["sample_fixed_column"],
        min_lib_size=config["qc"]["min_lib_size"],
        min_base_cov=config["qc"]["min_base_cov"],
        min_f_dedup_cov=config["qc"]["min_f_dedup_cov"],
        rebuild="--rebuild" if config["qc"]["rebuild"] else "",
    log:
        os.path.join(VCFS_DIR, "logs", "qc_statistics.log"),
    shell:
        r"""
        python {params.script} \
            --output-dir {params.output_dir:q} \
            --seq-run-name {params.seq_run:q} \
            --panel {params.panel:q} \
            --runs-base-dir {params.runs_base_dir:q} \
            --sample-fixed-column {params.sample_fixed_column:q} \
            --min-lib-size {params.min_lib_size} \
            --min-base-cov {params.min_base_cov} \
            --min-f-dedup-cov {params.min_f_dedup_cov} \
            {params.rebuild} > {log:q} 2>&1
        """


checkpoint select_samples:
    """Pick the samples to process. This is the checkpoint the per-sample fan-out hangs off."""
    input:
        classified=ALL_CLASSIFIED,
    output:
        SAMPLES_FOR_PROCESSING,
    params:
        min_base_cov=config["qc"]["min_base_cov"],
        min_f_dedup_cov_to_process=config["qc"]["min_f_dedup_cov_to_process"],
    log:
        os.path.join(VCFS_DIR, "logs", "select_samples.log"),
    script:
        "../scripts/select_samples.py"
