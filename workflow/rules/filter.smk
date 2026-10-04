# Stage 3: the mutation filtering chain.
#   Blood 2023 filtering -> recurrence filtering -> VAF filter + manual exclusions.


rule blood2023_filter:
    """Apply the Blood-2023-paper filtering (utils/filter_mutations.py)."""
    input:
        merged=MERGED_VARSCAN,
    output:
        listed=LISTED,
        not_listed=NOT_LISTED,
    params:
        utils_dir=UTILS_DIR,
        output_dir=VCFS_DIR,
        run_name=FILTER_RUN_NAME,
        sample_column=config["filter_mutations"]["sample_column"],
        dups_input_file=config["filter_mutations"]["dups_input_file"],
    log:
        os.path.join(VCFS_DIR, "logs", "blood2023_filter.log"),
    resources:
        mem_mb=16000,
    script:
        "../scripts/blood2023_filter.py"


rule filter_recurrent:
    """Drop recurrent mutations lacking haematopoietic COSMIC / ClinVar support.

    Skipped (input copied through) when filter_recurrent_without_support is false.
    """
    input:
        listed=LISTED,
    output:
        filtered=NO_RECURRENT,
    params:
        enabled=config["filter_recurrent_without_support"],
        force_include=config["force_include_variants"],
    log:
        os.path.join(VCFS_DIR, "logs", "filter_recurrent.log"),
    script:
        "../scripts/filter_recurrent.py"


rule vaf_filter:
    """Apply the final VAF threshold and drop known artefact / germline calls."""
    input:
        mutations=NO_RECURRENT,
    output:
        master=MASTER,
    params:
        min_final_vaf=config["min_final_vaf"],
        gene_min_vaf_overrides=config["gene_min_vaf_overrides"],
        exclusion_keys=config["exclusion_keys"],
        excluded_mutations=config["excluded_mutations"],
    log:
        os.path.join(VCFS_DIR, "logs", "vaf_filter.log"),
    script:
        "../scripts/vaf_filter.py"


rule qc_plots:
    """Mutations-per-sample, recurrence-vs-VAF and per-gene plots from the notebook."""
    input:
        listed=LISTED,
        not_listed=NOT_LISTED,
        samples=SAMPLES_FOR_PROCESSING,
    output:
        pdf=QC_PLOTS,
    log:
        os.path.join(VCFS_DIR, "logs", "qc_plots.log"),
    script:
        "../scripts/qc_plots.py"
