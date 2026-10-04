# Stage 6: email alert on target mutations and high-VAF CALR calls.
#
# This rule sends real email. It is part of the default run (matching the notebook), but
# `snakemake mutations` stops before it, and alert.enabled: false removes it from `all`.


rule mutation_alert:
    """Send the alert email for target mutations and CALR hits above their VAF cutoffs."""
    input:
        mutations=FINAL_MUTATIONS,
        calr_type1=CALR_TYPE1_RESULTS,
        calr_type2=CALR_TYPE2_RESULTS,
    output:
        sentinel=touch(os.path.join(VCFS_DIR, f"mutation_alert_{SEQ_RUN}.sent")),
    params:
        script=os.path.join(UTILS_DIR, "mutation_alert.py"),
        alert_config=from_repo(config["alert"]["config"]),
        seq_run=SEQ_RUN,
    log:
        os.path.join(VCFS_DIR, "logs", "mutation_alert.log"),
    shell:
        r"""
        python {params.script:q} \
            --mutations-file {input.mutations:q} \
            --config {params.alert_config:q} \
            --run-name {params.seq_run:q} \
            --calr-type1-file {input.calr_type1:q} \
            --calr-type2-file {input.calr_type2:q} > {log:q} 2>&1
        """
