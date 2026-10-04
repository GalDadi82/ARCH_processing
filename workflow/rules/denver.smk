# Stage 5: cross-reference the master mutation table against the Denver normals cohort.
#
# A mutation passes if it was never seen in Denver, or was seen and is on Denver's list of
# valid SNVs / indels (separate lists), or has a narrow listed match.


rule denver_query:
    """Ask the Denver variant index whether each master-table mutation was ever observed."""
    input:
        master=MASTER,
    output:
        is_in_denver=IS_IN_DENVER,
    params:
        query=os.path.join(config["denver"]["base_dir"], "variant_index", "query.py"),
    log:
        os.path.join(VCFS_DIR, "logs", "denver_query.log"),
    resources:
        mem_mb=16000,
    shell:
        r"""
        python {params.query:q} \
            --tsv {input.master:q} \
            -o {output.is_in_denver:q} > {log:q} 2>&1
        """


rule denver_filter:
    """Merge the Denver verdicts onto the master table and keep the valid mutations."""
    input:
        master=MASTER,
        is_in_denver=IS_IN_DENVER,
        filtered_snvs=from_repo(config["denver"]["filtered_snvs"]),
        filtered_indels=from_repo(config["denver"]["filtered_indels"]),
    output:
        with_denver=MASTER_WITH_DENVER,
        filtered=DENVER_FILTERED,
    log:
        os.path.join(VCFS_DIR, "logs", "denver_filter.log"),
    script:
        "../scripts/denver_filter.py"
