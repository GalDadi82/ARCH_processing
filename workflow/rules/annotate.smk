# Stage 2: per-sample ANNOVAR annotation, per-sample reformatting, then the merge into
# a single mutation table.
#
# The notebook had two interchangeable cells here -- an in-process loop and a
# "submit one bsub per sample" cell. Snakemake owns the scheduling now, so there is one
# rule and the LSF profile decides whether it runs locally or as a cluster job.


rule annovar:
    """Drop non-variant rows from a sample's varscan VCF, then annotate it with ANNOVAR."""
    input:
        vcf=sample_vcf,
    output:
        txt=MULTIANNO_TXT,
        vcf=MULTIANNO_VCF,
    params:
        script=os.path.join(workflow.basedir, "scripts", "annovar_one.sh"),
        out_prefix=os.path.join(VCFS_DIR, "{sample}", config["annovar"]["out_basename"]),
        filtered_vcf=lambda w: sample_vcf(w).replace(".vcf.gz", ".filtered.vcf.gz"),
        perl=config["annovar"]["perl"],
        table_annovar=config["annovar"]["table_annovar"],
        humandb=config["annovar"]["humandb"],
        buildver=config["annovar"]["buildver"],
        protocol=config["annovar"]["protocol"],
        operation=config["annovar"]["operation"],
        bcftools_module=config["annovar"]["bcftools_module"],
    log:
        os.path.join(VCFS_DIR, "logs", "annovar", "{sample}.log"),
    resources:
        mem_mb=8000,
        runtime=240,
    shell:
        r"""
        bash {params.script:q} \
            --in-vcf {input.vcf:q} \
            --filtered-vcf {params.filtered_vcf:q} \
            --out-prefix {params.out_prefix:q} \
            --perl {params.perl:q} \
            --table-annovar {params.table_annovar:q} \
            --humandb {params.humandb:q} \
            --buildver {params.buildver:q} \
            --protocol {params.protocol:q} \
            --operation {params.operation:q} \
            --bcftools-module {params.bcftools_module:q} > {log:q} 2>&1
        """


rule format_multianno:
    """Flatten one ANNOVAR multianno table: keep the annotations plus FILTER/Depth/VAF."""
    input:
        txt=MULTIANNO_TXT,
        vcf=MULTIANNO_VCF,
    output:
        formatted=MULTIANNO_FORMATTED,
    log:
        os.path.join(VCFS_DIR, "logs", "format_multianno", "{sample}.log"),
    script:
        "../scripts/format_multianno.py"


rule multianno_report:
    """Per-sample multianno file sizes -- the notebook's QC spot-check for empty outputs."""
    input:
        txt=multianno_raw,
        samples=SAMPLES_FOR_PROCESSING,
    output:
        report=MULTIANNO_REPORT,
    log:
        os.path.join(VCFS_DIR, "logs", "multianno_report.log"),
    script:
        "../scripts/multianno_report.py"


rule merge_formatted:
    """Concatenate every sample's formatted table into one merged varscan annotation file."""
    input:
        formatted=multianno_formatted,
        samples=SAMPLES_FOR_PROCESSING,
    output:
        merged=MERGED_VARSCAN,
    log:
        os.path.join(VCFS_DIR, "logs", "merge_formatted.log"),
    resources:
        mem_mb=16000,
    script:
        "../scripts/merge_formatted.py"
