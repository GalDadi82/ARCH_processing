#!/bin/bash
# Process one sample: filter VCF (bcftools) + ANNOVAR.
# Usage: run_annovar_one_lsf.sh BASE_DIR SAMPLE_NAME REBUILD
#   REBUILD: 0 = skip if output exists, 1 = run anyway
set -e
BASE_DIR="${1:?Usage: $0 BASE_DIR SAMPLE_NAME REBUILD}"
SAMPLE_NAME="${2:?Usage: $0 BASE_DIR SAMPLE_NAME REBUILD}"
REBUILD="${3:-0}"
vcf_file="${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.varscan_C8_AF0.001_Q15_norm.dedup.vcf.gz"
filtered_vcf_file="${BASE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.varscan_C8_AF0.001_Q15_norm.dedup.filtered.vcf.gz"
output_anvr_varscan="${BASE_DIR}/${SAMPLE_NAME}/varscan_annot_dedup_filt"

if [[ -f "${output_anvr_varscan}.hg38_multianno.txt" && "$REBUILD" != "1" ]]; then
  echo "${SAMPLE_NAME}: Skipping ANNOVAR, output already exists."
  exit 0
fi

echo "${SAMPLE_NAME}: Filtering VAF for actual mutations..."
module load BCFtools
bcftools view -i 'ALT != "."' "$vcf_file" -Oz -o "$filtered_vcf_file"

# perl="/apps/RH7U2/gnu/perl/5.24.0/bin/perl"
# annovar="/apps/RH7U2/general/annovar/2017Jun01/table_annovar.pl"
# human_db="/home/projects/shlush/yardenab/humandb/"
# db_list="refGene,gnomad211_genome,kaviar_20150923,cosmic70,avsnp150,clinvar_20170905"

perl="/usr/bin/perl"
annovar="/home/projects/shlush/shared/annovar/2022-08-02/table_annovar.pl"
human_db="/home/projects/shlush/shared/annovar/humandb/"
db_list="refGene,gnomad41_genome,kaviar_20150923,cosmic103_mapped,avsnp151,clinvar_20250721,dbnsfp47a"

hg="hg38"

operation="g,f,f,f,f,f,f"

echo "${SAMPLE_NAME}: Running ANNOVAR..."
$perl $annovar $filtered_vcf_file $human_db -buildver $hg -out $output_anvr_varscan -remove -protocol $db_list -operation $operation -nastring . -vcfinput

echo "${SAMPLE_NAME}: Done."
