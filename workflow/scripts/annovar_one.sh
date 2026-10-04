#!/bin/bash
# Process one sample: drop non-variant rows from its varscan VCF (bcftools), then annotate
# with ANNOVAR.
#
# Parameterised version of notebooks/run_annovar_one_lsf.sh -- the settings come from
# workflow/config/config.yaml instead of being hard-coded, and there is no REBUILD check
# because Snakemake decides whether the rule needs to run.
set -euo pipefail

bcftools_module=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-vcf)          in_vcf="$2"; shift 2 ;;
    --filtered-vcf)    filtered_vcf="$2"; shift 2 ;;
    --out-prefix)      out_prefix="$2"; shift 2 ;;
    --perl)            perl_exe="$2"; shift 2 ;;
    --table-annovar)   table_annovar="$2"; shift 2 ;;
    --humandb)         humandb="$2"; shift 2 ;;
    --buildver)        buildver="$2"; shift 2 ;;
    --protocol)        protocol="$2"; shift 2 ;;
    --operation)       operation="$2"; shift 2 ;;
    --bcftools-module) bcftools_module="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

: "${in_vcf:?--in-vcf is required}"
: "${filtered_vcf:?--filtered-vcf is required}"
: "${out_prefix:?--out-prefix is required}"

if [[ -n "$bcftools_module" ]]; then
  # `module` is a shell function, so it is not defined in a non-interactive shell.
  if ! declare -F module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
      # shellcheck disable=SC1090
      [[ -f "$init" ]] && source "$init" && break
    done
  fi
  if declare -F module >/dev/null 2>&1; then
    module load "$bcftools_module"
  else
    echo "WARNING: environment modules unavailable; relying on bcftools from PATH" >&2
  fi
fi

echo "Filtering VAF for actual mutations..."
bcftools view -i 'ALT != "."' "$in_vcf" -Oz -o "$filtered_vcf"

echo "Running ANNOVAR..."
"$perl_exe" "$table_annovar" \
  "$filtered_vcf" \
  "$humandb" \
  -buildver "$buildver" \
  -out "$out_prefix" \
  -remove \
  -protocol "$protocol" \
  -operation "$operation" \
  -nastring . \
  -vcfinput

echo "Done."
