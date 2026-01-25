#!/bin/bash
#===============================================================================
# Normalize VCF for benchmarking
# Usage: normalize_vcf.sh <input_vcf> <output_vcf> <reference_fasta>
#===============================================================================

set -euo pipefail

input_vcf="${1:-}"
output_vcf="${2:-}"
reference_fasta="${3:-}"

if [[ -z "${input_vcf}" || -z "${output_vcf}" || -z "${reference_fasta}" ]]; then
    echo "[ERROR] Usage: normalize_vcf.sh <input_vcf> <output_vcf> <reference_fasta>" >&2
    exit 1
fi

if ! command -v bcftools &>/dev/null; then
    echo "[ERROR] Tool not found: bcftools" >&2
    exit 1
fi

if ! command -v tabix &>/dev/null; then
    echo "[ERROR] Tool not found: tabix" >&2
    exit 1
fi

tmp_vcf="${output_vcf}.tmp"
bcftools norm -f "${reference_fasta}" -m -both "${input_vcf}" -Oz -o "${tmp_vcf}"
bcftools sort "${tmp_vcf}" -Oz -o "${output_vcf}"
rm -f "${tmp_vcf}"
tabix -f -p vcf "${output_vcf}"
