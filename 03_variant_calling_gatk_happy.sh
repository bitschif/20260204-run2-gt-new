#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking - GATK hap.py
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="gatk"
log_info "===== STEP 07: Benchmarking - ${CALLER} hap.py ====="

check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"
NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
check_file "${NORMALIZED_VCF}" || exit 1

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

log_info "Switching to happy-py27 environment for hap.py..."
if command -v conda &> /dev/null; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate happy-py27
else
    log_warn "conda not found; ensure happy-py27 is active before running hap.py"
fi

log_info "Benchmarking ${CALLER} with hap.py..."

BENCH_CALLER="${BENCH_DIR}/${CALLER}"
ensure_dir "${BENCH_CALLER}/happy"
HAPPY_PREFIX="${BENCH_CALLER}/happy/${PREFIX}_${CALLER}"
run_happy "${TRUTH_NORM}" "${NORMALIZED_VCF}" "${HAPPY_PREFIX}"

SUMMARY_CSV="${HAPPY_PREFIX}.summary.csv"
METRICS_TSV="${BENCH_CALLER}/${PREFIX}_${CALLER}_metrics.tsv"
if [[ -f "${SUMMARY_CSV}" ]]; then
    write_happy_metrics "${SUMMARY_CSV}" "${CALLER}" "${HAPPY_PREFIX}" "${METRICS_TSV}"
else
    log_warn "hap.py summary not found for ${CALLER}"
fi
