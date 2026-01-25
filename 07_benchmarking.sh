#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking
# Compare variant calls against truth set
# Tools: hap.py (Illumina)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Benchmarking ====="
start_timer

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")

#-------------------------------------------------------------------------------
# Run hap.py per caller
#-------------------------------------------------------------------------------
for caller in "${CALLERS[@]}"; do
    case "${caller}" in
        gatk) HAPPY_SCRIPT="${SCRIPT_DIR}/03_variant_calling_gatk_happy.sh" ;;
        deepvariant) HAPPY_SCRIPT="${SCRIPT_DIR}/04_variant_calling_deepvariant_happy.sh" ;;
        strelka2) HAPPY_SCRIPT="${SCRIPT_DIR}/05_variant_calling_strelka2_happy.sh" ;;
        freebayes) HAPPY_SCRIPT="${SCRIPT_DIR}/06_variant_calling_freebayes_happy.sh" ;;
        *) HAPPY_SCRIPT="" ;;
    esac

    if [[ -n "${HAPPY_SCRIPT}" && -f "${HAPPY_SCRIPT}" ]]; then
        log_info "Running hap.py for ${caller}..."
        bash "${HAPPY_SCRIPT}"
    else
        log_warn "hap.py script not found for ${caller}"
    fi
done

#-------------------------------------------------------------------------------
# Main benchmarking summary
#-------------------------------------------------------------------------------
SUMMARY="${BENCH_DIR}/benchmark_summary.tsv"
echo -e "Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1\tROC_AUC\tRuntimeSeconds" > "${SUMMARY}"

declare -A RUNTIME_MAP=()
if [[ -f "${LOG_DIR}/runtime.csv" ]]; then
    while IFS=',' read -r step duration; do
        case "${step}" in
            03_gatk) RUNTIME_MAP[gatk]="${duration}" ;;
            04_deepvariant) RUNTIME_MAP[deepvariant]="${duration}" ;;
            05_strelka2) RUNTIME_MAP[strelka2]="${duration}" ;;
            06_freebayes) RUNTIME_MAP[freebayes]="${duration}" ;;
        esac
    done < "${LOG_DIR}/runtime.csv"
else
    log_warn "runtime.csv not found; runtime metrics will be NA."
fi

for caller in "${CALLERS[@]}"; do
    METRICS_TSV="${BENCH_DIR}/${caller}/${PREFIX}_${caller}_metrics.tsv"
    if [[ -f "${METRICS_TSV}" ]]; then
        update_benchmark_summary "${caller}" "${METRICS_TSV}" "${RUNTIME_MAP[${caller}]:-NA}"
    else
        log_warn "  Metrics not found for ${caller}. Run variant calling script first."
    fi
done

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
log_info "===== Benchmarking Summary ====="
if command -v column &>/dev/null; then
    column -t -s$'\t' "${SUMMARY}"
else
    cat "${SUMMARY}"
fi

end_timer "07_benchmarking"
log_info "===== Benchmarking Complete ====="
