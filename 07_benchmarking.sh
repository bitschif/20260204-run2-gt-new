#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking
# Compare variant calls against truth set
# Tools: bcftools isec (fallback), hap.py, RTG Tools
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Benchmarking ====="
start_timer

check_tool "bcftools" || exit 1

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")
TRUTH_SNP="${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
TRUTH_INDEL="${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

check_file "${TRUTH_VCF}" || exit 1

#-------------------------------------------------------------------------------
# Function:  Benchmark with bcftools isec
#-------------------------------------------------------------------------------
benchmark_bcftools() {
    local caller=$1
    local query_vcf=$2
    local truth_vcf=$3
    local out_dir=$4
    local vtype=$5
    
    ensure_dir "${out_dir}"
    
    # Intersection
    bcftools isec -p "${out_dir}" -Oz "${truth_vcf}" "${query_vcf}" 2>/dev/null || true
    
    # Count:  0000=FN, 0001=FP, 0002=TP(truth), 0003=TP(query)
    local fn=$(bcftools view -H "${out_dir}/0000.vcf.gz" 2>/dev/null | wc -l || echo 0)
    local fp=$(bcftools view -H "${out_dir}/0001.vcf.gz" 2>/dev/null | wc -l || echo 0)
    local tp=$(bcftools view -H "${out_dir}/0002.vcf.gz" 2>/dev/null | wc -l || echo 0)
    
    # Calculate metrics
    local precision=0
    local recall=0
    local f1=0

    precision=$(awk -v tp="${tp}" -v fp="${fp}" 'BEGIN {
        denom = tp + fp;
        if (denom > 0) { printf "%.6f", tp / denom; } else { printf "0"; }
    }')

    recall=$(awk -v tp="${tp}" -v fn="${fn}" 'BEGIN {
        denom = tp + fn;
        if (denom > 0) { printf "%.6f", tp / denom; } else { printf "0"; }
    }')

    f1=$(awk -v p="${precision}" -v r="${recall}" 'BEGIN {
        denom = p + r;
        if (denom > 0) { printf "%.6f", 2 * p * r / denom; } else { printf "0"; }
    }')
    
    # Write summary
    cat > "${out_dir}/summary.txt" << EOF
Caller: ${caller}
VariantType: ${vtype}
TP: ${tp}
FP: ${fp}
FN:  ${fn}
Precision: ${precision}
Recall: ${recall}
F1: ${f1}
EOF
    
    echo -e "${caller}\t${vtype}\t${tp}\t${fp}\t${fn}\t${precision}\t${recall}\t${f1}"
}

#-------------------------------------------------------------------------------
# Main benchmarking
#-------------------------------------------------------------------------------
SUMMARY="${BENCH_DIR}/benchmark_summary.tsv"
echo -e "Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1" > "${SUMMARY}"

for caller in "${CALLERS[@]}"; do
    log_info "Benchmarking ${caller}..."
    
    QUERY_VCF="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.vcf.gz"
    QUERY_SNP="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_snp.vcf.gz"
    QUERY_INDEL="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_indel.vcf.gz"
    
    if [[ !  -f "${QUERY_VCF}" ]]; then
        log_warn "  VCF not found, skipping..."
        continue
    fi
    
    BENCH_CALLER="${BENCH_DIR}/${caller}"
    ensure_dir "${BENCH_CALLER}"
    
    # Benchmark ALL variants
    log_info "  ALL variants..."
    result=$(benchmark_bcftools "${caller}" "${QUERY_VCF}" "${TRUTH_VCF}" \
        "${BENCH_CALLER}/all" "ALL")
    echo "${result}" >> "${SUMMARY}"
    
    # Benchmark SNPs
    if [[ -f "${QUERY_SNP}" ]] && [[ -f "${TRUTH_SNP}" ]]; then
        log_info "  SNPs..."
        result=$(benchmark_bcftools "${caller}" "${QUERY_SNP}" "${TRUTH_SNP}" \
            "${BENCH_CALLER}/snp" "SNP")
        echo "${result}" >> "${SUMMARY}"
    fi
    
    # Benchmark INDELs
    if [[ -f "${QUERY_INDEL}" ]] && [[ -f "${TRUTH_INDEL}" ]]; then
        log_info "  INDELs..."
        result=$(benchmark_bcftools "${caller}" "${QUERY_INDEL}" "${TRUTH_INDEL}" \
            "${BENCH_CALLER}/indel" "INDEL")
        echo "${result}" >> "${SUMMARY}"
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
