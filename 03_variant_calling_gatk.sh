#!/bin/bash
#===============================================================================
# STEP 03: Variant Calling - GATK HaplotypeCaller
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="gatk"
log_info "===== STEP 03: ${CALLER} HaplotypeCaller ====="
start_timer

# Input
source "${PREPROC_DIR}/bam_path.sh"
check_file "${FINAL_BAM}" || exit 1
check_tool bcftools || exit 1
check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"

#-------------------------------------------------------------------------------
# 1. Run HaplotypeCaller
#-------------------------------------------------------------------------------
log_info "Running GATK HaplotypeCaller..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"

gatk HaplotypeCaller \
    --java-options "${JAVA_OPTS}" \
    -R "${REF_FASTA}" \
    -I "${FINAL_BAM}" \
    -O "${RAW_VCF}" \
    --standard-min-confidence-threshold-for-calling "${GATK_STAND_CALL_CONF}" \
    --native-pair-hmm-threads "${THREADS}" \
    2>&1 | tee "${LOG_DIR}/${CALLER}.log"

check_exit "HaplotypeCaller"

#-------------------------------------------------------------------------------
# 2. Hard filtering (GATK Best Practices)
#-------------------------------------------------------------------------------
log_info "Applying hard filters..."

# SNPs
gatk SelectVariants -V "${RAW_VCF}" -select-type SNP -O "${OUT_DIR}/snps_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/snps_raw.vcf.gz" \
    -O "${OUT_DIR}/snps_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
    --filter-name "SNP_FILTER"

# INDELs
gatk SelectVariants -V "${RAW_VCF}" -select-type INDEL -O "${OUT_DIR}/indels_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/indels_raw.vcf.gz" \
    -O "${OUT_DIR}/indels_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
    --filter-name "INDEL_FILTER"

# Merge
gatk MergeVcfs \
    -I "${OUT_DIR}/snps_filtered.vcf.gz" \
    -I "${OUT_DIR}/indels_filtered.vcf.gz" \
    -O "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz"

#-------------------------------------------------------------------------------
# 3. Extract PASS variants
#-------------------------------------------------------------------------------
log_info "Extracting PASS variants..."

PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"
bcftools view -f PASS "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz" -Oz -o "${PASS_VCF}"
tabix -p vcf "${PASS_VCF}"

# Split by type for benchmarking
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Normalize for benchmarking
#-------------------------------------------------------------------------------
log_info "Normalizing variants for benchmarking..."

NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
"${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

#-------------------------------------------------------------------------------
# 5. Benchmarking with hap.py
#-------------------------------------------------------------------------------
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
    update_benchmark_summary "${CALLER}" "${METRICS_TSV}"
else
    log_warn "hap.py summary not found for ${CALLER}"
fi

#-------------------------------------------------------------------------------
# 6. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)

log_info "Results:  $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

end_timer "03_${CALLER}"
log_info "===== ${CALLER} Complete ====="
