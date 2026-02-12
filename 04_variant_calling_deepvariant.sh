#!/bin/bash
#===============================================================================
# STEP 04: Variant Calling - DeepVariant (via Docker)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="deepvariant"
log_info "===== STEP 04: ${CALLER} (Docker) ====="
start_timer

# Check Docker
check_tool docker || exit 1

# Input
source "${PREPROC_DIR}/bam_path.sh"
check_file "${FINAL_BAM}" || exit 1
check_tool bcftools || exit 1
check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"
ensure_dir "${OUT_DIR}/intermediate"

#-------------------------------------------------------------------------------
# 1. Prepare paths for Docker
#-------------------------------------------------------------------------------
# Docker needs absolute paths
ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")
WES_BED_BASENAME=$(basename "${WES_BED}")

#-------------------------------------------------------------------------------
# 2. Run DeepVariant via Docker
#-------------------------------------------------------------------------------
log_info "Running DeepVariant via Docker..."
log_info "  Image: ${DEEPVARIANT_IMAGE}"
log_info "  Model: WES"

DOCKER_USER_NAME="${RUN_AS_USER:-${USER:-hxt}}"
DOCKER_USER="$(id -u "${DOCKER_USER_NAME}")":"$(id -g "${DOCKER_USER_NAME}")"

docker run \
    --rm \
    --user "${DOCKER_USER}" \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    ${DEEPVARIANT_IMAGE} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref="/ref/${REF_BASENAME}" \
    --reads="/input/${BAM_BASENAME}" \
    --output_vcf="/output/${PREFIX}_${CALLER}_raw.vcf.gz" \
    --output_gvcf="/output/${PREFIX}_${CALLER}.g.vcf.gz" \
    --intermediate_results_dir="/output/intermediate" \
    --num_shards="${THREADS}" \
    --regions="/ref/${WES_BED_BASENAME}" \
    2>&1 | tee "${LOG_DIR}/${CALLER}.log"

check_exit "DeepVariant"

#-------------------------------------------------------------------------------
# 3. Process output
#-------------------------------------------------------------------------------
log_info "Processing DeepVariant output..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"
PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"

# Sanity check output
check_file "${RAW_VCF}" || exit 1
if [[ ! -s "${RAW_VCF}" ]]; then
    log_error "DeepVariant output is empty: ${RAW_VCF}"
    exit 1
fi

# Helper: attempt to index, with a bgzip fallback if the file is plain gzip.
index_vcf() {
    local vcf_path="$1"
    if tabix -f -p vcf "${vcf_path}"; then
        return 0
    fi
    log_warn "tabix indexing failed for ${vcf_path}; attempting bgzip re-compression and retry."
    if check_tool bgzip; then
        local tmp_bgz="${vcf_path}.bgz"
        gunzip -c "${vcf_path}" | bgzip -c > "${tmp_bgz}"
        mv -f "${tmp_bgz}" "${vcf_path}"
        tabix -f -p vcf "${vcf_path}"
        return $?
    fi
    log_error "bgzip not available; cannot re-compress ${vcf_path} for tabix indexing."
    return 1
}

# Ensure sorted/bgzip for tabix (DeepVariant output can be unsorted)
SORT_TMP="${OUT_DIR}/intermediate/bcftools_sort"
ensure_dir "${SORT_TMP}"
SORTED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.sorted.vcf.gz"
bcftools sort "${RAW_VCF}" -Oz -o "${SORTED_VCF}" -T "${SORT_TMP}"
mv -f "${SORTED_VCF}" "${RAW_VCF}"

# Index
index_vcf "${RAW_VCF}"

# DeepVariant outputs mostly PASS, but filter to be safe
bcftools view -f "PASS,." "${RAW_VCF}" -Oz -o "${PASS_VCF}"
index_vcf "${PASS_VCF}"

# Split by type
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
index_vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
index_vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

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
# 5. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)

log_info "Results: $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

# Cleanup intermediate files
rm -rf "${OUT_DIR}/intermediate"

end_timer "04_${CALLER}"
log_info "===== ${CALLER} Complete ====="
