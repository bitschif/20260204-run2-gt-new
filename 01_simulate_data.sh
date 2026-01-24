#!/bin/bash
#===============================================================================
# STEP 01: Simulate Data
# Tools: simutator (iqbal-lab-org/simutator), art_illumina
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 01: Simulate Data ====="
start_timer

#-------------------------------------------------------------------------------
# 1. Check tools and reference
#-------------------------------------------------------------------------------
check_tool simutator || exit 1
check_tool art_illumina || exit 1
check_file "${REF_FASTA}" || exit 1

#-------------------------------------------------------------------------------
# 2. Generate mutations with simutator
#-------------------------------------------------------------------------------
log_info "Generating mutations with simutator..."
log_info "  SNPs: every ${SNP_DIST} bp"
log_info "  Deletions: ${DEL_LEN} bp every ${DEL_DIST} bp"
log_info "  Insertions: ${INS_LEN} bp every ${INS_DIST} bp"

SIM_PREFIX="${SIM_DIR}/${PREFIX}"

simutator mutate_fasta \
    --snps ${SNP_DIST} \
    --dels ${DEL_DIST}:${DEL_LEN} \
    --ins ${INS_DIST}:${INS_LEN} \
    --seed ${SEED} \
    "${REF_FASTA}" \
    "${SIM_PREFIX}"

check_exit "simutator"

#-------------------------------------------------------------------------------
# 3. Collect ALL simutator *.original.vcf (each mutation type is independent)
#    Then merge into a single truth VCF (SNP + DEL + INS)
#    And build a combined mutated FASTA using bcftools consensus
#-------------------------------------------------------------------------------

log_info "Collecting simutator outputs (NOTE: simutator writes one set per mutation type)..."

# Grab all truth VCFs where reference is the ORIGINAL fasta
shopt -s nullglob
ORIGINAL_VCFS=( "${SIM_PREFIX}"*.original.vcf )
shopt -u nullglob

if (( ${#ORIGINAL_VCFS[@]} == 0 )); then
    log_error "No simutator *.original.vcf found with prefix: ${SIM_PREFIX}"
    ls -la "${SIM_DIR}/"
    exit 1
fi

log_info "Found ${#ORIGINAL_VCFS[@]} simutator truth VCF(s):"
for v in "${ORIGINAL_VCFS[@]}"; do
    log_info "  - ${v}"
done

# Temporary working directory for merge steps
TMP_MERGE_DIR="$(mktemp -d -p "${SIM_DIR}" "${PREFIX}.merge_truth.XXXXXX")"
cleanup_tmp() { rm -rf "${TMP_MERGE_DIR}"; }
trap cleanup_tmp EXIT

# Helper: add missing GT FORMAT header only if absent
fix_gt_header_if_missing() {
    local in_vcf="$1"
    local out_vcf="$2"

    if grep -q '^##FORMAT=<ID=GT,' "${in_vcf}"; then
        cp -f "${in_vcf}" "${out_vcf}"
        return 0
    fi

    awk '
    /^##fileformat/ {
        print $0
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        next
    }
    { print $0 }
    ' "${in_vcf}" > "${out_vcf}"
}

# Fix headers (if needed), then combine bodies under one header, then sort
FIXED_VCFS=()
for vcf in "${ORIGINAL_VCFS[@]}"; do
    fixed="${TMP_MERGE_DIR}/$(basename "${vcf}").fixed.vcf"
    fix_gt_header_if_missing "${vcf}" "${fixed}"
    FIXED_VCFS+=( "${fixed}" )
done

# Build a single VCF stream: header from first + records from all, then sort+bgzip
RAW_MERGED_VCF="${TMP_MERGE_DIR}/${PREFIX}.truth.raw.sorted.vcf.gz"

{
    bcftools view -h "${FIXED_VCFS[0]}"
    for f in "${FIXED_VCFS[@]}"; do
        bcftools view -H "${f}"
    done
} | bcftools sort -Oz -o "${RAW_MERGED_VCF}" -

# Normalize vs reference, de-duplicate, then write final TRUTH_VCF
# --check-ref x: drop records whose REF doesn't match the reference (helps if rare overlaps happen)
bcftools norm \
    --check-ref x \
    -f "${REF_FASTA}" \
    -m -both \
    -d both \
    "${RAW_MERGED_VCF}" \
    -Ou \
| bcftools sort -Oz -o "${TRUTH_VCF}" -

tabix -p vcf "${TRUTH_VCF}"

log_info "Combined truth VCF written: ${TRUTH_VCF}"

# Create the combined mutated FASTA from REF_FASTA + TRUTH_VCF
# This FASTA will now contain SNP + DEL + INS simultaneously.
MUTATED_FASTA="${SIM_DIR}/${PREFIX}_mutated_combined.fa"

SAMPLE_NAME="$(bcftools query -l "${TRUTH_VCF}" | head -n 1)"
if [[ -z "${SAMPLE_NAME}" ]]; then
    log_error "No sample found in TRUTH_VCF (cannot run bcftools consensus)"
    exit 1
fi

bcftools consensus -s "${SAMPLE_NAME}" -f "${REF_FASTA}" "${TRUTH_VCF}" > "${MUTATED_FASTA}"
samtools faidx "${MUTATED_FASTA}"

log_info "Combined mutated FASTA written: ${MUTATED_FASTA}"

# Create separate SNP and INDEL files
TRUTH_TYPED_VCF="${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz"
bcftools +fill-tags "${TRUTH_VCF}" -Oz -o "${TRUTH_TYPED_VCF}" -- -t TYPE
tabix -p vcf "${TRUTH_TYPED_VCF}"

bcftools view -i 'TYPE="snp"' "${TRUTH_TYPED_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
bcftools view -i 'TYPE="indel"' "${TRUTH_TYPED_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

# Count variants
TOTAL_VARS=$(bcftools view -H "${TRUTH_VCF}" | wc -l)
SNP_COUNT=$(bcftools view -H "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz" | wc -l)
INDEL_COUNT=$(bcftools view -H "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz" | wc -l)

log_info "  Total variants: ${TOTAL_VARS}"
log_info "  SNPs: ${SNP_COUNT}"
log_info "  Indels: ${INDEL_COUNT}"

#-------------------------------------------------------------------------------
# 5. Generate reads with ART
#-------------------------------------------------------------------------------
log_info "Generating reads with ART Illumina..."
log_info "  Coverage: ${COVERAGE}x"
log_info "  Read length: ${READ_LENGTH} bp"
log_info "  Fragment: ${FRAGMENT_MEAN} +/- ${FRAGMENT_SD} bp"

READ_PREFIX="${SIM_DIR}/${PREFIX}"

art_illumina \
    -ss ${ART_PLATFORM} \
    -i "${MUTATED_FASTA}" \
    -p \
    -l ${READ_LENGTH} \
    -f ${COVERAGE} \
    -m ${FRAGMENT_MEAN} \
    -s ${FRAGMENT_SD} \
    -rs ${SEED} \
    -o "${READ_PREFIX}_" \
    -na

check_exit "art_illumina"

# Rename and compress output files
mv "${READ_PREFIX}_1.fq" "${READ_PREFIX}_R1.fastq"
mv "${READ_PREFIX}_2.fq" "${READ_PREFIX}_R2.fastq"
gzip -f "${READ_PREFIX}_R1.fastq"
gzip -f "${READ_PREFIX}_R2.fastq"

log_info "  R1: ${READ_PREFIX}_R1.fastq.gz"
log_info "  R2: ${READ_PREFIX}_R2.fastq.gz"

#-------------------------------------------------------------------------------
# 6. Create callable regions BED
#-------------------------------------------------------------------------------
log_info "Creating callable regions BED..."

awk -v OFS='\t' '{print $1, 0, $2}' "${REF_FAI}" > "${HIGH_CONF_BED}"

log_info "  BED: ${HIGH_CONF_BED}"

end_timer "01_simulate_data"
log_info "===== Simulation Complete ====="
