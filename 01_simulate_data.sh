#!/bin/bash
#===============================================================================
# STEP 01: Simulate Data
# Tools: simutator (iqbal-lab-org/simutator), art_illumina
# Simutator tạo mutations, ART tạo reads
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config. sh"
source "${SCRIPT_DIR}/scripts/helper_functions. sh"

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
# Syntax: simutator mutate_fasta [options] in. fasta out_prefix
# Options: --snps DIST, --dels DIST: LEN, --ins DIST: LEN, --seed
#-------------------------------------------------------------------------------
log_info "Generating mutations with simutator..."
log_info "  SNPs: every ${SNP_DIST} bp"
log_info "  Deletions: ${DEL_LEN} bp every ${DEL_DIST} bp"
log_info "  Insertions: ${INS_LEN} bp every ${INS_DIST} bp"

SIM_PREFIX="${SIM_DIR}/${PREFIX}"

# Run simutator với tất cả mutation types trong 1 lệnh
simutator mutate_fasta \
    --snps ${SNP_DIST} \
    --dels ${DEL_DIST}:${DEL_LEN} \
    --ins ${INS_DIST}:${INS_LEN} \
    --seed ${SEED} \
    "${REF_FASTA}" \
    "${SIM_PREFIX}"

check_exit "simutator"

#-------------------------------------------------------------------------------
# 3. Find output files from simutator
# Output format: prefix. snp.dist-X.del. dist-Y-len-Z.ins.dist-A-len-B.fa
#-------------------------------------------------------------------------------
log_info "Processing simutator output..."

# Find mutated FASTA (có thể có nhiều patterns)
MUTATED_FASTA=$(ls ${SIM_PREFIX}*.fa 2>/dev/null | grep -v ".original" | head -1)
ORIGINAL_VCF=$(ls ${SIM_PREFIX}*.original.vcf 2>/dev/null | head -1)

if [[ -z "${MUTATED_FASTA}" ]] || [[ -z "${ORIGINAL_VCF}" ]]; then
    log_error "simutator output files not found"
    ls -la "${SIM_DIR}/"
    exit 1
fi

log_info "  Mutated FASTA: ${MUTATED_FASTA}"
log_info "  Truth VCF: ${ORIGINAL_VCF}"

# Index mutated FASTA
samtools faidx "${MUTATED_FASTA}"

#-------------------------------------------------------------------------------
# 4. Process truth VCF
#-------------------------------------------------------------------------------
log_info "Processing truth VCF..."

# Sort, compress, index
bcftools sort "${ORIGINAL_VCF}" -Oz -o "${TRUTH_VCF}"
tabix -p vcf "${TRUTH_VCF}"

# Create separate SNP and INDEL files for benchmarking
bcftools view -v snps "${TRUTH_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
bcftools view -v indels "${TRUTH_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

# Count variants
N_SNP=$(bcftools view -H -v snps "${TRUTH_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${TRUTH_VCF}" | wc -l)
N_TOTAL=$((N_SNP + N_INDEL))

log_info "Truth variants: ${N_TOTAL} total (${N_SNP} SNPs, ${N_INDEL} INDELs)"

#-------------------------------------------------------------------------------
# 5. Simulate reads with ART Illumina
#-------------------------------------------------------------------------------
log_info "Simulating Illumina reads with ART..."
log_info "  Coverage: ${COVERAGE}x"
log_info "  Read length: ${READ_LENGTH} bp"
log_info "  Fragment:  ${FRAGMENT_MEAN} ± ${FRAGMENT_SD} bp"
log_info "  Platform: ${ART_PLATFORM}"

READ_PREFIX="${SIM_DIR}/${PREFIX}"

art_illumina \
    -ss "${ART_PLATFORM}" \
    -i "${MUTATED_FASTA}" \
    -p \
    -l "${READ_LENGTH}" \
    -f "${COVERAGE}" \
    -m "${FRAGMENT_MEAN}" \
    -s "${FRAGMENT_SD}" \
    -o "${READ_PREFIX}_" \
    -na \
    --rndSeed "${SEED}" \
    2>&1 | tee "${LOG_DIR}/art_illumina. log"

check_exit "ART Illumina"

# Rename and compress
mv "${READ_PREFIX}_1.fq" "${READ_PREFIX}_R1.fastq"
mv "${READ_PREFIX}_2.fq" "${READ_PREFIX}_R2.fastq"

log_info "Compressing FASTQ files..."
gzip -f "${READ_PREFIX}_R1.fastq"
gzip -f "${READ_PREFIX}_R2.fastq"

#-------------------------------------------------------------------------------
# 6. Create callable regions BED
#-------------------------------------------------------------------------------
log_info "Creating callable regions BED..."

# Get chromosome length and create full region BED
awk -v OFS='\t' '{print $1, 0, $2}' "${REF_FAI}" > "${HIGH_CONF_BED}"

#-------------------------------------------------------------------------------
# 7. Summary
#-------------------------------------------------------------------------------
R1_COUNT=$(zcat "${READ_PREFIX}_R1.fastq.gz" | wc -l)
R1_READS=$((R1_COUNT / 4))

cat > "${SIM_DIR}/simulation_info.txt" << EOF
=== SIMULATION INFO ===
Reference: ${REF_FASTA}
Chromosome: ${CHR_TO_USE}
Mutated FASTA: ${MUTATED_FASTA}
Truth VCF: ${TRUTH_VCF}

Mutations:
  Total:  ${N_TOTAL}
  SNPs: ${N_SNP}
  INDELs: ${N_INDEL}

Simutator parameters:
  SNP distance: ${SNP_DIST} bp
  Deletion:  ${DEL_LEN} bp every ${DEL_DIST} bp
  Insertion: ${INS_LEN} bp every ${INS_DIST} bp

Read simulation:
  Read pairs: ${R1_READS}
  Coverage: ${COVERAGE}x
  Read length: ${READ_LENGTH} bp
  Fragment size: ${FRAGMENT_MEAN} ± ${FRAGMENT_SD} bp
  Platform: ${ART_PLATFORM}

Output files:
  R1: ${READ_PREFIX}_R1.fastq.gz
  R2: ${READ_PREFIX}_R2.fastq.gz
EOF

cat "${SIM_DIR}/simulation_info.txt"

end_timer "01_simulate_data"
log_info "===== Simulation Complete ====="
