#!/bin/bash
# STEP 00: Setup Environment

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 00: Setup Environment ====="

# 1. Create directories
log_info "Creating directories..."
ALL_DIRS=(
    "${DATA_DIR}" "${REF_DIR}" "${SIM_DIR}" "${RESULTS_DIR}"
    "${LOG_DIR}" "${PREPROC_DIR}" "${VARIANT_DIR}" "${BENCH_DIR}"
    "${FIGURE_DIR}" "${METRICS_DIR}"
)
VARIANT_CALLERS=(gatk deepvariant strelka2 freebayes)

for dir in "${ALL_DIRS[@]}"; do ensure_dir "$dir"; done
for caller in "${VARIANT_CALLERS[@]}"; do
    ensure_dir "${VARIANT_DIR}/${caller}"
    ensure_dir "${BENCH_DIR}/${caller}"
done

echo "step,duration_seconds" > "${LOG_DIR}/runtime.csv"

# 2. Check required tools
log_info "Checking required tools..."
REQUIRED_TOOLS=(bwa samtools bcftools gatk fastp fastqc bgzip tabix art_illumina freebayes simutator)

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "$tool" &> /dev/null; then
        log_info "  ✓ $tool"
    else
        log_warn "  ✗ $tool (missing)"
    fi
done

if command -v docker &> /dev/null; then
    log_info "  ✓ docker (for DeepVariant, Strelka2)"
else
    log_warn "  ✗ docker (DeepVariant, Strelka2 will not run)"
fi

# 3. Download reference genome
log_info "Preparing reference genome..."
if [[ ! -f "${REF_FASTA}" ]]; then
    log_info "Downloading ${CHR_TO_USE}..."
    cd "${REF_DIR}"

    UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHR_TO_USE}.fa.gz"
    BACKUP_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${CHR_TO_USE}.fa.gz"

    if ! wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa.gz" "${UCSC_URL}" \
        && ! wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa.gz" "${BACKUP_URL}"; then
        log_error "Failed to download reference genome! Please download manually and place the file in ${REF_DIR}."
        exit 1
    fi

    log_info "Decompressing reference genome..."
    gunzip -f "${CHR_TO_USE}.fa.gz"

    if [[ ! -f "${REF_FASTA}" ]]; then
        log_error "Decompression failed!"
        exit 1
    fi

    log_info "Indexing reference genome..."
    samtools faidx "${REF_FASTA}"
    bwa index "${REF_FASTA}"
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}" 2>/dev/null || \
    samtools dict "${REF_FASTA}" > "${REF_DICT}"
else
    log_info "Reference already exists: ${REF_FASTA}"
fi

# 4. Verify reference files
log_info "Verifying reference files..."
ALL_REF_FILES=("${REF_FASTA}" "${REF_FAI}" "${REF_DICT}")
MISSING_FILES=false

for file in "${ALL_REF_FILES[@]}"; do
    if [[ ! -f "${file}" ]]; then
        log_error "  ✗ $(basename ${file}) - MISSING"
        MISSING_FILES=true
    else
        log_info "  ✓ $(basename ${file})"
    fi
done

if [[ "${MISSING_FILES}" == true ]]; then
    log_error "Some reference files are missing!"
    exit 1
fi

REF_SIZE=$(stat -c%s "${REF_FASTA}")
REF_SEQS=$(grep -c ">" "${REF_FASTA}")
log_info "  Reference size: ${REF_SIZE} bytes"
log_info "  Number of sequences: ${REF_SEQS}"

log_info "===== Setup Complete ====="