#!/bin/bash
#===============================================================================
# STEP 00: Setup Environment
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 00: Setup Environment ====="

#-------------------------------------------------------------------------------
# 1. Create directories
#-------------------------------------------------------------------------------
for dir in "${DATA_DIR}" "${REF_DIR}" "${SIM_DIR}" "${RESULTS_DIR}" \
           "${LOG_DIR}" "${PREPROC_DIR}" "${VARIANT_DIR}" "${BENCH_DIR}" \
           "${FIGURE_DIR}" "${METRICS_DIR}"; do
    ensure_dir "$dir"
done

for caller in gatk deepvariant strelka2 freebayes; do
    ensure_dir "${VARIANT_DIR}/${caller}"
    ensure_dir "${BENCH_DIR}/${caller}"
done

# Runtime log header
echo "step,duration_seconds" > "${LOG_DIR}/runtime. csv"

#-------------------------------------------------------------------------------
# 2. Check required tools
#-------------------------------------------------------------------------------
log_info "Checking required tools..."

for tool in bwa samtools bcftools gatk fastp fastqc bgzip tabix art_illumina freebayes simutator; do
    if command -v "$tool" &> /dev/null; then
        log_info "  ✓ $tool"
    else
        log_warn "  ✗ $tool (missing)"
    fi
done

# Check Docker
if command -v docker &> /dev/null; then
    log_info "  ✓ docker (for DeepVariant, Strelka2)"
else
    log_warn "  ✗ docker (DeepVariant, Strelka2 sẽ không chạy được)"
fi

#-------------------------------------------------------------------------------
# 3. Download reference genome - FIX ERROR HANDLING
#-------------------------------------------------------------------------------
if [[ !  -f "${REF_FASTA}" ]]; then
    log_info "Downloading ${CHR_TO_USE} from UCSC..."
    
    cd "${REF_DIR}"
    
    # URL cho hg38 chromosomes từ UCSC
    UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHR_TO_USE}. fa.gz"
    
    # Alternative URLs nếu UCSC không hoạt động
    NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr22.fna.gz"
    ENSEMBL_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens. GRCh38.dna.chromosome.22.fa.gz"
    
    DOWNLOAD_SUCCESS=false
    
    # Try UCSC first
    log_info "  Trying UCSC..."
    if wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa. gz" "${UCSC_URL}" 2>/dev/null; then
        DOWNLOAD_SUCCESS=true
        log_info "  Downloaded from UCSC"
    else
        log_warn "  UCSC download failed, trying Ensembl..."
        
        # Try Ensembl
        if wget -q --show-progress --timeout=60 -O "${CHR_TO_USE}.fa.gz" "${ENSEMBL_URL}" 2>/dev/null; then
            DOWNLOAD_SUCCESS=true
            log_info "  Downloaded from Ensembl"
        else
            log_warn "  Ensembl download failed, trying alternative method..."
            
            # Try with curl
            if curl -L -f -o "${CHR_TO_USE}.fa.gz" "${UCSC_URL}" 2>/dev/null; then
                DOWNLOAD_SUCCESS=true
                log_info "  Downloaded with curl"
            fi
        fi
    fi
    
    # Check if download was successful
    if [[ "${DOWNLOAD_SUCCESS}" == false ]] || [[ ! -f "${CHR_TO_USE}.fa.gz" ]]; then
        log_error "Failed to download reference genome!"
        log_info ""
        log_info "Please download manually:"
        log_info "  Option 1 (UCSC):"
        log_info "    wget ${UCSC_URL}"
        log_info ""
        log_info "  Option 2 (Ensembl):"
        log_info "    wget ${ENSEMBL_URL}"
        log_info ""
        log_info "  Then place the file in:  ${REF_DIR}/"
        log_info "  And run: gunzip ${CHR_TO_USE}.fa.gz"
        exit 1
    fi
    
    # Check file size (should be > 1MB for chr22)
    FILE_SIZE=$(stat -f%z "${CHR_TO_USE}.fa.gz" 2>/dev/null || stat -c%s "${CHR_TO_USE}.fa.gz" 2>/dev/null || echo "0")
    if [[ "${FILE_SIZE}" -lt 1000000 ]]; then
        log_error "Downloaded file is too small (${FILE_SIZE} bytes). Download may have failed."
        rm -f "${CHR_TO_USE}.fa.gz"
        exit 1
    fi
    
    log_info "  File size: $(ls -lh ${CHR_TO_USE}.fa.gz | awk '{print $5}')"
    
    # Decompress
    log_info "  Decompressing..."
    gunzip -f "${CHR_TO_USE}.fa.gz"
    
    if [[ ! -f "${REF_FASTA}" ]]; then
        log_error "Decompression failed!"
        exit 1
    fi
    
    # Index reference
    log_info "  Indexing reference with samtools..."
    samtools faidx "${REF_FASTA}"
    
    log_info "  Indexing reference with BWA..."
    bwa index "${REF_FASTA}"
    
    log_info "  Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}" 2>/dev/null || \
    samtools dict "${REF_FASTA}" > "${REF_DICT}"
    
    log_info "  Reference preparation complete!"
    
else
    log_info "Reference already exists:  ${REF_FASTA}"
    
    # Check if indexes exist
    if [[ ! -f "${REF_FAI}" ]]; then
        log_info "  Creating fasta index..."
        samtools faidx "${REF_FASTA}"
    fi
    
    if [[ ! -f "${REF_FASTA}.bwt" ]]; then
        log_info "  Creating BWA index..."
        bwa index "${REF_FASTA}"
    fi
    
    if [[ ! -f "${REF_DICT}" ]]; then
        log_info "  Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}" 2>/dev/null || \
        samtools dict "${REF_FASTA}" > "${REF_DICT}"
    fi
fi

#-------------------------------------------------------------------------------
# 4. Verify reference files
#-------------------------------------------------------------------------------
log_info "Verifying reference files..."

MISSING_FILES=false

for file in "${REF_FASTA}" "${REF_FAI}"; do
    if [[ -f "${file}" ]]; then
        log_info "  ✓ $(basename ${file})"
    else
        log_error "  ✗ $(basename ${file}) - MISSING"
        MISSING_FILES=true
    fi
done

if [[ "${MISSING_FILES}" == true ]]; then
    log_error "Some reference files are missing!"
    exit 1
fi

# Show reference info
REF_SIZE=$(ls -lh "${REF_FASTA}" | awk '{print $5}')
REF_SEQS=$(grep -c "^>" "${REF_FASTA}")
log_info "  Reference size: ${REF_SIZE}"
log_info "  Number of sequences: ${REF_SEQS}"

log_info "===== Setup Complete ====="
