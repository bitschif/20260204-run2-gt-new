#!/bin/bash
#===============================================================================
# Helper script: Rename chromosomes from "22" to "chr22" format
# Usage: ./rename_chromosomes.sh input.vcf.gz output.vcf.gz
#===============================================================================

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"

if [[ -z "${INPUT_VCF}" || -z "${OUTPUT_VCF}" ]]; then
    echo "Usage: $0 <input.vcf.gz> <output.vcf.gz>"
    exit 1
fi

# Create chromosome mapping file
CHR_MAP=$(mktemp)
for i in {1..22} X Y M MT; do
    echo -e "${i}\tchr${i}" >> "${CHR_MAP}"
done
# Handle MT -> chrM
echo -e "MT\tchrM" >> "${CHR_MAP}"

# Check if renaming is needed
FIRST_CHR=$(bcftools view -H "${INPUT_VCF}" | head -1 | cut -f1)

if [[ "${FIRST_CHR}" == chr* ]]; then
    echo "Chromosomes already in 'chr' format, copying file..."
    cp "${INPUT_VCF}" "${OUTPUT_VCF}"
else
    echo "Renaming chromosomes from '${FIRST_CHR}' format to 'chr' format..."
    bcftools annotate --rename-chrs "${CHR_MAP}" "${INPUT_VCF}" -Oz -o "${OUTPUT_VCF}"
fi

# Index output
tabix -f -p vcf "${OUTPUT_VCF}"

# Cleanup
rm -f "${CHR_MAP}"

echo "Done: ${OUTPUT_VCF}"
