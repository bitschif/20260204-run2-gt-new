#!/bin/bash
#===============================================================================
# CONFIG: Cấu hình chung cho pipeline variant calling benchmarking
# Theo GATK Best Practices (nf-core/sarek)
#===============================================================================

#-------------------------------------------------------------------------------
# SYSTEM RESOURCES (RAM 16GB, 4 CPU)
#-------------------------------------------------------------------------------
export THREADS=4
export MAX_MEMORY="14G"
export JAVA_OPTS="-Xmx12G -XX:ParallelGCThreads=2"

#-------------------------------------------------------------------------------
# DIRECTORY PATHS
#-------------------------------------------------------------------------------
export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
export DATA_DIR="${PROJECT_DIR}/data"
export REF_DIR="${DATA_DIR}/reference"
export SIM_DIR="${DATA_DIR}/simulated"
export RESULTS_DIR="${PROJECT_DIR}/results"
export LOG_DIR="${PROJECT_DIR}/logs"

export PREPROC_DIR="${RESULTS_DIR}/preprocessing"
export VARIANT_DIR="${RESULTS_DIR}/variants"
export BENCH_DIR="${RESULTS_DIR}/benchmarks"
export FIGURE_DIR="${RESULTS_DIR}/plots"
export METRICS_DIR="${RESULTS_DIR}/final_metrics"

#-------------------------------------------------------------------------------
# REFERENCE GENOME
#-------------------------------------------------------------------------------
export GENOME_VERSION="hg38"
export CHR_TO_USE="chr22"
export REF_FASTA="${REF_DIR}/${CHR_TO_USE}.fa"
export REF_DICT="${REF_DIR}/${CHR_TO_USE}.dict"
export REF_FAI="${REF_DIR}/${CHR_TO_USE}.fa.fai"

# Known sites for BQSR (GATK Bundle)
export DBSNP="${REF_DIR}/dbsnp_146.hg38.chr22.vcf.gz"
export KNOWN_INDELS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz"
export KNOWN_SNPS="${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz"

#-------------------------------------------------------------------------------
# SIMULATION PARAMETERS (simuG + ART)
#-------------------------------------------------------------------------------
# simuG parameters for mutation simulation
export SNP_COUNT=7000          # Number of SNPs to simulate
export INDEL_COUNT=3500        # Number of indels to simulate (combined del+ins)
export INDEL_MIN_LEN=1         # Minimum indel length
export INDEL_MAX_LEN=5         # Maximum indel length

# ART Illumina parameters
export COVERAGE=60
export READ_LENGTH=150
export FRAGMENT_MEAN=350
export FRAGMENT_SD=50
export ART_PLATFORM="HS25"
export SEED=42

#-------------------------------------------------------------------------------
# QUALITY THRESHOLDS
#-------------------------------------------------------------------------------
export MIN_BASE_QUALITY=20
export MIN_MAPPING_QUALITY=20
export MIN_READ_LENGTH=50

#-------------------------------------------------------------------------------
# SAMPLE INFO
#-------------------------------------------------------------------------------
export SAMPLE_NAME="SIMULATED_SAMPLE"
export PREFIX="${SAMPLE_NAME}_${CHR_TO_USE}"
export READ_GROUP="@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"

#-------------------------------------------------------------------------------
# DOCKER IMAGES
#-------------------------------------------------------------------------------
export DEEPVARIANT_VERSION="1.6.1"
export DEEPVARIANT_IMAGE="google/deepvariant:${DEEPVARIANT_VERSION}"
export STRELKA2_IMAGE="quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"

#-------------------------------------------------------------------------------
# VARIANT CALLER PARAMETERS
#-------------------------------------------------------------------------------
# GATK HaplotypeCaller
export GATK_STAND_CALL_CONF=30

# FreeBayes
export FB_MIN_ALT_COUNT=3
export FB_MIN_ALT_FRACTION=0.2

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------
export TRUTH_VCF="${SIM_DIR}/${PREFIX}_truth.vcf.gz"
export HIGH_CONF_BED="${REF_DIR}/chr22_non_N_regions.bed"

echo "[CONFIG] Loaded successfully"
echo "[CONFIG] Project: ${PROJECT_DIR}"
echo "[CONFIG] Chromosome: ${CHR_TO_USE}, Coverage: ${COVERAGE}x"
