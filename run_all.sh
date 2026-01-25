#!/bin/bash
#===============================================================================
# MASTER SCRIPT:  Run entire variant calling benchmarking pipeline
# Usage: ./run_all.sh [OPTIONS]
#
# Options:
#   --skip-setup         Skip environment setup
#   --skip-simulation    Skip data simulation
#   --skip-preprocessing Skip preprocessing
#   --skip-calling       Skip variant calling
#   --callers LIST       Comma-separated callers (gatk,deepvariant,strelka2,freebayes)
#   -h, --help           Show help
#===============================================================================

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="${SCRIPT_DIR}"

# Default options
SKIP_SETUP=false
SKIP_SIMULATION=false
SKIP_PREPROCESSING=false
SKIP_CALLING=false
RUN_GATK=true
RUN_DEEPVARIANT=true
RUN_STRELKA2=true
RUN_FREEBAYES=true

#-------------------------------------------------------------------------------
# Parse arguments
#-------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-setup)
            SKIP_SETUP=true
            shift
            ;;
        --skip-simulation)
            SKIP_SIMULATION=true
            shift
            ;;
        --skip-preprocessing)
            SKIP_PREPROCESSING=true
            shift
            ;;
        --skip-calling)
            SKIP_CALLING=true
            shift
            ;;
        --callers)
            RUN_GATK=false
            RUN_DEEPVARIANT=false
            RUN_STRELKA2=false
            RUN_FREEBAYES=false
            IFS=',' read -ra CALLERS <<< "$2"
            for c in "${CALLERS[@]}"; do
                case "${c}" in
                    gatk) RUN_GATK=true ;;
                    deepvariant) RUN_DEEPVARIANT=true ;;
                    strelka2) RUN_STRELKA2=true ;;
                    freebayes) RUN_FREEBAYES=true ;;
                    *) echo "Unknown caller: ${c}"; exit 1 ;;
                esac
            done
            shift 2
            ;;
        -h|--help)
            head -20 "$0" | tail -15
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

#-------------------------------------------------------------------------------
# Load configuration
#-------------------------------------------------------------------------------
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

#-------------------------------------------------------------------------------
# Print banner
#-------------------------------------------------------------------------------
echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║     VARIANT CALLING BENCHMARKING PIPELINE                          ║"
echo "║     GATK | DeepVariant | Strelka2 | FreeBayes                      ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Start time: $(date)"
echo "Project directory: ${PROJECT_DIR}"
echo ""

TOTAL_START=$(date +%s)

#-------------------------------------------------------------------------------
# Step 0: Setup Environment
#-------------------------------------------------------------------------------
if [[ "${SKIP_SETUP}" == false ]]; then
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│ STEP 0/6: Setting up environment                                 │"
    echo "└──────────────────────────────────────────────────────────────────┘"
    bash "${SCRIPT_DIR}/00_setup_environment.sh"
else
    echo "[SKIP] Step 0: Setup environment"
fi

#-------------------------------------------------------------------------------
# Step 1: Simulate Data
#-------------------------------------------------------------------------------
if [[ "${SKIP_SIMULATION}" == false ]]; then
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│ STEP 1/6: Simulating mutations and reads                         │"
    echo "└──────────────────────────────────────────────────────────────────┘"
    bash "${SCRIPT_DIR}/01_simulate_data.sh"
else
    echo "[SKIP] Step 1: Simulation"
fi

#-------------------------------------------------------------------------------
# Step 2: Preprocessing
#-------------------------------------------------------------------------------
if [[ "${SKIP_PREPROCESSING}" == false ]]; then
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│ STEP 2/6: Preprocessing (GATK Best Practices)                    │"
    echo "└──────────────────────────────────────────────────────────────────┘"
    bash "${SCRIPT_DIR}/02_preprocessing.sh"
else
    echo "[SKIP] Step 2: Preprocessing"
fi

#-------------------------------------------------------------------------------
# Steps 3-6: Variant Calling
#-------------------------------------------------------------------------------
if [[ "${SKIP_CALLING}" == false ]]; then
    
    # GATK
    if [[ "${RUN_GATK}" == true ]]; then
        echo ""
        echo "┌──────────────────────────────────────────────────────────────────┐"
        echo "│ STEP 3/6: GATK HaplotypeCaller                                   │"
        echo "└──────────────────────────────────────────────────────────────────┘"
        bash "${SCRIPT_DIR}/03_variant_calling_gatk.sh"
    else
        echo "[SKIP] Step 3: GATK"
    fi
    
    # DeepVariant
    if [[ "${RUN_DEEPVARIANT}" == true ]]; then
        echo ""
        echo "┌──────────────────────────────────────────────────────────────────┐"
        echo "│ STEP 4/6: DeepVariant (Docker)                                   │"
        echo "└──────────────────────────────────────────────────────────────────┘"
        bash "${SCRIPT_DIR}/04_variant_calling_deepvariant.sh"
    else
        echo "[SKIP] Step 4: DeepVariant"
    fi
    
    # Strelka2
    if [[ "${RUN_STRELKA2}" == true ]]; then
        echo ""
        echo "┌──────────────────────────────────────────────────────────────────┐"
        echo "│ STEP 5/6: Strelka2 Germline (Docker)                             │"
        echo "└──────────────────────────────────────────────────────────────────┘"
        bash "${SCRIPT_DIR}/05_variant_calling_strelka2.sh"
    else
        echo "[SKIP] Step 5: Strelka2"
    fi
    
    # FreeBayes
    if [[ "${RUN_FREEBAYES}" == true ]]; then
        echo ""
        echo "┌──────────────────────────────────────────────────────────────────┐"
        echo "│ STEP 6/6: FreeBayes                                              │"
        echo "└──────────────────────────────────────────────────────────────────┘"
        bash "${SCRIPT_DIR}/06_variant_calling_freebayes.sh"
    else
        echo "[SKIP] Step 6: FreeBayes"
    fi
    
else
    echo "[SKIP] Steps 3-6: Variant calling"
fi

#-------------------------------------------------------------------------------
# Step 9: R Visualization
#-------------------------------------------------------------------------------
echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│ BONUS: Generating visualizations with R                          │"
echo "└──────────────────────────────────────────────────────────────────┘"

if command -v Rscript &> /dev/null; then
    Rscript "${SCRIPT_DIR}/R/visualize_benchmark.R" "${PROJECT_DIR}" 2>&1 || {
        echo "[WARN] R visualization failed, but pipeline completed"
    }
else
    echo "[WARN] Rscript not found.  Run manually:"
    echo "  Rscript R/visualize_benchmark.R ${PROJECT_DIR}"
fi

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
TOTAL_END=$(date +%s)
TOTAL_DURATION=$((TOTAL_END - TOTAL_START))
HOURS=$((TOTAL_DURATION / 3600))
MINUTES=$(((TOTAL_DURATION % 3600) / 60))
SECONDS=$((TOTAL_DURATION % 60))

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                     PIPELINE COMPLETE!                              ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""
echo "Output locations:"
echo "  ├── Simulated data:   ${SIM_DIR}/"
echo "  ├── Preprocessed:     ${PREPROC_DIR}/"
echo "  ├── Variants:        ${VARIANT_DIR}/"
echo "  ├── Benchmarks:      ${BENCH_DIR}/"
echo "  └── Figures:         ${FIGURE_DIR}/"
echo ""
echo "Key files:"
echo "  ├── Metrics TSV:     ${RESULTS_DIR}/metrics_long.tsv"
echo "  ├── Runtime TSV:     ${RESULTS_DIR}/runtime.tsv"
echo "  └── Runtime log:     ${LOG_DIR}/runtime.csv"
echo ""
echo "End time: $(date)"
echo ""

#-------------------------------------------------------------------------------
# Print final benchmark results
#-------------------------------------------------------------------------------
echo "═══════════════════════════════════════════════════════════════════════"
echo "FINAL BENCHMARK RESULTS"
echo "═══════════════════════════════════════════════════════════════════════"
if [[ -f "${RESULTS_DIR}/metrics_long.tsv" ]]; then
    column -t -s$'\t' "${RESULTS_DIR}/metrics_long.tsv"
fi
echo "═══════════════════════════════════════════════════════════════════════"
