#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking
# Compare variant calls against truth set
# Tools: hap.py (Illumina)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Benchmarking ====="
start_timer

check_tool "hap.py" || exit 1
if command -v python3 &>/dev/null; then
    PYTHON_BIN="python3"
elif command -v python &>/dev/null; then
    PYTHON_BIN="python"
else
    log_error "Python not found for parsing hap.py summary"
    exit 1
fi

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")
TRUTH_SNP="${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
TRUTH_INDEL="${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1
check_file "${REF_FASTA}" || exit 1

#-------------------------------------------------------------------------------
# Function: Benchmark with hap.py
#-------------------------------------------------------------------------------
benchmark_happy() {
    local caller=$1
    local query_vcf=$2
    local truth_vcf=$3
    local out_prefix=$4
    
    ensure_dir "$(dirname "${out_prefix}")"
    hap.py "${truth_vcf}" "${query_vcf}" \
        -f "${HIGH_CONF_BED}" \
        -r "${REF_FASTA}" \
        -o "${out_prefix}"
}

parse_happy_summary() {
    local summary_csv=$1
    local caller=$2

    "${PYTHON_BIN}" - "${summary_csv}" "${caller}" <<'PY'
import csv
import re
import sys

summary_csv, caller = sys.argv[1:3]

def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())

with open(summary_csv, newline="") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames is None:
        sys.exit(0)
    normalized = {norm(name): name for name in reader.fieldnames}

    def pick(*names):
        for name in names:
            key = norm(name)
            if key in normalized:
                return normalized[key]
        return None

    type_key = pick("type")
    tp_key = pick("tp")
    fp_key = pick("fp")
    fn_key = pick("fn")
    precision_key = pick("precision", "prec")
    recall_key = pick("recall")
    f1_key = pick("f1", "f1score")

    if not all([type_key, tp_key, fp_key, fn_key, precision_key, recall_key, f1_key]):
        sys.exit(0)

    wanted = {
        "ALL": {"all", "total", "overall"},
        "SNP": {"snp"},
        "INDEL": {"indel"},
    }

    for row in reader:
        row_type = row.get(type_key, "").strip().lower()
        for out_type, labels in wanted.items():
            if row_type in labels:
                print(
                    f"{caller}\t{out_type}\t"
                    f"{row.get(tp_key, '')}\t{row.get(fp_key, '')}\t{row.get(fn_key, '')}\t"
                    f"{row.get(precision_key, '')}\t{row.get(recall_key, '')}\t{row.get(f1_key, '')}"
                )
                break
PY
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
    
    log_info "  hap.py..."
    HAPPY_PREFIX="${BENCH_CALLER}/happy/${PREFIX}_${caller}"
    benchmark_happy "${caller}" "${QUERY_VCF}" "${TRUTH_VCF}" "${HAPPY_PREFIX}"

    SUMMARY_CSV="${HAPPY_PREFIX}.summary.csv"
    if [[ -f "${SUMMARY_CSV}" ]]; then
        parse_happy_summary "${SUMMARY_CSV}" "${caller}" >> "${SUMMARY}"
    else
        log_warn "  hap.py summary not found for ${caller}"
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
