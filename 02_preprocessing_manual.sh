#!/bin/bash
#===============================================================================
# STEP 02: Preprocessing (GATK Best Practices - nf-core/sarek)
# Manual commands without variable expansion.
# Pipeline: FastQC → fastp → BWA-MEM → MarkDuplicates → BQSR
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# 1. FastQC - Raw reads quality control
#-------------------------------------------------------------------------------
mkdir -p results/preprocessing/fastqc_raw
fastqc -t 4 \
  -o results/preprocessing/fastqc_raw \
  data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  2>&1 | tee logs/fastqc_raw.log

#-------------------------------------------------------------------------------
# 2. fastp - Trimming and filtering
#-------------------------------------------------------------------------------
fastp \
  -i data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  -I data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  -o results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R1.fastq.gz \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R2.fastq.gz \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --cut_front --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --thread 4 \
  --json results/preprocessing/SIMULATED_SAMPLE_chr22_fastp.json \
  --html results/preprocessing/SIMULATED_SAMPLE_chr22_fastp.html \
  2>&1 | tee logs/fastp.log

#-------------------------------------------------------------------------------
# 3. FastQC - Trimmed reads
#-------------------------------------------------------------------------------
mkdir -p results/preprocessing/fastqc_trimmed
fastqc -t 4 \
  -o results/preprocessing/fastqc_trimmed \
  results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R1.fastq.gz \
  results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R2.fastq.gz \
  2>&1 | tee logs/fastqc_trimmed.log

#-------------------------------------------------------------------------------
# 4. BWA-MEM - Alignment
#-------------------------------------------------------------------------------
bwa mem \
  -t 4 \
  -R "@RG\tID:SIMULATED_SAMPLE\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  -M \
  data/reference/chr22.fa \
  results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R1.fastq.gz \
  results/preprocessing/SIMULATED_SAMPLE_chr22_trimmed_R2.fastq.gz \
  2> logs/bwa_mem.log | \
  samtools sort -@ 4 -m 2G -o results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam

#-------------------------------------------------------------------------------
# 5. GATK MarkDuplicates
#-------------------------------------------------------------------------------
gatk MarkDuplicates \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  -M results/preprocessing/SIMULATED_SAMPLE_chr22_dup_metrics.txt \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX true \
  2>&1 | tee logs/markduplicates.log

#-------------------------------------------------------------------------------
# 6. GATK BaseRecalibrator & ApplyBQSR
#-------------------------------------------------------------------------------
# Ensure known-sites VCFs match the "chr" naming convention if needed.
# If your VCFs use "22" instead of "chr22", run:
# scripts/rename_chromosomes.sh <input.vcf.gz> <output.vcf.gz>

gatk BaseRecalibrator \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --known-sites data/reference/dbsnp_146.hg38.chr22.vcf.gz \
  --known-sites data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  2>&1 | tee logs/baserecalibrator.log

gatk ApplyBQSR \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --bqsr-recal-file results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam \
  2>&1 | tee logs/applybqsr.log

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam

#-------------------------------------------------------------------------------
# 7. Filter by mapping quality
#-------------------------------------------------------------------------------
samtools view \
  -@ 4 \
  -b \
  -q 20 \
  -F 1796 \
  results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam | \
  samtools sort -@ 4 -o results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam

#-------------------------------------------------------------------------------
# 8. Alignment statistics (samtools stats, mosdepth)
#-------------------------------------------------------------------------------
samtools stats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_stats.txt
samtools flagstat results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_flagstat.txt
samtools idxstats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_idxstats.txt

# Optional: mosdepth for coverage (if installed)
# mosdepth -t 4 --by 1000 \
#   results/preprocessing/SIMULATED_SAMPLE_chr22_coverage \
#   results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam

#-------------------------------------------------------------------------------
# 9. Export for downstream scripts
#-------------------------------------------------------------------------------
echo "FINAL_BAM=results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam" \
  > results/preprocessing/bam_path.sh
