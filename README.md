# Variant Calling Benchmarking Pipeline 

Pipeline so sánh hiệu suất 4 variant caller: GATK, DeepVariant, Strelka2, FreeBayes.

## Quy trình thực hiện

### Phần I. Tạo cấu trúc folder

```bash
git clone https://github.com/bitschif/variant-benchmarking.git
cd variant-benchmarking/
```

```bash
mkdir -p data/reference
mkdir -p data/simulated
mkdir -p results/preprocessing
mkdir -p results/variants/gatk
mkdir -p results/variants/deepvariant
mkdir -p results/variants/strelka2
mkdir -p results/variants/freebayes
mkdir -p logs
```

### Phần II. Download và chuẩn bị dữ liệu

#### 2.1. Download Reference Genome (chr22 - hg38)

```bash
#pwd: variant-benchmarking/data/reference

# download từ UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# giải nén
gunzip chr22.fa.gz

# index reference
samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

#### 2.1b. (Tuỳ chọn) Tạo BED file exome cho WES

```bash
# pwd: variant-benchmarking/data/reference

# Download GENCODE annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# Extract chr22 exons → BED → merge → padding 100bp
grep "^chr22" gencode.v44.annotation.gtf | \
  awk '$3=="exon" {print $1"\t"$4-1"\t"$5}' | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - > chr22_exome_targets.bed

bedtools slop -i chr22_exome_targets.bed -g chr22.fa.fai -b 100 | \
  bedtools merge -i - > chr22_exome_targets_padded.bed

# Cho Strelka2 (cần bgzip + tabix)
bgzip -c chr22_exome_targets_padded.bed > chr22_exome_targets_padded.bed.gz
tabix -p bed chr22_exome_targets_padded.bed.gz
```

#### 2.2. Download Known Sites cho BQSR

Known sites giúp cải thiện chất lượng Base Quality Score Recalibration.

```bash
# pwd: variant-benchmarking/data/reference

# === dbSNP ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Extract chr22 và đảm bảo chromosome naming là "chr22"
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz

bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

# === Mills and 1000G Gold Standard Indels ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

# === 1000G Phase 1 SNPs ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 1000G_phase1.snps.high_confidence.hg38.vcf.gz -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz

# xoá file không sử dụng đến
rm -f Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx
rm -f Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
rm -f 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

**Một số lưu ý:**

```bash
# pwd = variant-benchmarking/data/reference

# reference genome từ UCSC sử dụng format "chr22"
# tất cả VCF files phải có tên chromosome khớp với tên ở file reference

# lệnh kiểm tra tên chromosome trong các file VCF
bcftools view -h file.vcf.gz | grep "^##contig"
tabix -l gile.vcf.gz # ở đồ án và repository này lấy tên là "chr22" làm chuẩn

# nếu gặp lỗi "chromosome not found", kiểm tra và đổi lại tên chromosome
# cách đổi tên "22" thành "chr22" nếu cần
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz

# tạo file BED chứa vùng non-N của chromosome
python3 << 'EOF'
from Bio import SeqIO
import re

fasta = "chr22.fa"
output = "chr22_non_N_regions.bed"

with open(output, 'w') as out:
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq).upper()
        chrom = record.id
        for m in re.finditer(r'[ATGC]+', seq):
            start = m.start()
            end = m.end()
            if end - start >= 1000:
                out.write(f"{chrom}\t{start}\t{end}\n")

print("Done")
EOF
```

#### 2.3. Tạo đột biến SNPs và INDELs với simuG

```bash
# pwd: variant-benchmarking/data

# clone simuG về folder <data>
git clone https://github.com/yjx1217/simuG.git

# dùng simuG tạo đột biến
perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 7000 \
  -indel_count 3500 \
  -indel_min_len 1 \
  -indel_max_len 5 \
  -prefix simulated/SIMULATED_SAMPLE_chr22
```

#### 2.4. Merge các file VCF output của simuG tạo truth VCF

```bash
# pwd: variant-benchmarking/data

# đảm bảo có .fai từ reference
samtools faidx reference/chr22.fa

# reheader từng file VCF với contig info
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf

bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf

# concat sau khi đã có header chuẩn
bcftools concat \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf | \
bcftools sort -Oz -o simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

tabix -p vcf simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

```

#### 2.5. Tạo file fastq giả lập bằng ART Illumina

```bash
#pwd: variant-benchmarking/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

art_illumina \
    -ss HS25 \
    -i "${MUTATED_FASTA}" \
    -p \
    -l 150 \
    -f 60 \
    -m 350 \
    -s 50 \
    -rs 42 \
    -o "${SIM_DIR}/${PREFIX}_" \
    -na

# đổi tên và nén
mv "${SIM_DIR}/${PREFIX}_1.fq" "${SIM_DIR}/${PREFIX}_R1.fastq"
mv "${SIM_DIR}/${PREFIX}_2.fq" "${SIM_DIR}/${PREFIX}_R2.fastq"
gzip "${SIM_DIR}/${PREFIX}_R1.fastq"
gzip "${SIM_DIR}/${PREFIX}_R2.fastq"
```

### Phần III. Tiền xử lý và gọi biến thể

Sau khi hoàn thành PHẦN 1 và PHẦN 2, chạy các script theo thứ tự:

```bash
#pwd: variant-benchmarking

bash 02_preprocessing.sh
bash 03_variant_calling_gatk.sh
bash 04_variant_calling_deepvariant.sh
bash 05_variant_calling_strelka2.sh
bash 06_variant_calling_freebayes.sh
```

## Cấu trúc thư mục

```
.
├── config/
│   └── config.sh
├── scripts/
│   └── helper_functions.sh
├── data/
│   ├── reference/
│   │   ├── chr22.fa
│   │   └── chr22_non_N_regions.bed
│   └── simulated/
│       ├── *_truth.vcf.gz
│       ├── *_R1.fastq.gz
│       └── *_R2.fastq.gz
├── results/
│   ├── preprocessing/
│   └── variants/
│       ├── gatk/
│       ├── deepvariant/
│       ├── strelka2/
│       └── freebayes/
├── 02_preprocessing.sh
├── 03_variant_calling_gatk.sh
├── 04_variant_calling_deepvariant.sh
├── 05_variant_calling_strelka2.sh
├── 06_variant_calling_freebayes.sh
└── README.md
```

## Các tools được sử dụng:

- fastq, bwa, samtools, bcftools
- simuG, art_illumina
- gatk, freebayes, docker (cho DeepVariant và Strelka2)
