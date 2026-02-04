# Variant Calling Benchmarking Pipeline 

Pipeline so sánh hiệu suất 4 variant caller: GATK, DeepVariant, Strelka2, FreeBayes

## Quy trình thực hiện

### PHẦN 1: Tạo cấu trúc folder

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

### PHẦN 2: Download và chuẩn bị dữ liệu

#### 2.1. Download Reference Genome (chr22 - hg38)

```bash
#pwd: variant-benchmarking/data/reference

# Download từ UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# Hoặc dùng curl
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# Giải nén
gunzip chr22.fa.gz

# Index reference
samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

#### 2.2. Download Known Sites cho BQSR (Tùy chọn nhưng khuyến nghị)

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

# Extract chr22
bcftools view -r chr22 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

# === 1000G Phase 1 SNPs ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# Extract chr22
bcftools view -r chr22 1000G_phase1.snps.high_confidence.hg38.vcf.gz -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz

# Dọn dẹp file gốc (tùy chọn)
rm -f Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx
rm -f Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
rm -f 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

**Một số lưu ý:**

```bash
# pwd = variant-benchmarking/data/reference
# Reference genome từ UCSC sử dụng format "chr22"
# Tất cả VCF files phải có tên chromosome khớp với tên ở file reference
# Kiểm tra tên chromosome trong VCF
bcftools view -h file.vcf.gz | grep "^##contig"
tabix -l gile.vcf.gz # ở đồ án và repository này lấy tên là "chr22" làm chuẩn

# Nếu gặp lỗi "chromosome not found", kiểm tra và đổi lại chromosome.
# Cách đổi tên "22" thành "chr22" nếu cần
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz

# Tạo file BED chứa vùng non-N của chromosome
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

print("Done!")
EOF
```

#### 2.3. Tạo dữ liệu giả lập với Simutator

```bash
# pwd: variant-benchmarking/data

REF_FASTA="/reference/chr22.fa"
SIM_DIR="/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"

simutator mutate_fasta \
    --snps 7000 \
    --dels 2000:3 \
    --ins 2000:2 \
    --seed 42 \
    "${REF_FASTA}" \
    "${SIM_DIR}/${PREFIX}"
```

#### 2.4. Merge và chuẩn bị Truth VCF

```bash
#pwd: variant-benchmarking/data

REF_FASTA="/reference/chr22.fa"
SIM_DIR="/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
TRUTH_VCF="${SIM_DIR}/${PREFIX}_truth.vcf.gz"

# Merge tất cả VCF
bcftools concat ${SIM_DIR}/${PREFIX}*.original.vcf | bcftools sort -Oz -o ${TRUTH_VCF}
tabix -p vcf ${TRUTH_VCF}

# Tạo mutated FASTA
SAMPLE_NAME=$(bcftools query -l ${TRUTH_VCF} | head -n 1)
bcftools consensus -s "${SAMPLE_NAME}" -f "${REF_FASTA}" "${TRUTH_VCF}" > "${SIM_DIR}/${PREFIX}_mutated_combined.fa"
samtools faidx "${SIM_DIR}/${PREFIX}_mutated_combined.fa"

# Tạo file SNP và INDEL riêng
bcftools +fill-tags ${TRUTH_VCF} -Oz -o "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -- -t TYPE
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz"

bcftools view -i 'TYPE="snp"' "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -Oz -o "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
bcftools view -i 'TYPE="indel"' "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -Oz -o "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
```

#### 2.5. Tạo reads với ART Illumina

```bash
#pwd: variant-benchmarking/data

SIM_DIR="/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}_mutated_combined.fa"

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

# Đổi tên và nén
mv "${SIM_DIR}/${PREFIX}_1.fq" "${SIM_DIR}/${PREFIX}_R1.fastq"
mv "${SIM_DIR}/${PREFIX}_2.fq" "${SIM_DIR}/${PREFIX}_R2.fastq"
gzip "${SIM_DIR}/${PREFIX}_R1.fastq"
gzip "${SIM_DIR}/${PREFIX}_R2.fastq"
```

### PHẦN 3: Chạy Pipeline

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

## Yêu cầu công cụ

- fastq, bwa, samtools, bcftools
- simutator, art_illumina
- gatk, freebayes, docker (cho DeepVariant và Strelka2)
