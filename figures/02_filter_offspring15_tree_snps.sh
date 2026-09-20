#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"

RAW_VCF="${PROJECT}/figure/offspring_NJ/01_direct_call/offspring15.direct_diploid_SB.raw.vcf.gz"

OUT_DIR="${PROJECT}/07_phylogeny/figure1_offspring_NJ/03_filtered"
TMP_DIR="${PROJECT}/07_phylogeny/figure1_offspring_NJ/tmp"
LOG_DIR="${PROJECT}/logs/figure1_offspring_NJ"

mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${LOG_DIR}"

RAW_SNP="${OUT_DIR}/offspring15.biallelic.snps.raw.vcf.gz"
FILTERED_SNP="${OUT_DIR}/offspring15.biallelic.snps.hardfiltered.vcf.gz"
PASS_SNP="${OUT_DIR}/offspring15.biallelic.snps.PASS.vcf.gz"
MASKED_SNP="${OUT_DIR}/offspring15.biallelic.snps.PASS.masked.vcf.gz"
FINAL_VCF="${OUT_DIR}/offspring15.tree_ready.complete_polymorphic_snps.vcf.gz"

JAVA_MEM="60g"
THREADS=32

MIN_GQ=30
MIN_DP=10
MAX_DP=80


# ============================================================
# 1. Extract biallelic SNPs
# ============================================================

gatk \
    --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" \
    SelectVariants \
    -R "${REF}" \
    -V "${RAW_VCF}" \
    -O "${RAW_SNP}" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    > "${LOG_DIR}/offspring15.SelectVariants.log" 2>&1


# ============================================================
# 2. Site-level hard filtering
# ============================================================

gatk \
    --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" \
    VariantFiltration \
    -R "${REF}" \
    -V "${RAW_SNP}" \
    -O "${FILTERED_SNP}" \
    --filter-name "QUAL30" \
    --filter-expression "QUAL < 30.0" \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS60" \
    --filter-expression "FS > 60.0" \
    --filter-name "SOR3" \
    --filter-expression "SOR > 3.0" \
    --filter-name "MQ40" \
    --filter-expression "MQ < 40.0" \
    --filter-name "MQRankSum-12.5" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum-8" \
    --filter-expression "ReadPosRankSum < -8.0" \
    > "${LOG_DIR}/offspring15.VariantFiltration.log" 2>&1


# ============================================================
# 3. Keep PASS SNPs
# ============================================================

bcftools view \
    --threads "${THREADS}" \
    -f PASS \
    -Oz \
    -o "${PASS_SNP}" \
    "${FILTERED_SNP}"

bcftools index \
    --threads "${THREADS}" \
    -f -t "${PASS_SNP}"


# ============================================================
# 4. Mask unreliable genotypes
# ============================================================

# Mask:
#   heterozygous genotype
#   GQ < 30
#   DP < 10
#   DP > 80
#
# Because offspring are biologically haploid, 0/1 is treated
# as possible mapping, CNV or sequencing artefact.

bcftools +setGT \
    "${PASS_SNP}" \
    -Oz \
    -o "${MASKED_SNP}" \
    -- \
    -t q \
    -n . \
    -i "GT=\"het\" || FMT/GQ<${MIN_GQ} || FMT/DP<${MIN_DP} || FMT/DP>${MAX_DP}"

bcftools index -f -t "${MASKED_SNP}"


# ============================================================
# 5. Recalculate AC/AN/missingness and retain complete sites
# ============================================================

# Requirements:
#   no missing genotype among 15 samples
#   both reference and alternative alleles are represented
#   therefore the site remains polymorphic after masking

bcftools +fill-tags \
    "${MASKED_SNP}" \
    -Ou \
    -- \
    -t AC,AN,NS \
|
bcftools view \
    -i 'INFO/AN=30 && INFO/NS=15 && INFO/AC>0 && INFO/AC<INFO/AN' \
    -Oz \
    -o "${FINAL_VCF}"

bcftools index \
    --threads "${THREADS}" \
    -f -t "${FINAL_VCF}"


# ============================================================
# 6. Summaries
# ============================================================

N_RAW=$(bcftools view -H "${RAW_VCF}" | wc -l)
N_BIALLELIC=$(bcftools view -H "${RAW_SNP}" | wc -l)
N_PASS=$(bcftools view -H "${PASS_SNP}" | wc -l)
N_FINAL=$(bcftools view -H "${FINAL_VCF}" | wc -l)
N_HET=$(bcftools view -H -g het "${FINAL_VCF}" | wc -l)
N_MISSING=$(bcftools view -H -g miss "${FINAL_VCF}" | wc -l)
N_SAMPLE=$(bcftools query -l "${FINAL_VCF}" | wc -l)

SUMMARY="${OUT_DIR}/offspring15.tree_SNP_filtering_summary.txt"

{
    echo "Raw variants:                  ${N_RAW}"
    echo "Raw biallelic SNPs:            ${N_BIALLELIC}"
    echo "PASS biallelic SNPs:           ${N_PASS}"
    echo "Final complete polymorphic:    ${N_FINAL}"
    echo "Sites containing heterozygote: ${N_HET}"
    echo "Sites containing missing GT:   ${N_MISSING}"
    echo "Samples:                       ${N_SAMPLE}"
    echo
    echo "Final VCF:"
    echo "${FINAL_VCF}"
} | tee "${SUMMARY}"

bcftools stats \
    "${FINAL_VCF}" \
    > "${OUT_DIR}/offspring15.tree_ready.bcftools_stats.txt"
