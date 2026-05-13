#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Structural variant liftover and benchmarking pipeline
#
# This script:
#   1. Lifts SVs from CHM13v2 to GRCh38 using CrossMap
#   2. Generates REF/ALT-swapped VCF records
#   3. Performs liftover on swapped alleles
#   4. Benchmarks both lifted VCFs against gnomAD-SV using Truvari
#
# Requirements:
#   - CrossMap
#   - bcftools
#   - bgzip/tabix (htslib)
#   - truvari
#
# Usage:
#   bash sv_crossmap_and_benchmark_GnomAD.sh \
#       chm13v2-grch38.chain \
#       GRCh38.p14.fasta \
#       511.candidate.sv.chm13.vcf.gz \
#       gnomad.v4.1.sv.sites.vcf.gz
#
###############################################################################

CHAIN=$1
REFERENCE=$2
INPUT_VCF=$3
GNOMAD_VCF=$4

PREFIX=$(basename ${INPUT_VCF} .vcf.gz)

###############################################################################
# Step 1. Direct liftover
###############################################################################

echo "[INFO] Running direct liftover..."

CrossMap vcf \
    ${CHAIN} \
    ${INPUT_VCF} \
    ${REFERENCE} \
    ${PREFIX}.lifted.hg38.vcf

bcftools view ${PREFIX}.lifted.hg38.vcf -Oz -o ${PREFIX}.lifted.hg38.vcf.gz -W=tbi

###############################################################################
# Step 2. Generate REF/ALT-swapped VCF
###############################################################################

echo "[INFO] Generating REF/ALT-swapped VCF..."

bcftools view -h ${INPUT_VCF} > ${PREFIX}.header.txt

bcftools view -H ${INPUT_VCF} | \
awk 'BEGIN{OFS="\t"}{
    tmp=$4;
    $4=$5;
    $5=tmp;
    print
}' | \
cat ${PREFIX}.header.txt - | \
bgzip -c > ${PREFIX}.swap.vcf.gz

tabix -p vcf ${PREFIX}.swap.vcf.gz

###############################################################################
# Step 3. Liftover swapped VCF
###############################################################################

echo "[INFO] Running swapped liftover..."

CrossMap vcf \
    ${CHAIN} \
    ${PREFIX}.swap.vcf.gz \
    ${REFERENCE} \
    ${PREFIX}.swap.lifted.hg38.vcf

bcftools view ${PREFIX}.swap.lifted.hg38.vcf -Oz -o ${PREFIX}.swap.lifted.hg38.vcf.gz -W=tbi

###############################################################################
# Step 4. Truvari benchmarking (direct liftover)
###############################################################################

echo "[INFO] Benchmarking direct liftover against gnomAD-SV..."

truvari bench \
    --base ${GNOMAD_VCF} \
    --comp ${PREFIX}.lifted.hg38.vcf.gz \
    --pctovl 0.5 \
    --pctsize 0.7 \
    --refdist 1000 \
    -o truvari_out_hg38

###############################################################################
# Step 5. Truvari benchmarking (swap liftover)
###############################################################################

echo "[INFO] Benchmarking swap-liftover against gnomAD-SV..."

truvari bench \
    --base ${GNOMAD_VCF} \
    --comp ${PREFIX}.swap.lifted.hg38.vcf.gz \
    --pctovl 0.5 \
    --pctsize 0.7 \
    --refdist 1000 \
    -o truvari_out_hg38_swap

###############################################################################
# Done
###############################################################################

echo "[INFO] Pipeline completed successfully."
