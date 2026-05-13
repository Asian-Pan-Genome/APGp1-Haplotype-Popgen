#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Filter population variants and generate MAF-filtered VCFs
#
# This script:
#   1. Subsets variants by sample list
#   2. Calculates allele frequency statistics
#   3. Generates MAF-filtered VCF files
#
# Requirements:
#   - bcftools
#
# Usage:
#   bash 00.prep_input_vcf.sh \
#       chr1 \
#       mask_name \
#       sample_list.txt
#
###############################################################################

CHR=$1
MASK=$2
SAMPLE_LIST=$3

###############################################################################
# Input configuration
###############################################################################

INPUT_DIR="vcfwave_merge/4.polarize/${CHR}"

VCF=CHM13-APGp1-HPRCp1-HGSVCp3_MC.ex_tc.final.SNPs_SVs.biallelic.addtag.filter.unrelated.${MASK}.vcf.gz

OUTDIR=vcf/${CHR}

PREFIX=${CHR}.$(basename ${SAMPLE_LIST} .txt).${MASK}

THREADS=16

###############################################################################
# Start
###############################################################################

date
echo "[INFO] Processing chromosome: ${CHR}"
echo "[INFO] Sample list: ${SAMPLE_LIST}"
echo "[INFO] Mask: ${MASK}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Subset samples and calculate allele-frequency tags
###############################################################################

echo "[INFO] Generating sample-filtered VCF..."

bcftools +fill-tags \
    --threads ${THREADS} \
    ${INPUT_DIR}/${VCF} -- \
    -t AN,AC,AF,MAF,F_MISSING,MAC:1=MAC,HWE | \
bcftools view \
    -S ${SAMPLE_LIST} \
    -Oz \
    --threads ${THREADS} \
    -o ${OUTDIR}/${PREFIX}.vcf.gz \
    -W=tbi

###############################################################################
# Step 2. Generate MAF >= 1% VCF
###############################################################################

echo "[INFO] Filtering variants with MAF >= 1%..."

bcftools view \
    -i 'MAF>=0.01' \
    ${OUTDIR}/${PREFIX}.vcf.gz \
    -Oz \
    --threads ${THREADS} \
    -o ${OUTDIR}/${PREFIX}.maf1pct.vcf.gz \
    -W=tbi

###############################################################################
# Step 3. Generate MAF >= 5% VCF
###############################################################################

echo "[INFO] Filtering variants with MAF >= 5%..."

bcftools view \
    -i 'MAF>=0.05' \
    ${OUTDIR}/${PREFIX}.vcf.gz \
    -Oz \
    --threads ${THREADS} \
    -o ${OUTDIR}/${PREFIX}.maf5pct.vcf.gz \
    -W=tbi

###############################################################################
# Done
###############################################################################

echo "[INFO] Variant filtering completed successfully."

date
