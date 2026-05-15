#!/usr/bin/env bash

###############################################################################
# Pangenome variant preprocessing and filtering pipeline
#
# Description:
#   This script preprocesses Minigraph-Cactus VCF files for downstream
#   population genetic analyses, including:
#     1. Sample subsetting
#     2. Diploid genotype conversion
#     3. Reference normalization and sorting
#     4. Biallelic variant extraction
#     5. Variant annotation and filtering
#
# Usage:
#   bash 01.variant_filter.sh <var_type> <chr>
#
# Example:
#   bash 01.variant_filter.sh [SNP/Indel/MNP/SV] chr1
#
# Requirements:
#   - bcftools >= 1.16
#   - python3
#   - hap2dip.py
#
# Input:
#   CHM13-APGp1-HPRCp1-HGSVCp3_MC.remove_telo_cent.final.<var_type>.vcf.gz
#   CHM13-APGp1_MC.remove_telo_cent.final.<var_type>.vcf.gz
#
# Output:
#   Final filtered VCF:
#   *.biallelic.addtag.filter.vcf.gz
#
###############################################################################

set -euo pipefail

############################
# Arguments
############################

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 <var_type> <chr>"
    exit 1
fi

VAR_TYPE="$1"
CHR="$2"

############################
# Configuration
############################

THREADS=16

INPUT_DIR="pangenome/CHM13-APGp1-HPRCp1-HGSVCp3"
OUT_DIR="536_hap/${CHR}"

REF="chm13v2.0.fa"

VCF="CHM13-APGp1-HPRCp1-HGSVCp3_MC.remove_telo_cent.final.${VAR_TYPE}.vcf.gz"

SAMPLE_LIST="536_hap.txt"
FAM_FILE="268.sample.dip.fam"

MISSING_SCRIPT="./hap2dip.py"

mkdir -p "${OUT_DIR}"

PREFIX="${OUT_DIR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${VAR_TYPE}"

echo "[$(date)] Starting SV preprocessing pipeline..."
echo "Chromosome : ${CHR}"
echo "Variant type : ${VAR_TYPE}"
echo

###############################################################################
# Step 1. Sample filtering
###############################################################################

echo "[$(date)] Step 1: Sample filtering"

bcftools view \
    -S "${SAMPLE_LIST}" \
    -r "${CHR}" \
    --min-ac 1 \
    "${INPUT_DIR}/${VCF}" \
    -Oz \
    -o "${PREFIX}.vcf.gz" \
    --threads "${THREADS}" \
    -W=tbi

echo "[$(date)] Step 1 completed"
echo

###############################################################################
# Step 2. Convert to diploid format
###############################################################################

echo "[$(date)] Step 2: Diploid genotype conversion"

bcftools view \
    "${PREFIX}.vcf.gz" \
    -o "${PREFIX}.vcf"

python3 "${MISSING_SCRIPT}" \
    -dip \
    -fam "${FAM_FILE}" \
    -vcf "${PREFIX}.vcf" \
    | bcftools view \
        -Oz \
        -o "${PREFIX}.dip.vcf.gz" \
        -W=tbi

rm -f "${PREFIX}.vcf"

echo "[$(date)] Step 2 completed"
echo

###############################################################################
# Step 3. Normalize, sort, and clean INFO fields
###############################################################################

echo "[$(date)] Step 3: Normalization and sorting"

bcftools norm \
    --check-ref w \
    -f "${REF}" \
    "${PREFIX}.dip.vcf.gz" \
    --threads "${THREADS}" \
    | bcftools sort - \
    | bcftools annotate \
        -x INFO/AF \
        --threads "${THREADS}" \
    | bcftools view \
        --trim-alt-alleles \
        -Oz \
        -o "${PREFIX}.filtered.vcf.gz" \
        --threads "${THREADS}" \
        -W=tbi

echo "[$(date)] Step 3 completed"
echo

###############################################################################
# Step 4. Extract biallelic variants
###############################################################################

echo "[$(date)] Step 4: Biallelic filtering"

bcftools view \
    -m2 -M2 \
    "${PREFIX}.filtered.vcf.gz" \
    -Oz \
    -o "${PREFIX}.biallelic.vcf.gz" \
    --threads "${THREADS}" \
    -W=tbi

echo "[$(date)] Step 4 completed"
echo

###############################################################################
# Step 5. Recalculate INFO tags
###############################################################################

echo "[$(date)] Step 5: INFO tag annotation"

bcftools annotate \
    -x INFO/AN,INFO/AC \
    "${PREFIX}.biallelic.vcf.gz" \
    --threads "${THREADS}" \
    | bcftools +fill-tags \
        --threads "${THREADS}" \
        -- \
        -t AN,AC,AF,MAF,F_MISSING,MAC:1=MAC,HWE \
    | bcftools view \
        -Oz \
        -o "${PREFIX}.biallelic.addtag.vcf.gz" \
        --threads "${THREADS}" \
        -W=tbi

echo "[$(date)] Step 5 completed"
echo

###############################################################################
# Step 6. Final variant filtering
###############################################################################

echo "[$(date)] Step 6: Variant filtering"

# Filtering criteria:
#   - Remove symbolic alleles (*)
#   - Missing rate > 5%
#   - Monomorphic variants
#   - Hardy-Weinberg equilibrium P < 1e-10 (for non-PAR in chrX and chrY, do not use this flag)

bcftools view \
    -e "ALT=='*' || F_MISSING > 0.05 || INFO/MAC==0 || INFO/HWE < 1e-10" \
    "${PREFIX}.biallelic.addtag.vcf.gz" \
    -Oz \
    -o "${PREFIX}.biallelic.addtag.filter.vcf.gz" \
    --threads "${THREADS}" \
    -W=tbi

echo "[$(date)] Step 6 completed"
echo

###############################################################################
# Done
###############################################################################

echo "[$(date)] All steps completed successfully."
echo
echo "Final output:"
echo "  ${PREFIX}.biallelic.addtag.filter.vcf.gz"
echo
echo "Pipeline finished."
