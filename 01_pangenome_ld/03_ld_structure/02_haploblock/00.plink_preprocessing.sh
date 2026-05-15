#!/usr/bin/env bash

###############################################################################
# PLINK preprocessing pipeline for haplotype block analysis
#
# Usage:
#   bash 00.plink_preprocessing.sh <chromosome> <variant_type>
#
# Example:
#   bash 00.plink_preprocessing.sh chr1 SNPs
#
# Description:
#   This script performs:
#       1. VCF → PLINK BED conversion
#       2. Sample filtering
#       3. Variant/sample quality control
#       4. LD pruning
#
# Requirements:
#   - PLINK v1.9 or newer
#   - bgzip-compressed VCF (*.vcf.gz)
#
# Inputs:
#   - CHM13-APGp1 pangenome VCF
#   - Sample keep list
#
# Output:
#   block/<chr>/
#       ├── *.bed / *.bim / *.fam
#       ├── *.qc_passed.*
#       ├── *.prune.*
#       └── *.pruned.*
###############################################################################

set -euo pipefail

# ------------------------------- Arguments ---------------------------------- #

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 <chromosome> <variant_type>"
    exit 1
fi

CHR=$1
VAR_TYPE=$2

# ------------------------------ Configuration ------------------------------- #

VCF="${HOME}/1.haplo_block/00.call_variants/03.pangenome/04.CHM13-APGp1/4.flag_paper/9.vcfwave_merge/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.dip.${VAR_TYPE}.vcf.gz"

KEEP_LIST="155.unrelated.APGp1.sample.list"

OUTDIR="block/${CHR}"

THREADS=4

# QC thresholds
MAF=0.01
GENO=0.1
MIND=0.1
HWE=1e-6

mkdir -p "${OUTDIR}"

PREFIX=$(basename "${VCF}" .vcf.gz)

# --------------------------------- Logging ---------------------------------- #

echo "============================================================"
echo "PLINK preprocessing pipeline"
echo "Chromosome : ${CHR}"
echo "Variant type : ${VAR_TYPE}"
echo "Start time : $(date)"
echo "============================================================"

# ----------------------------- Input checking ------------------------------- #

[[ -f "${VCF}" ]] || {
    echo "[ERROR] VCF file not found: ${VCF}"
    exit 1
}

[[ -f "${KEEP_LIST}" ]] || {
    echo "[ERROR] Sample keep list not found: ${KEEP_LIST}"
    exit 1
}

# --------------------------- VCF → BED conversion --------------------------- #

echo
echo "[Step 1] Converting VCF to PLINK BED format..."

plink \
    --vcf "${VCF}" \
    --make-bed \
    --vcf-half-call m \
    --keep "${KEEP_LIST}" \
    --double-id \
    --threads "${THREADS}" \
    --out "${OUTDIR}/${PREFIX}"

# ----------------------------------- QC ------------------------------------- #

echo
echo "[Step 2] Performing variant/sample QC..."

plink \
    --bfile "${OUTDIR}/${PREFIX}" \
    --maf "${MAF}" \
    --geno "${GENO}" \
    --mind "${MIND}" \
    --hwe "${HWE}" \
    --make-bed \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --out "${OUTDIR}/${PREFIX}.qc_passed"

# ------------------------------- LD pruning --------------------------------- #

echo
echo "[Step 3] LD pruning..."

plink \
    --bfile "${OUTDIR}/${PREFIX}.qc_passed" \
    --indep-pairwise 50 5 0.2 \
    --threads "${THREADS}" \
    --out "${OUTDIR}/${PREFIX}.qc_passed.prune"

plink \
    --bfile "${OUTDIR}/${PREFIX}.qc_passed" \
    --extract "${OUTDIR}/${PREFIX}.qc_passed.prune.prune.in" \
    --make-bed \
    --threads "${THREADS}" \
    --out "${OUTDIR}/${PREFIX}.qc_passed.pruned"

# -------------------------------- Completion -------------------------------- #

echo
echo "============================================================"
echo "PLINK preprocessing completed successfully"
echo "Output directory : ${OUTDIR}"
echo "Finish time      : $(date)"
echo "============================================================"
