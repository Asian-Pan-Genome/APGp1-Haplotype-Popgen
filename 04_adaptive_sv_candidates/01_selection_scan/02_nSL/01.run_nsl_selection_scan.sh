#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run nSL selection scans within a population
#
# This script:
#   1. Loads population-specific phased VCFs
#   2. Runs nSL scans using selscan
#   3. Performs scans for both MAF >= 1% and MAF >= 5% datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 01.run_nsl_selection_scan.sh \
#       chr1 \
#       mask_name \
#       population
#
###############################################################################

CHR=$1
MASK=$2
POP=$3

###############################################################################
# Input configuration
###############################################################################

INPUT_DIR=vcf/${CHR}

VCF_MAF1=${INPUT_DIR}/${CHR}.${POP}.${MASK}.maf1pct.vcf.gz
VCF_MAF5=${INPUT_DIR}/${CHR}.${POP}.${MASK}.maf5pct.vcf.gz

OUTDIR=${POP}/${CHR}

THREADS=8

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Population: ${POP}"
echo "[INFO] Mask: ${MASK}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. nSL scan (MAF >= 1%)
###############################################################################

echo "[INFO] Running nSL scan (MAF >= 1%)..."

selscan \
    --nsl \
    --vcf ${VCF_MAF1} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf1pct \
    --threads ${THREADS}

###############################################################################
# Step 2. nSL scan (MAF >= 5%)
###############################################################################

echo "[INFO] Running nSL scan (MAF >= 5%)..."

selscan \
    --nsl \
    --vcf ${VCF_MAF5} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf5pct \
    --threads ${THREADS}

###############################################################################
# Done
###############################################################################

echo "[INFO] nSL scans completed successfully."

date
