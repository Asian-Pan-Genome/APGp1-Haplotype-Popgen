#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run XP-nSL selection scans between focal and outgroup populations
#
# This script:
#   1. Loads population-specific MAF-filtered VCFs
#   2. Runs XP-nSL scans using selscan
#   3. Performs scans for both MAF >= 1% and MAF >= 5% datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 01.run_xpnsl_selection_scan.sh \
#       chr1 \
#       mask_name \
#       focal_population \
#       outgroup_population
#
###############################################################################

CHR=$1
MASK=$2
FOCAL_POP=$3
OUTGROUP_POP=$4

###############################################################################
# Input configuration
###############################################################################

INPUT_DIR=vcf/${CHR}

FOCAL_VCF_MAF1=${INPUT_DIR}/${CHR}.${FOCAL_POP}.${MASK}.maf1pct.vcf.gz
FOCAL_VCF_MAF5=${INPUT_DIR}/${CHR}.${FOCAL_POP}.${MASK}.maf5pct.vcf.gz

OUTGROUP_VCF_MAF1=${INPUT_DIR}/${CHR}.${OUTGROUP_POP}.${MASK}.maf1pct.vcf.gz
OUTGROUP_VCF_MAF5=${INPUT_DIR}/${CHR}.${OUTGROUP_POP}.${MASK}.maf5pct.vcf.gz

OUTDIR=${FOCAL_POP}_vs_${OUTGROUP_POP}/${CHR}

THREADS=8

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Focal population: ${FOCAL_POP}"
echo "[INFO] Outgroup population: ${OUTGROUP_POP}"
echo "[INFO] Mask: ${MASK}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. XP-nSL scan (MAF >= 1%)
###############################################################################

echo "[INFO] Running XP-nSL scan (MAF >= 1%)..."

selscan \
    --xpnsl \
    --vcf ${FOCAL_VCF_MAF1} \
    --vcf-ref ${OUTGROUP_VCF_MAF1} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf1pct \
    --threads ${THREADS}

###############################################################################
# Step 2. XP-nSL scan (MAF >= 5%)
###############################################################################

echo "[INFO] Running XP-nSL scan (MAF >= 5%)..."

selscan \
    --xpnsl \
    --vcf ${FOCAL_VCF_MAF5} \
    --vcf-ref ${OUTGROUP_VCF_MAF5} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf5pct \
    --threads ${THREADS}

###############################################################################
# Done
###############################################################################

echo "[INFO] XP-nSL scans completed successfully."

date
