#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Normalize genome-wide XP-nSL scores
#
# This script:
#   1. Collects chromosome-level XP-nSL output files
#   2. Performs genome-wide normalization using selscan norm
#   3. Normalizes both MAF >= 1% and MAF >= 5% scans
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 02.normalize_xpnsl_scores.sh \
#       mask_name \
#       focal_population \
#       outgroup_population
#
###############################################################################

MASK=$1
POP1=$2
POP2=$3

###############################################################################
# Configuration
###############################################################################

COMPARISON_DIR=${POP1}_vs_${POP2}

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Focal population: ${POP1}"
echo "[INFO] Outgroup population: ${POP2}"
echo "[INFO] Mask: ${MASK}"

###############################################################################
# Step 1. Normalize XP-nSL scores (MAF >= 1%)
###############################################################################

echo "[INFO] Normalizing XP-nSL scores (MAF >= 1%)..."

MAF1_FILES=$(for chr in chr{1..22}; do
    echo "./${COMPARISON_DIR}/${chr}/${chr}.${MASK}.maf1pct.xpnsl.out"
done)

norm \
    --xpnsl \
    --files ${MAF1_FILES} \
    --bp-win

###############################################################################
# Step 2. Normalize XP-nSL scores (MAF >= 5%)
###############################################################################

echo "[INFO] Normalizing XP-nSL scores (MAF >= 5%)..."

MAF5_FILES=$(for chr in chr{1..22}; do
    echo "./${COMPARISON_DIR}/${chr}/${chr}.${MASK}.maf5pct.xpnsl.out"
done)

norm \
    --xpnsl \
    --files ${MAF5_FILES} \
    --bp-win

###############################################################################
# Done
###############################################################################

echo "[INFO] XP-nSL normalization completed successfully."

date
