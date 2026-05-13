#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Normalize genome-wide nSL scores
#
# This script:
#   1. Collects chromosome-level nSL scan results
#   2. Performs genome-wide normalization using selscan norm
#   3. Normalizes both MAF >= 1% and MAF >= 5% scans
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 02.normalize_nsl_scores.sh \
#       population \
#       mask_name
#
###############################################################################

POP=$1
MASK=$2

###############################################################################
# Configuration
###############################################################################

SCAN_DIR=${POP}

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Population: ${POP}"
echo "[INFO] Mask: ${MASK}"

###############################################################################
# Step 1. Normalize nSL scores (MAF >= 1%)
###############################################################################

echo "[INFO] Normalizing nSL scores (MAF >= 1%)..."

MAF1_FILES=$(for chr in chr{1..22}; do
    echo "./${SCAN_DIR}/${chr}/${chr}.${MASK}.maf1pct.nsl.out"
done)

norm \
    --nsl \
    --files ${MAF1_FILES} \
    --bp-win

###############################################################################
# Step 2. Normalize nSL scores (MAF >= 5%)
###############################################################################

echo "[INFO] Normalizing nSL scores (MAF >= 5%)..."

MAF5_FILES=$(for chr in chr{1..22}; do
    echo "./${SCAN_DIR}/${chr}/${chr}.${MASK}.maf5pct.nsl.out"
done)

norm \
    --nsl \
    --files ${MAF5_FILES} \
    --bp-win

###############################################################################
# Done
###############################################################################

echo "[INFO] nSL normalization completed successfully."

date
