#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Normalize genome-wide iHH12 scores
#
# This script:
#   1. Collects chromosome-level iHH12 scan results
#   2. Performs genome-wide normalization using selscan norm
#   3. Supports normalization of MAF-filtered datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 02.normalize_ihh12_scores.sh \
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
# Step 1. Normalize iHH12 scores (MAF >= 5%)
###############################################################################

echo "[INFO] Normalizing iHH12 scores (MAF >= 5%)..."

MAF5_FILES=$(for chr in chr{1..22}; do
    echo "./${SCAN_DIR}/${chr}/${chr}.${MASK}.maf5pct.ihh12.out"
done)

norm \
    --ihh12 \
    --files ${MAF5_FILES} \
    --bp-win

###############################################################################
# Optional: Normalize iHH12 scores (MAF >= 1%)
###############################################################################
#
# Uncomment the following block if MAF >= 1% normalization is required.
#
# MAF1_FILES=$(for chr in chr{1..22}; do
#     echo "./${SCAN_DIR}/${chr}/${chr}.${MASK}.maf1pct.ihh12.out"
# done)
#
# norm \
#     --ihh12 \
#     --files ${MAF1_FILES} \
#     --bp-win
#
###############################################################################

###############################################################################
# Done
###############################################################################

echo "[INFO] iHH12 normalization completed successfully."

date
