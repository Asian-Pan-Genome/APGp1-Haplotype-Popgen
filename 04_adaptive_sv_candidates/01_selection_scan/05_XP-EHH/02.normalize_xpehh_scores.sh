#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Normalize genome-wide XP-EHH scores
#
# This script:
#   1. Collects chromosome-level XP-EHH scan results
#   2. Performs genome-wide normalization using selscan norm
#   3. Supports normalization of MAF-filtered datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 02.normalize_xpehh_scores.sh \
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
# Step 1. Normalize XP-EHH scores (MAF >= 5%)
###############################################################################

echo "[INFO] Normalizing XP-EHH scores (MAF >= 5%)..."

MAF5_FILES=$(for chr in chr{1..22}; do
    echo "./${COMPARISON_DIR}/${chr}/${chr}.${MASK}.maf5pct.xpehh.out"
done)

norm \
    --xpehh \
    --files ${MAF5_FILES} \
    --bp-win

###############################################################################
# Optional: Normalize XP-EHH scores (MAF >= 1%)
###############################################################################
#
# Uncomment the following block if MAF >= 1% normalization is required.
#
# MAF1_FILES=$(for chr in chr{1..22}; do
#     echo "./${COMPARISON_DIR}/${chr}/${chr}.${MASK}.maf1pct.xpehh.out"
# done)
#
# norm \
#     --xpehh \
#     --files ${MAF1_FILES} \
#     --bp-win
#
###############################################################################

###############################################################################
# Done
###############################################################################

echo "[INFO] XP-EHH normalization completed successfully."

date
