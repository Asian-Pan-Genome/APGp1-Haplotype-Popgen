#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Predict genetic map positions for query variants
#
# This script:
#   1. Loads a reference recombination map
#   2. Predicts genetic map positions for query variants
#      using predictGMAP
#
# Requirements:
#   - predictGMAP
#
# Usage:
#   bash 01.predict_genetic_map_positions.sh \
#       chr1 \
#       mask_name
#
###############################################################################

CHR=$1
MASK=$2

###############################################################################
# Input configuration
###############################################################################

REF_MAP="chm13v2.0/1kgp_3202/genetic_map/masked/averaged/${CHR}.no_chr.masked.map.txt"

QUERY_POS=../04_XP-nSL/vcf/${CHR}/${CHR}.CHA.${MASK}.maf5pct.pos.txt

OUTDIR=./recom_map

PREDICT_GMAP=predictGMAP

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Mask: ${MASK}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Predict genetic map positions
###############################################################################

echo "[INFO] Predicting genetic map positions..."

${PREDICT_GMAP} \
    --out ${OUTDIR}/${CHR}.${MASK}.predicted.map.txt \
    --query ${QUERY_POS} \
    --ref ${REF_MAP}

###############################################################################
# Done
###############################################################################

echo "[INFO] Genetic map prediction completed successfully."

date
