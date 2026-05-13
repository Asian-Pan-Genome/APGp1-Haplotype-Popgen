#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run iHH12 selection scans within a population
#
# This script:
#   1. Loads population-specific phased VCFs
#   2. Uses predicted genetic maps for recombination-aware scans
#   3. Runs iHH12 scans using selscan
#   4. Supports scans on MAF-filtered datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 01.run_ihh12_selection_scan.sh \
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

INPUT_DIR="../10.nSL/vcf/${CHR}"

VCF_MAF1=${INPUT_DIR}/${CHR}.${POP}.${MASK}.maf1pct.vcf.gz
VCF_MAF5=${INPUT_DIR}/${CHR}.${POP}.${MASK}.maf5pct.vcf.gz

GENETIC_MAP=./recom_map/${CHR}.${MASK}.predicted.map.txt

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
# Step 1. iHH12 scan (MAF >= 5%)
###############################################################################

echo "[INFO] Running iHH12 scan (MAF >= 5%)..."

selscan \
    --ihh12 \
    --vcf ${VCF_MAF5} \
    --map ${GENETIC_MAP} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf5pct \
    --threads ${THREADS}

###############################################################################
# Optional: iHH12 scan (MAF >= 1%)
###############################################################################
#
# Uncomment the following block if MAF >= 1% scans are required.
#
# selscan \
#     --ihh12 \
#     --vcf ${VCF_MAF1} \
#     --map ${GENETIC_MAP} \
#     --out ${OUTDIR}/${CHR}.${MASK}.maf1pct \
#     --threads ${THREADS}
#
###############################################################################

###############################################################################
# Done
###############################################################################

echo "[INFO] iHH12 scans completed successfully."

date
