#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run XP-EHH selection scans between focal and outgroup populations
#
# This script:
#   1. Loads population-specific phased VCFs
#   2. Uses predicted genetic maps for recombination-aware scans
#   3. Runs XP-EHH scans using selscan
#   4. Performs scans on MAF-filtered datasets
#
# Requirements:
#   - selscan
#
# Usage:
#   bash 01.run_xpehh_selection_scan.sh \
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

INPUT_DIR="../04_XP-nSL/vcf/${CHR}"

FOCAL_VCF_MAF1=${INPUT_DIR}/${CHR}.${FOCAL_POP}.${MASK}.maf1pct.vcf.gz
FOCAL_VCF_MAF5=${INPUT_DIR}/${CHR}.${FOCAL_POP}.${MASK}.maf5pct.vcf.gz

OUTGROUP_VCF_MAF1=${INPUT_DIR}/${CHR}.${OUTGROUP_POP}.${MASK}.maf1pct.vcf.gz
OUTGROUP_VCF_MAF5=${INPUT_DIR}/${CHR}.${OUTGROUP_POP}.${MASK}.maf5pct.vcf.gz

GENETIC_MAP=./recom_map/${CHR}.${MASK}.predicted.map.txt

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
# Step 1. XP-EHH scan (MAF >= 5%)
###############################################################################

echo "[INFO] Running XP-EHH scan (MAF >= 5%)..."

selscan \
    --xpehh \
    --vcf ${FOCAL_VCF_MAF5} \
    --vcf-ref ${OUTGROUP_VCF_MAF5} \
    --map ${GENETIC_MAP} \
    --out ${OUTDIR}/${CHR}.${MASK}.maf5pct \
    --threads ${THREADS}

###############################################################################
# Optional: XP-EHH scan (MAF >= 1%)
###############################################################################
#
# Uncomment the following block if MAF >= 1% scans are required.
#
# selscan \
#     --xpehh \
#     --vcf ${FOCAL_VCF_MAF1} \
#     --vcf-ref ${OUTGROUP_VCF_MAF1} \
#     --map ${GENETIC_MAP} \
#     --out ${OUTDIR}/${CHR}.${MASK}.maf1pct \
#     --threads ${THREADS}
#
###############################################################################

###############################################################################
# Done
###############################################################################

echo "[INFO] XP-EHH scans completed successfully."

date
