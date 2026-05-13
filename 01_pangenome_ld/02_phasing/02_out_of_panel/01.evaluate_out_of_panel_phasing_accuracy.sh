#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Evaluate out-of-panel phasing accuracy using WhatsHap
#
# This script:
#   1. Compares assembly-derived phased variants against
#      out-of-panel statistically phased variants
#   2. Computes switch error statistics using WhatsHap
#   3. Outputs switch-error BED files and phasing summaries
#
# Requirements:
#   - whatshap
#
# Usage:
#   bash 01.evaluate_out_of_panel_phasing_accuracy.sh \
#       chr1 \
#       SAMPLE_ID
#
###############################################################################

CHR=$1
SAMPLE=$2

###############################################################################
# Input configuration
###############################################################################

INPUT_DIR=./Assembly

ASSEMBLY_VCF=${INPUT_DIR}/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.SNPs_INDELs_SVs.filtered.norm.rename.bcf

OUT_OF_PANEL_VCF=${INPUT_DIR}/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.SNPs_INDELs_SVs.filtered.norm.rename.1kg.out_of_panel.phased.full.bcf

OUTDIR=${CHR}

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Sample: ${SAMPLE}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Evaluate phasing accuracy
###############################################################################

echo "[INFO] Running WhatsHap compare..."

whatshap compare \
    --sample ${SAMPLE} \
    --switch-error-bed ${OUTDIR}/${SAMPLE}.assem.1kgp.out_of_panel.ser.bed \
    --longest-block-tsv ${OUTDIR}/${SAMPLE}.assem.1kgp.out_of_panel.longest-block.tsv \
    --tsv-pairwise ${OUTDIR}/${SAMPLE}.assem.1kgp.out_of_panel.pairwise.tsv \
    ${ASSEMBLY_VCF} \
    ${OUT_OF_PANEL_VCF}

###############################################################################
# Done
###############################################################################

echo "[INFO] Out-of-panel phasing evaluation completed successfully."

date
