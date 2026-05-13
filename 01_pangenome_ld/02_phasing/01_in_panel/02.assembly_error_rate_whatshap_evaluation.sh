#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Evaluate phasing switch errors using WhatsHap
#
# This script:
#   1. Compares phased assembly variants against
#      trio-phased benchmark callsets
#   2. Computes switch and flip error statistics using WhatsHap
#   3. Outputs switch- & flip-error BED files and phasing summaries
#
# Requirements:
#   - whatshap
#
# Usage:
#   bash 02.assembly_error_rate_whatshap_evaluation.sh \
#       chr1 \
#       SAMPLE_ID
#
###############################################################################

CHR=$1
SAMPLE=$2

###############################################################################
# Input configuration
###############################################################################

DATA_DIR=./intersected_site_vcf

BENCHMARK_VCF="trio_ngs/${CHR}/${CHR}.phased.recalibrated.addtag.filter.masked.vcf.gz"

DATA_VCF=${DATA_DIR}/${CHR}/isec/0001.bcf

OUTDIR=${CHR}

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Sample: ${SAMPLE}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Evaluate switch errors
###############################################################################

echo "[INFO] Running WhatsHap compare..."

whatshap compare \
    --sample ${SAMPLE} \
    --switch-error-bed ${OUTDIR}/${SAMPLE}.overlap.assem.1kgp.ser.bed \
    --longest-block-tsv ${OUTDIR}/${SAMPLE}.overlap.assem.1kgp.longest-block.tsv \
    --tsv-pairwise ${OUTDIR}/${SAMPLE}.overlap.assem.1kgp.pairwise.tsv \
    ${BENCHMARK_VCF} \
    ${DATA_VCF}

###############################################################################
# Done
###############################################################################

echo "[INFO] Switch error evaluation completed successfully."

date
