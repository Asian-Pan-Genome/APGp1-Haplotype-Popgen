#!/usr/bin/env bash

set -euo pipefail

########################################
# BAM postprocessing
#
# Includes:
#   1. Duplicate removal
#   2. Coordinate sorting
#   3. BAM indexing
#
# Usage:
#   bash 02.postprocess_bam.sh \
#       config/hg38.config \
#       SAMPLE
########################################

CONFIG=$1
SAMPLE=$2

source "${CONFIG}"

OUTDIR="bam/${SAMPLE}"
TMPDIR="${OUTDIR}/tmp"

mkdir -p "${TMPDIR}" logs

INPUT_BAM="${OUTDIR}/${SAMPLE}.filtered.bam"
RMdup_BAM="${OUTDIR}/${SAMPLE}.rmdup.bam"
SORT_BAM="${OUTDIR}/${SAMPLE}.rmdup.sorted.bam"

date
echo "[INFO] Processing BAM: ${SAMPLE}"

# duplicate removal
sambamba markdup \
    -r \
    -p \
    -t 10 \
    --tmpdir="${TMPDIR}" \
    "${INPUT_BAM}" \
    "${RMdup_BAM}"

# coordinate sorting
samtools sort \
    -@ 10 \
    -o "${SORT_BAM}" \
    "${RMdup_BAM}"

# BAM index
samtools index \
    -@ 10 \
    "${SORT_BAM}"

date
echo "[INFO] BAM postprocessing completed."
