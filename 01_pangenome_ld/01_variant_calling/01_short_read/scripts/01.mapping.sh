#!/usr/bin/env bash

set -euo pipefail

########################################
# Read mapping
#
# Usage:
#   bash 01.mapping.sh \
#       config/hg38.config \
#       SAMPLE \
#       sample_R1.fastq.gz \
#       sample_R2.fastq.gz
########################################

CONFIG=$1
SAMPLE=$2
R1=$3
R2=$4

source "${CONFIG}"

OUTDIR="bam/${SAMPLE}"
TMPDIR="${OUTDIR}/tmp"

mkdir -p "${OUTDIR}" "${TMPDIR}" logs

date
echo "[INFO] Mapping sample: ${SAMPLE}"

bwa mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:WGS\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" | \
samtools view \
    -@ ${THREADS} \
    -bh \
    -q 30 \
    -F 256 \
    -F 2048 \
    -o "${OUTDIR}/${SAMPLE}.filtered.bam"

date
echo "[INFO] Mapping completed."
