#!/usr/bin/env bash

set -euo pipefail

########################################
# BQSR + HaplotypeCaller
#
# Usage:
#   bash 03.bqsr_haplotypecaller.sh \
#       config/hg38.config \
#       SAMPLE
########################################

CONFIG=$1
SAMPLE=$2

source "${CONFIG}"

BAMDIR="bam/${SAMPLE}"
GVCFDIR="gvcf/${SAMPLE}"
TMPDIR="tmp/${SAMPLE}"

mkdir -p "${GVCFDIR}" "${TMPDIR}" logs

INPUT_BAM="${BAMDIR}/${SAMPLE}.rmdup.sorted.bam"
BQSR_TABLE="${BAMDIR}/${SAMPLE}.bqsr.table"
BQSR_BAM="${BAMDIR}/${SAMPLE}.bqsr.bam"
OUTPUT_GVCF="${GVCFDIR}/${SAMPLE}.g.vcf.gz"

date
echo "[INFO] Running BQSR: ${SAMPLE}"

########################################
# BaseRecalibrator
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${TMPDIR}" \
    BaseRecalibrator \
    -R "${REFERENCE}" \
    -I "${INPUT_BAM}" \
    --known-sites "${KNOWN_SITES}" \
    -O "${BQSR_TABLE}"

########################################
# ApplyBQSR
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${TMPDIR}" \
    ApplyBQSR \
    -R "${REFERENCE}" \
    -I "${INPUT_BAM}" \
    --bqsr-recal-file "${BQSR_TABLE}" \
    -O "${BQSR_BAM}"

########################################
# HaplotypeCaller
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${TMPDIR}" \
    HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${BQSR_BAM}" \
    -O "${OUTPUT_GVCF}" \
    --emit-ref-confidence GVCF \
    -stand-call-conf 30

tabix -p vcf "${OUTPUT_GVCF}"

date
echo "[INFO] HaplotypeCaller completed."
