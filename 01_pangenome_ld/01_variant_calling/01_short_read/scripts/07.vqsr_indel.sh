#!/usr/bin/env bash

set -euo pipefail

########################################
# INDEL VQSR
#
# Usage:
#   bash 07.vqsr_indel.sh \
#       config/hg38.config
########################################

CONFIG=$1

source "${CONFIG}"

mkdir -p vqsr logs tmp

INPUT="vqsr/cohort.snp.vqsr.vcf.gz"

date
echo "[INFO] Running INDEL VQSR"

########################################
# VariantRecalibrator
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=tmp" \
    VariantRecalibrator \
    -R "${REFERENCE}" \
    -V "${INPUT}" \
    -mode INDEL \
    -O vqsr/cohort.indel.recal \
    --tranches-file vqsr/cohort.indel.tranches \
    --resource:mills,known=true,training=true,truth=true,prior=12.0 "${MILLS}" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${KNOWN_INDELS}" \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an ReadPosRankSum \
    -an MQRankSum \
    --max-gaussians 2

########################################
# ApplyVQSR
########################################

${GATK} \
    --java-options "-Xmx50G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=tmp" \
    ApplyVQSR \
    -R "${REFERENCE}" \
    -V "${INPUT}" \
    -mode INDEL \
    --recal-file vqsr/cohort.indel.recal \
    --tranches-file vqsr/cohort.indel.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -O vqsr/cohort.final.vcf.gz

tabix -p vcf vqsr/cohort.final.vcf.gz

date
echo "[INFO] INDEL VQSR completed."
