#!/usr/bin/env bash

set -euo pipefail

########################################
# SNP VQSR
#
# Usage:
#   bash 06.vqsr_snp.sh \
#       config/hg38.config
########################################

CONFIG=$1

source "${CONFIG}"

mkdir -p vqsr logs tmp

INPUT="joint_calling/cohort.raw.vcf.gz"

date
echo "[INFO] Running SNP VQSR"

########################################
# VariantRecalibrator
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=tmp" \
    VariantRecalibrator \
    -R "${REFERENCE}" \
    -V "${INPUT}" \
    -mode SNP \
    -O vqsr/cohort.snp.recal \
    --tranches-file vqsr/cohort.snp.tranches \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${HAPMAP}" \
    --resource:omini,known=false,training=true,truth=false,prior=12.0 "${OMNI}" \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${PHASE1}" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${DBSNP}" \
    -an QD \
    -an FS \
    -an MQ \
    -an ReadPosRankSum \
    -an MQRankSum \
    -an SOR \
    -an DP \
    --max-gaussians 2

########################################
# ApplyVQSR
########################################

${GATK} \
    --java-options "-Xmx100G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=tmp" \
    ApplyVQSR \
    -R "${REFERENCE}" \
    -V "${INPUT}" \
    -mode SNP \
    --recal-file vqsr/cohort.snp.recal \
    --tranches-file vqsr/cohort.snp.tranches \
    --truth-sensitivity-filter-level 99.8 \
    -O vqsr/cohort.snp.vqsr.vcf.gz

tabix -p vcf vqsr/cohort.snp.vqsr.vcf.gz

date
echo "[INFO] SNP VQSR completed."
