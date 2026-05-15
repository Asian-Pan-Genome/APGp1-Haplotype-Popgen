#!/usr/bin/env bash

set -euo pipefail

########################################
# Joint genotyping
#
# Usage:
#   bash 05.genotype_gvcf.sh \
#       config/hg38.config
########################################

CONFIG=$1

source "${CONFIG}"

mkdir -p joint_calling logs

INPUT="joint_calling/cohort.g.vcf.gz"
OUTPUT="joint_calling/cohort.raw.vcf.gz"

date
echo "[INFO] Running joint genotyping"

${GATK} GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V "${INPUT}" \
    -O "${OUTPUT}"

tabix -p vcf "${OUTPUT}"

date
echo "[INFO] Joint genotyping completed."
