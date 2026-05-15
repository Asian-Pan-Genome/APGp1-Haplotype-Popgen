#!/usr/bin/env bash

set -euo pipefail

########################################
# Combine GVCFs
#
# Usage:
#   bash 04.combine_gvcf.sh \
#       config/hg38.config \
#       gvcf.list
########################################

CONFIG=$1
GVCF_LIST=$2

source "${CONFIG}"

mkdir -p joint_calling logs

OUTPUT="joint_calling/cohort.g.vcf.gz"

date
echo "[INFO] Combining GVCFs"

${GATK} CombineGVCFs \
    -R "${REFERENCE}" \
    $(awk '{print "-V "$1}' "${GVCF_LIST}") \
    -O "${OUTPUT}"

date
echo "[INFO] CombineGVCFs completed."
