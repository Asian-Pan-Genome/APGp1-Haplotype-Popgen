#!/usr/bin/env bash

set -euo pipefail

########################################
# Prepare reference genome
#
# Usage:
#   bash 00.prepare_reference.sh config/hg38.config
########################################

CONFIG=$1
source "${CONFIG}"

mkdir -p logs

echo "[INFO] Preparing reference: ${REFERENCE}"
date

# samtools index
if [[ ! -f "${REFERENCE}.fai" ]]; then
    samtools faidx "${REFERENCE}"
fi

# bwa index
if [[ ! -f "${REFERENCE}.bwt" ]]; then
    bwa index "${REFERENCE}"
fi

# gatk dictionary
DICT="${REFERENCE%.*}.dict"

if [[ ! -f "${DICT}" ]]; then
    ${GATK} CreateSequenceDictionary \
        -R "${REFERENCE}" \
        -O "${DICT}"
fi

date
echo "[INFO] Reference preparation completed."
