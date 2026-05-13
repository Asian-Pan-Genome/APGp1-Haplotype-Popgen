#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Benchmark lifted structural variants against the 1KCP SV resource
#
# This script:
#   1. Benchmarks lifted SVs using Truvari
#   2. Extracts true-positive matched SVs
#   3. Generates paired SV coordinates between benchmark sets
#
# Requirements:
#   - truvari
#   - bcftools
#   - awk
#   - sed
#
# Usage:
#   bash 01.benchmark_lifted_svs_against_1kcp.sh \
#       511.candidate.sv.chm13.lifted.hg38.vcf.gz \
#       1kcp.sv.vcf.gz
#
###############################################################################

BASE_VCF=$1
COMP_VCF=$2

OUTDIR="truvari_out_hg38"

###############################################################################
# Step 1. Truvari benchmarking
###############################################################################

echo "[INFO] Running Truvari benchmark..."

truvari bench \
    --base ${BASE_VCF} \
    --comp ${COMP_VCF} \
    --pctovl 0.5 \
    --pctsize 0.7 \
    --refdist 1000 \
    -o ${OUTDIR}

###############################################################################
# Step 2. Extract matched SV records from comparison set
###############################################################################

echo "[INFO] Extracting true-positive SV records..."

bcftools view -H ${OUTDIR}/tp-comp.vcf.gz | \
sed 's/^chr//g' | \
awk 'BEGIN{OFS="\t"}{
    print $3, $1, $2, $4, $5
}' > tp.sv.tsv

###############################################################################
# Step 3. Generate paired SV coordinates
###############################################################################

echo "[INFO] Generating paired SV table..."

paste \
    <(bcftools view -H ${OUTDIR}/tp-base.vcf.gz | cut -f1-3) \
    <(bcftools view -H ${OUTDIR}/tp-comp.vcf.gz | cut -f1-3) | \
sort -k1,2V > tp.paired.sv.tsv

###############################################################################
# Done
###############################################################################

echo "[INFO] Benchmark completed successfully."
