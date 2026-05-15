#!/usr/bin/env bash

###############################################################################
# PopLDdecay wrapper for population-specific LD decay analysis
#
# Usage:
#   bash popLDdecay.subset.sh <vcf.gz> <output_dir> <subpop.txt> <max_dist_kb>
#
# Example:
#   bash popLDdecay.subset.sh \
#       input.vcf.gz \
#       results/ \
#       EUR.samples.txt \
#       500
#
# Requirements:
#   - PopLDdecay
#   - bgzip-compressed VCF (*.vcf.gz)
#
# Inputs:
#   1. VCF file (bgzip compressed)
#   2. Output directory
#   3. Sample list for subpopulation analysis
#   4. Maximum LD distance (kb)
#
# Output:
#   <output_dir>/<vcf_prefix>.<max_dist_kb>k.LDdecay.stat.gz
###############################################################################

set -euo pipefail

# ------------------------------- Arguments ---------------------------------- #

if [[ $# -ne 4 ]]; then
    echo "Usage: bash $0 <vcf.gz> <output_dir> <subpop.txt> <max_dist_kb>"
    exit 1
fi

VCF=$1
OUTDIR=$2
SUBPOP=$3
MAX_DIST=$4

# ----------------------------- Input checking ------------------------------- #

[[ -f "${VCF}" ]] || { echo "[ERROR] VCF file not found: ${VCF}"; exit 1; }
[[ -f "${SUBPOP}" ]] || { echo "[ERROR] Subpopulation file not found: ${SUBPOP}"; exit 1; }

mkdir -p "${OUTDIR}"

PREFIX=$(basename "${VCF}" .vcf.gz)
OUTFILE="${OUTDIR}/${PREFIX}.${MAX_DIST}k.LDdecay.stat.gz"

# --------------------------------- Logging --------------------------------- #

echo "============================================================"
echo "PopLDdecay analysis started"
echo "Start time : $(date)"
echo "VCF        : ${VCF}"
echo "SubPop     : ${SUBPOP}"
echo "MaxDist    : ${MAX_DIST} kb"
echo "Output     : ${OUTFILE}"
echo "============================================================"

# ------------------------------ Run analysis -------------------------------- #

PopLDdecay \
    -InVCF "${VCF}" \
    -SubPop "${SUBPOP}" \
    -MaxDist "${MAX_DIST}" \
    -MAF 0.01 \
    -OutType 2 \
    -OutStat "${OUTFILE}"

# -------------------------------- Completion -------------------------------- #

echo
echo "============================================================"
echo "PopLDdecay analysis completed successfully"
echo "Finish time: $(date)"
echo "============================================================"
