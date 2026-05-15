#!/usr/bin/env bash

###############################################################################
# Haplotype block inference using PLINK
#
# Usage:
#   bash 01.plink_hapblock.sh <chromosome> <variant_type> <group>
#
# Example:
#   bash 01.plink_hapblock.sh chr1 SNPs EAS
#
# Description:
#   This script identifies haplotype blocks from QC-passed PLINK datasets
#   using PLINK's --blocks function.
#
# Requirements:
#   - PLINK v1.9 or newer
#
# Input:
#   block/<chr>/*.qc_passed.*
#
# Output:
#   block/<group>/
#       ├── *.blocks
#       ├── *.blocks.det
#       └── log files
###############################################################################

set -euo pipefail

# ------------------------------- Arguments ---------------------------------- #

if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 <chromosome> <variant_type> <group>"
    exit 1
fi

CHR=$1
VAR_TYPE=$2
GROUP=$3

# ------------------------------ Configuration ------------------------------- #

INPUT_PREFIX="block/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.dip.${VAR_TYPE}.qc_passed"

OUTDIR="block/${GROUP}"

THREADS=8

mkdir -p "${OUTDIR}"

# ----------------------------- Input checking ------------------------------- #

for ext in bed bim fam; do
    [[ -f "${INPUT_PREFIX}.${ext}" ]] || {
        echo "[ERROR] Missing input file: ${INPUT_PREFIX}.${ext}"
        exit 1
    }
done

# --------------------------------- Logging ---------------------------------- #

echo "============================================================"
echo "Haplotype block inference"
echo "Chromosome : ${CHR}"
echo "Variant type : ${VAR_TYPE}"
echo "Group : ${GROUP}"
echo "Start time : $(date)"
echo "============================================================"

# -------------------------- Haplotype block calling ------------------------- #

plink \
    --bfile "${INPUT_PREFIX}" \
    --chr "${CHR}" \
    --blocks-min-maf 0.05 \
    --blocks-max-kb 10000 \
    --blocks 'no-pheno-req' \
    --threads "${THREADS}" \
    --out "${OUTDIR}/${CHR}.${GROUP}.ld_block"

# -------------------------------- Completion -------------------------------- #

echo
echo "============================================================"
echo "Haplotype block inference completed successfully"
echo "Output prefix : ${OUTDIR}/${CHR}.${GROUP}.ld_block"
echo "Finish time   : $(date)"
echo "============================================================"
