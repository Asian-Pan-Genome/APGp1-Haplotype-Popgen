#!/usr/bin/env bash

# Calculate pairwise LD around SV sites using PLINK

set -euo pipefail

CHR=$1
WINDOW_KB=$2

VCF="./plink/CHM13-APGp1_MC.autosome.ex_tc.final.dip.SNPs_INDELs_SVs.qc_passed"

echo "[INFO] $(date)"
echo "[INFO] Calculating LD for ${CHR}"

awk -v chr="${CHR}" '$1==chr' sv.site.dir.txt | \
while read -r chr num site; do

    OUTDIR="ld/${chr}/${num}"
    PREFIX="${OUTDIR}/${chr}_${num}_${WINDOW_KB}k"

    mkdir -p "${OUTDIR}"

    plink \
        --bfile "${VCF}" \
        --r2 \
        --ld-snps "${site}" \
        --ld-window-kb "${WINDOW_KB}" \
        --ld-window 99999 \
        --ld-window-r2 0 \
        --threads 16 \
        --out "${PREFIX}"

done

echo "[INFO] Finished: $(date)"
