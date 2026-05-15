#!/usr/bin/env bash

# Annotate non-orphan SVs in difficult/easy regions

set -euo pipefail

VCF="CHM13-APGp1_MC.autosome.ex_tc.final.dip.SVs.vcf.gz"

GIAB="CHM13_alldifficultregions.bed.gz"
T2T="hs1.combined_mask.bed"
PANMASK="chm13.pm151a-v2.easy.bed"

echo -e "r2\ttotal\tgiab_diff\tgiab_easy\tt2tmask_diff\tpanmask_diff\tpanmask_easy"

for r2 in 0.2 0.4; do

    FILE="sv.non_orphan.${r2}.1000000k.tsv"

    TOTAL=$(wc -l < "${FILE}")

    awk '{print $3}' "${FILE}" > query.tmp

    GIAB_DIFF=$(bcftools view -H -T "${GIAB}" -i "ID=@query.tmp" "${VCF}" | wc -l)
    T2T_DIFF=$(bcftools view -H -T ^"${T2T}" -i "ID=@query.tmp" "${VCF}" | wc -l)
    PAN_DIFF=$(bcftools view -H -T ^"${PANMASK}" -i "ID=@query.tmp" "${VCF}" | wc -l)

    echo -e "${r2}\t${TOTAL}\t${GIAB_DIFF}\t$((TOTAL-GIAB_DIFF))\t${T2T_DIFF}\t${PAN_DIFF}\t$((TOTAL-PAN_DIFF))"

    rm -f query.tmp

done
