#!/usr/bin/env bash

# Count short variants tagged by SVs

set -euo pipefail

CHR=$1
LD=$2
WINDOW_KB=1000000

VCF="CHM13-APGp1.chr${CHR}.recalibrated.addtag.filter.ex_tc.vcf.gz"

awk -v chr="${CHR}" '$1==chr' "sv.${WINDOW_KB}k.${LD}.short.tagged.site.dir.txt" | \
while read -r chr num site; do

    FILE="ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.lite.${LD}.txt"

    awk -v ld="${LD}" '$4>ld' \
        "ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.lite.txt" \
        > "${FILE}"

    cut -f1-2 "${FILE}" > "${FILE}.chr_pos"

    COUNT=$(bcftools view -H \
        -R "${FILE}.chr_pos" \
        "${VCF}" | wc -l)

    echo -e "${chr}\t${num}\t${site}\t${COUNT}"

done > "ld/${CHR}.sv.${WINDOW_KB}k.ld.short.${LD}.stat.sr_tag.txt"
