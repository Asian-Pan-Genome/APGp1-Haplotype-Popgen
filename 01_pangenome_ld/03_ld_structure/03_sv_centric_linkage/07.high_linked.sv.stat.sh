#!/usr/bin/env bash

# Get top LD-linked short variants for each SV

set -euo pipefail

CHR=$1
WINDOW_KB=$2

OUT="ld/${CHR}.sv.${WINDOW_KB}k.ld.short.txt"

awk -v chr="${CHR}" '$1==chr' sv.site.dir.txt | \
while read -r chr num site; do

    awk 'NR==FNR{a[$1];next}!($6 in a)' sv.site.id.txt \
        "ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.txt" | \
        sort -k7,7gr | head -n 1

done > "${OUT}"
