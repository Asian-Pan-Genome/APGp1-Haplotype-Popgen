#!/usr/bin/env bash

# Extract lightweight LD tables

set -euo pipefail

CHR=$1
WINDOW_KB=$2

echo "[INFO] $(date)"

awk -v chr="${CHR}" '$1==chr' sv.site.dir.txt | \
while read -r chr num site; do

    INPUT="ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.txt"
    OUTPUT="ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.lite.txt"

    awk '{print $4,$5,$6,$7,$8,$9}' OFS="\t" "${INPUT}" > "${OUTPUT}"

done

echo "[INFO] Finished"
