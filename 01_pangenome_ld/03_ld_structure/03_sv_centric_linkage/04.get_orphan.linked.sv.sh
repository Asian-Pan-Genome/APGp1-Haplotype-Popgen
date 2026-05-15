#!/usr/bin/env bash

# Identify orphan SVs with exactly one high-LD SNP

set -euo pipefail

OUT="sv.orphan_high_ld_snp.tsv"

> "${OUT}"

while read -r chr num site; do

    FILE="ld/${chr}/${num}/${chr}_${num}_500k.ld"

    COUNT=$(sed '1d' "${FILE}" | awk '$3!=$6 && $NF>0.8' | wc -l)

    if [[ ${COUNT} -eq 1 ]]; then

        sed '1d' "${FILE}" | \
        awk '$3!=$6 && $NF>0.8' | \
        awk '{
            d=$5-$2;
            if(d<0)d=-d;
            print $0"\t"d+1
        }'

    fi

done < sv.site.dir.txt > "${OUT}"
