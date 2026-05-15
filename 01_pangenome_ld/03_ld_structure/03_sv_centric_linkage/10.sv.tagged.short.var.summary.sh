#!/usr/bin/env bash

# Summarize GWAS/QTL annotations for SV-tagged variants

set -euo pipefail

CHR=$1
LD=$2
WINDOW_KB=1000000

GWAS="1.gwas/${LD}/sv.tagged.short.var.gwas.txt"
EQTL="2.qtl/${LD}/eQTL/sv.tagged.short.var.eQTL.txt"
SQTL="2.qtl/${LD}/sQTL/sv.tagged.short.var.sQTL.txt"

awk -v chr="${CHR}" '$1==chr' "sv.${WINDOW_KB}k.${LD}.short.tagged.site.dir.txt" | \
while read -r chr num site; do

    FILE="ld/${chr}/${num}/${chr}_${num}_${WINDOW_KB}k.ld.lite.${LD}.txt"

    GWAS_N=$(awk 'NR==FNR{a[$3];next}($1 in a)' "${FILE}" "${GWAS}" | wc -l)
    EQTL_N=$(awk 'NR==FNR{a[$3];next}($1 in a)' "${FILE}" "${EQTL}" | wc -l)
    SQTL_N=$(awk 'NR==FNR{a[$3];next}($1 in a)' "${FILE}" "${SQTL}" | wc -l)

    echo -e "${chr}\t${num}\t${site}\t${GWAS_N}\t${EQTL_N}\t${SQTL_N}"

done > "ld/${CHR}.sv.${WINDOW_KB}k.ld.short.${LD}.stat.gwas_qtl.txt"
