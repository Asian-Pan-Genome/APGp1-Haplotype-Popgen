#!/usr/bin/env bash

# Intersect SV-tagged short variants with QTL datasets

set -euo pipefail

LD=$1
CHR=$2
QTL=$3

VCF="sv.1000000k.${LD}.short.tagged.var.vcf.gz"
QTL_VCF="${CHR}.${QTL}.final.vcf.gz"

OUT="2.qtl/${LD}/${QTL}/${CHR}"

mkdir -p "${OUT}"

bcftools isec \
    -c none \
    "${VCF}" \
    "${QTL_VCF}" \
    -p "${OUT}" \
    -Ob \
    -n=2 \
    --threads 2 \
    -W
