#!/usr/bin/env bash

# Intersect SV-tagged short variants with GWAS catalog

set -euo pipefail

LD=$1

VCF="sv.1000000k.${LD}.short.tagged.var.vcf.gz"
GWAS="chm13v2.0_GWASv1.0rsids_e100_r2022-03-08.vcf.gz"

OUT="1.gwas/${LD}"

mkdir -p "${OUT}"

bcftools isec \
    -c all \
    "${VCF}" \
    "${GWAS}" \
    -p "${OUT}" \
    -Ob \
    --threads 4 \
    -W
