#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Identify alternatively represented structural variants
#
# This script:
#   1. Extracts false-negative (FN) and false-positive (FP) SVs
#      from Truvari benchmark results
#   2. Converts SV coordinates to BED format
#   3. Identifies nearby discordant SV pairs
#   4. Detects alternatively represented SVs based on
#      breakpoint proximity and SV length similarity
#
# Requirements:
#   - bcftools
#   - bedtools
#   - awk
#
# Usage:
#   bash 02.get_alternatively_represented_svs.sh \
#       truvari_out_hg38
#
###############################################################################

TRUVARI_DIR=$1

###############################################################################
# Step 1. Extract FN and FP SV information
###############################################################################

echo "[INFO] Extracting false-negative SVs..."

bcftools query \
    -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' \
    ${TRUVARI_DIR}/fn.vcf.gz > fn.txt

echo "[INFO] Extracting false-positive SVs..."

bcftools query \
    -f '%CHROM\t%POS\t%ID\t%INFO/EndDistance\n' \
    ${TRUVARI_DIR}/fp.vcf.gz > fp.txt

###############################################################################
# Step 2. Convert to BED format
###############################################################################

echo "[INFO] Converting SVs to BED format..."

awk 'BEGIN{OFS="\t"}{
    print $1, $2-1, $2, $3, $4, $5
}' fn.txt > fn.bed

awk 'BEGIN{OFS="\t"}{
    print $1, $2-1, $2, $3, $4
}' fp.txt > fp.bed

###############################################################################
# Step 3. Identify nearby discordant SVs
###############################################################################

echo "[INFO] Searching for nearby discordant SVs..."

bedtools window \
    -w 100 \
    -a fn.bed \
    -b fp.bed > nearby.discordant.sv.txt

###############################################################################
# Step 4. Detect alternatively represented SV pairs
###############################################################################

echo "[INFO] Identifying alternatively represented SVs..."

awk 'BEGIN{OFS="\t"}{
    a = ($6 < 0 ? -$6 : $6)
    b = ($11 < 0 ? -$11 : $11)
    diff = (a > b ? a-b : b-a)

    print $0, diff
}' nearby.discordant.sv.txt | \
awk '$NF < 100' > alternatively_represented.pair.sv.tsv

###############################################################################
# Done
###############################################################################

echo "[INFO] Analysis completed successfully."
