#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Extract eQTLs overlapping candidate structural variants
#
# This script:
#   1. Merges Truvari true-positive SVs and alternatively
#      represented SV candidates
#   2. Intersects candidate SVs with 1KCP eQTL summary statistics
#   3. Extracts unique overlapping SV identifiers
#   4. Retrieves the top-ranked eQTL record for each SV
#
# Requirements:
#   - awk
#   - sort
#   - zcat
#
# Usage:
#   bash 03.extract_overlapping_sv_eqtls.sh \
#       truvari_out_hg38 \
#       alternatively_represented.pair.sv.tsv \
#       sv.pos.sel.eqtl.summary.tsv
#
###############################################################################

TRUVARI_DIR=$1
ALT_REPRESENTED_SV=$2
EQTL_SUMMARY=$3

###############################################################################
# Step 1. Merge overlapping SV candidates
###############################################################################

echo "[INFO] Generating overlapping SV candidate list..."

cat \
    <(cut -f1-3 ${TRUVARI_DIR}/tp.sv.tsv) \
    <(
        awk 'BEGIN{OFS="\t"}{
            print $10, $7, $9
        }' ${ALT_REPRESENTED_SV} | \
        sed 's/\tchr/\t/g'
    ) | \
sort -k1,1V > all.overlapped.candidates.sv.tsv

###############################################################################
# Step 2. Intersect with chromosome-level 1KCP eQTL summaries
###############################################################################

echo "[INFO] Intersecting SVs with 1KCP eQTL summaries..."

for chr in chr{1..22}; do
    awk '
        NR==FNR{
            a[$1 FS $2 FS $3]
            next
        }
        ($1 FS $2 FS $3) in a
    ' all.overlapped.candidates.sv.tsv \
      <(zcat 1kcp.eqtl.${chr}.summary.gz)
done > all.overlapped.candidates.sv.pos.sel.eqtl.summary.tsv

###############################################################################
# Step 3. Extract unique overlapping SV identifiers
###############################################################################

echo "[INFO] Extracting unique SV identifiers..."

cut -f1 \
    all.overlapped.candidates.sv.pos.sel.eqtl.summary.tsv | \
sort -k1,1V -u \
> all.overlapped.candidates.sv.pos.sel.eqtl.record.txt

###############################################################################
# Step 4. Retrieve top-ranked eQTL record per SV
###############################################################################

echo "[INFO] Selecting top-ranked eQTL record for each SV..."

while read -r sv; do
    awk -v sv="${sv}" '
        $1 == sv
    ' ${EQTL_SUMMARY} | \
    sort -k14,14g | \
    head -n 1
done < all.overlapped.candidates.sv.pos.sel.eqtl.record.txt \
> all.overlapped.candidates.sv.pos.sel.eqtl.record.final.tsv

###############################################################################
# Done
###############################################################################

echo "[INFO] eQTL extraction completed successfully."
