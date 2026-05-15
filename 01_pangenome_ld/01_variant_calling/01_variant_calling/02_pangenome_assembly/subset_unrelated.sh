#!/usr/bin/env bash

###############################################################################
# Subset unrelated samples from filtered SV callset
#
# Description:
#   This script extracts unrelated individuals from the final filtered
#   biallelic SV callset and recalculates population genetic summary tags.
#
#   Two subset modes are supported:
#     1. unrelated       -> 258 unrelated samples
#     2. unrelated_CHA   -> 99 unrelated CHA samples
#
# Usage:
#   bash subset_unrelated.sh <var_type> <chr> <group>
#
# Arguments:
#   var_type   Variant type (e.g. ins, del, snv)
#   chr        Chromosome name (e.g. chr1)
#   group      Subset group:
#                - unrelated
#                - unrelated_CHA
#
# Example:
#   bash subset_unrelated.sh del chr1 unrelated
#   bash subset_unrelated.sh ins chr22 unrelated_CHA
#
# Requirements:
#   - bcftools >= 1.16
#
# Input:
#   *.biallelic.addtag.filter.vcf.gz
#
# Output:
#   *.unrelated.vcf.gz
#   *.unrelated.CHA.vcf.gz
#
###############################################################################

set -euo pipefail

############################
# Arguments
############################

if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 <var_type> <chr> <group>"
    echo
    echo "Groups:"
    echo "  unrelated"
    echo "  unrelated_CHA"
    exit 1
fi

VAR_TYPE="$1"
CHR="$2"
GROUP="$3"

############################
# Configuration
############################

THREADS=16
OUT_DIR="./${CHR}"

INPUT_VCF="${OUT_DIR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${VAR_TYPE}.biallelic.addtag.filter.vcf.gz"

############################
# Select sample subset
############################

case "${GROUP}" in
    unrelated)
        SAMPLE_LIST="./258.unrelated.sample.txt"
        SUFFIX="unrelated"
        ;;

    unrelated_CHA)
        SAMPLE_LIST="./99.unrelated.CHA.sample.txt"
        SUFFIX="unrelated.CHA"
        ;;

    *)
        echo "Error: invalid group '${GROUP}'"
        echo "Valid groups: unrelated | unrelated_CHA"
        exit 1
        ;;
esac

OUTPUT_VCF="${OUT_DIR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${VAR_TYPE}.biallelic.addtag.filter.${SUFFIX}.vcf.gz"

############################
# Start
############################

echo "[$(date)] Starting unrelated sample subsetting"
echo "Chromosome : ${CHR}"
echo "Variant type : ${VAR_TYPE}"
echo "Group : ${GROUP}"
echo

###############################################################################
# Subset samples and recalculate INFO tags
###############################################################################

bcftools view \
    -S "${SAMPLE_LIST}" \
    -c 1 \
    "${INPUT_VCF}" \
    --threads "${THREADS}" \
    | bcftools +fill-tags \
        --threads "${THREADS}" \
        -- \
        -t AN,AC,AF,MAF,F_MISSING,MAC:1=MAC,HWE \
    | bcftools view \
        -Oz \
        -o "${OUTPUT_VCF}" \
        --threads "${THREADS}" \
        -W=tbi

###############################################################################
# Done
###############################################################################

echo
echo "[$(date)] Subsetting completed successfully."
echo
echo "Output:"
echo "  ${OUTPUT_VCF}"
echo
echo "Finished."
echo
