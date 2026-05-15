#!/usr/bin/env bash

###############################################################################
# Variant statistics summary pipeline
#
# Description:
#   This script generates summary statistics for filtered biallelic variants,
#   including:
#
#     1. Total variant counts
#     2. Rare variant counts (AF/MAF thresholds)
#     3. Counts within easy genomic regions:
#          - Panmask easy regions
#          - GIAB easy regions
#
#   Supported variant classes:
#       - SNPs
#       - INDELs
#       - MNPs
#       - SVs
#
# Usage:
#
#   1. Generate primary site list
#      bash variant_stats.sh primary <var_type>
#
#   2. Generate MAF-based statistics
#      bash variant_stats.sh maf <chr>
#
#   3. Generate AF-based statistics
#      bash variant_stats.sh af <chr>
#
# Examples:
#
#   bash variant_stats.sh primary SVs
#   bash variant_stats.sh maf chr1
#   bash variant_stats.sh af chr22
#
# Requirements:
#   - bcftools >= 1.16
#
###############################################################################

set -euo pipefail

############################
# Arguments
############################

if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "  bash $0 primary <var_type>"
    echo "  bash $0 maf <chr>"
    echo "  bash $0 af <chr>"
    exit 1
fi

MODE="$1"
ARG="$2"

############################
# Configuration
############################

THREADS=16

PANMASK_BED=~/chm13v2_data/easy_region/panmask/chm13.pm151a-v2.easy.bed.gz
GIAB_BED=~/chm13v2_data/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz

VARIANT_TYPES=("SNPs" "INDELs" "MNPs" "SVs")

###############################################################################
# Function: count_variants
###############################################################################

count_variants () {

    local VCF="$1"
    local LABEL="$2"
    local FILTER_EXPR="$3"
    local REGION="$4"

    local CMD="bcftools view -H"

    if [[ -n "${REGION}" ]]; then
        CMD+=" -T ${REGION}"
    fi

    if [[ -n "${FILTER_EXPR}" ]]; then
        CMD+=" -i '${FILTER_EXPR}'"
    fi

    CMD+=" ${VCF}"

    eval "${CMD}" | wc -l | awk -v label="${LABEL}" '{print label"\t"$1}'
}

###############################################################################
# Mode 1. Generate primary site list
###############################################################################

if [[ "${MODE}" == "primary" ]]; then

    VAR_TYPE="${ARG}"

    echo "[$(date)] Generating primary site list"
    echo "Variant type : ${VAR_TYPE}"
    echo

    OUTPUT="${VAR_TYPE}.primary.site.txt"

    > "${OUTPUT}"

    for CHR in chr{1..22}; do

        INPUT_VCF="../${CHR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${VAR_TYPE}.biallelic.addtag.filter.vcf.gz"

        echo "Processing ${CHR}..."

        bcftools view -H "${INPUT_VCF}" >> "${OUTPUT}"

    done

    echo
    echo "[$(date)] Completed"
    echo "Output:"
    echo "  ${OUTPUT}"
    echo

    exit 0
fi

###############################################################################
# Mode 2 & 3. AF / MAF statistics
###############################################################################

CHR="${ARG}"

if [[ "${MODE}" == "maf" ]]; then
    TAG="MAF"
    OUTPUT="${CHR}.biallelic.stat.MAF.txt"

elif [[ "${MODE}" == "af" ]]; then
    TAG="AF"
    OUTPUT="${CHR}.biallelic.stat.AF.txt"

else
    echo "Error: invalid mode '${MODE}'"
    exit 1
fi

echo "[$(date)] Generating ${TAG}-based statistics"
echo "Chromosome : ${CHR}"
echo

###############################################################################
# Header
###############################################################################

echo -e "VariantType\tCategory\tCount" > "${OUTPUT}"

###############################################################################
# Statistics
###############################################################################

for TYPE in "${VARIANT_TYPES[@]}"; do

    VCF="../${CHR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${TYPE}.biallelic.addtag.filter.vcf.gz"

    echo "Processing ${TYPE}..."

    ########################################
    # Whole genome
    ########################################

    count_variants "${VCF}" "${TYPE}\t0.05"     "${TAG}<=0.05" "" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\t0.01"     "${TAG}<=0.01" "" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\tall"      ""              "" >> "${OUTPUT}"

    ########################################
    # Panmask easy regions
    ########################################

    count_variants "${VCF}" "${TYPE}\t0.05_panmask" "${TAG}<=0.05" "${PANMASK_BED}" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\t0.01_panmask" "${TAG}<=0.01" "${PANMASK_BED}" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\tall_panmask"  ""              "${PANMASK_BED}" >> "${OUTPUT}"

    ########################################
    # GIAB easy regions
    ########################################

    count_variants "${VCF}" "${TYPE}\t0.05_giab_easy" "${TAG}<=0.05" "${GIAB_BED}" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\t0.01_giab_easy" "${TAG}<=0.01" "${GIAB_BED}" >> "${OUTPUT}"
    count_variants "${VCF}" "${TYPE}\tall_giab_easy"  ""              "${GIAB_BED}" >> "${OUTPUT}"

done

###############################################################################
# Done
###############################################################################

echo
echo "[$(date)] Statistics generation completed."
echo
echo "Output:"
echo "  ${OUTPUT}"
echo
echo "Finished."
echo
