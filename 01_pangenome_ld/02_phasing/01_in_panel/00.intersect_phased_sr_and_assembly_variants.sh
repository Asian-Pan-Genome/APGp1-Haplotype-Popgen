#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Intersect phased and assembly-based variant callsets
#
# This script:
#   1. Compares phased SR-based variant calls against assembly-derived variants
#   2. Identifies shared variants using bcftools isec
#   3. Outputs intersected callsets in BCF format
#
# Requirements:
#   - bcftools
#
# Usage:
#   bash 00.intersect_phased_sr_and_assembly_variants.sh \
#       chr1
#
###############################################################################

CHR=$1

###############################################################################
# Input configuration
###############################################################################

PHASED_VCF="SR/02.chm13v2/ex_tc/${CHR}/CHM13-APGp1.${CHR}.recalibrated.addtag.filter.ex_tc.phased.bcf"

ASSEMBLY_VCF="Assembly/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.SNPs_INDELs_SVs.filtered.norm.rename.bcf"

OUTDIR=./intersected_site_vcf/${CHR}

THREADS=8

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Intersect phased and assembly-based variants
###############################################################################

echo "[INFO] Running bcftools isec..."

bcftools isec \
    ${PHASED_VCF} \
    ${ASSEMBLY_VCF} \
    -p ${OUTDIR}/isec \
    -c none \
    -n =2 \
    -Ob \
    --threads ${THREADS} \
    -W=tbi

###############################################################################
# Done
###############################################################################

echo "[INFO] Variant intersection completed successfully."

date
