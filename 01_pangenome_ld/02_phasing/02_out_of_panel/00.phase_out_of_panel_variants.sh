#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Phase out-of-panel variants using in-panel scaffolds
#
# This script:
#   1. Extracts overlapping in-panel scaffold variants
#   2. Removes missing sites from scaffold haplotypes
#   3. Performs statistical phasing using phase_common_static in SHAPEIT5
#   4. Extracts out-of-panel phased variants
#
# Requirements:
#   - bcftools
#   - phase_common_static (SHAPEIT5)
#
# Usage:
#   bash phase_out_of_panel_variants.sh \
#       chr1
#
###############################################################################

CHR=$1

###############################################################################
# Input configuration
###############################################################################

INPUT_VCF=Assembly/${CHR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.SNPs_INDELs_SVs.filtered.norm.rename.bcf

OVERLAPPED_SITES=../01_in_panel/intersected_site_vcf/${CHR}/isec/overlapped.sites.txt

REFERENCE_PANEL=1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz

GENETIC_MAP=1KGP_chm13v2/t2t_native_scaled_maps/${CHR}.t2t.scaled.gmap.gz

###############################################################################
# Output configuration
###############################################################################

PREFIX=$(basename ${INPUT_VCF} .bcf)

SCAFFOLD_VCF=./${CHR}/${PREFIX}.1kg.in_panel.no_missing.bcf

PHASED_VCF=./${CHR}/${PREFIX}.1kg.in_panel_scaffolding.phased.full.bcf

OUT_OF_PANEL_VCF=./${CHR}/${PREFIX}.1kg.out_of_panel.phased.full.bcf

THREADS=24

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"

###############################################################################
# Step 1. Extract in-panel scaffold variants
###############################################################################

echo "[INFO] Extracting scaffold variants..."

bcftools view \
    -R ${OVERLAPPED_SITES} \
    ${INPUT_VCF} \
    --threads ${THREADS} | \
bcftools +fill-tags \
    --threads ${THREADS} \
    - \
    -- -t AN,AC,AF,MAF,F_MISSING,MAC:1=MAC,HWE | \
bcftools view \
    -i 'F_MISSING==0' \
    -Ob \
    -o ${SCAFFOLD_VCF} \
    --threads ${THREADS} \
    -W=tbi

###############################################################################
# Step 2. Perform statistical phasing
###############################################################################

echo "[INFO] Running statistical phasing..."

phase_common_static \
    --input ${INPUT_VCF} \
    --scaffold ${SCAFFOLD_VCF} \
    --output ${PHASED_VCF} \
    --thread ${THREADS} \
    --map ${GENETIC_MAP} \
    --region ${CHR}

###############################################################################
# Step 3. Extract out-of-panel variants
###############################################################################

echo "[INFO] Extracting out-of-panel phased variants..."

bcftools view \
    -T ^${OVERLAPPED_SITES} \
    ${PHASED_VCF} \
    -Ob \
    -o ${OUT_OF_PANEL_VCF} \
    --threads ${THREADS} \
    -W=tbi

###############################################################################
# Done
###############################################################################

echo "[INFO] Out-of-panel phasing completed successfully."

date
