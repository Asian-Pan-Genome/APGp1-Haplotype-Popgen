#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run window-based PBS selection scans
#
# This script:
#   1. Loads chromosome-level polarized biallelic variants
#   2. Defines focal, sister, and outgroup populations
#   3. Performs sliding-window PBS scans using pbscan
#
# Requirements:
#   - pbscan
#
# Usage:
#   bash 02.run_pbs_selection_scan_window.sh \
#       chr1 \
#       variant_type \
#       polarization_label \
#       output_prefix \
#       focal_population \
#       sister_population \
#       outgroup_population
#
# Example:
#   bash 02.run_pbs_selection_scan_window.sh \
#       chr1 \
#       SNPs \
#       polarized.masked \
#       CHA \
#       CHA \
#       non_EAS_AFR \
#       AFR
#
###############################################################################

CHR=$1
VAR_TYPE=$2
POLARIZATION=$3
PREFIX=$4
FOCAL_POP=$5
SISTER_POP=$6
OUTGROUP_POP=$7

###############################################################################
# Input configuration
###############################################################################

INPUT_VCF="vcfwave_merge/4.polarize/${CHR}/CHM13-APGp1-HPRCp1-HGSVCp3_MC.${CHR}.ex_tc.final.${VAR_TYPE}.biallelic.addtag.filter.unrelated.${POLARIZATION}.vcf"

FOCAL_SAMPLES=./${FOCAL_POP}.txt
SISTER_SAMPLES=./${SISTER_POP}.txt
OUTGROUP_SAMPLES=./${OUTGROUP_POP}.txt

OUTDIR=./${VAR_TYPE}/${FOCAL_POP}/${CHR}

###############################################################################
# PBS scan parameters
###############################################################################

MAF=0.05
WINDOW_SIZE=100
STEP_SIZE=50
MIN_VARIANTS=5

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Chromosome: ${CHR}"
echo "[INFO] Variant type: ${VAR_TYPE}"
echo "[INFO] Polarization: ${POLARIZATION}"
echo "[INFO] Focal population: ${FOCAL_POP}"
echo "[INFO] Sister population: ${SISTER_POP}"
echo "[INFO] Outgroup population: ${OUTGROUP_POP}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Run sliding-window PBS scan
###############################################################################

echo "[INFO] Running windowed PBS scan..."

pbscan \
    -vcf ${INPUT_VCF} \
    -maf ${MAF} \
    -pop1 ${FOCAL_SAMPLES} \
    -pop2 ${SISTER_SAMPLES} \
    -pop3 ${OUTGROUP_SAMPLES} \
    -out ${OUTDIR}/${PREFIX} \
    -win ${WINDOW_SIZE} \
    -step ${STEP_SIZE} \
    -min ${MIN_VARIANTS}

###############################################################################
# Done
###############################################################################

echo "[INFO] Windowed PBS scan completed successfully."

date
