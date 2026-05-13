#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Run population branch statistic (PBS) selection scans
#
# This script:
#   1. Loads polarized biallelic variant datasets
#   2. Defines focal, sister, and outgroup populations
#   3. Runs PBS scans using pbscan
#
# Requirements:
#   - pbscan
#
# Usage:
#   bash 01.run_pbs_selection_scan_site.sh \
#       variant_type \
#       polarization_label \
#       output_prefix \
#       focal_population \
#       sister_population \
#       outgroup_population \
#       maf_threshold
#
# Example:
#   bash 01.run_pbs_selection_scan_site.sh \
#       SVs \
#       polarized.masked \
#       CHA \
#       CHA \
#       non_EAS_AFR \
#       AFR \
#       0.05
#
###############################################################################

VAR_TYPE=$1
POLARIZATION=$2
PREFIX=$3
FOCAL_POP=$4
SISTER_POP=$5
OUTGROUP_POP=$6
MAF=$7

###############################################################################
# Input configuration
###############################################################################

INPUT_VCF="vcfwave_merge/4.polarize/CHM13-APGp1-HPRCp1-HGSVCp3_MC.ex_tc.final.${VAR_TYPE}.biallelic.addtag.filter.unrelated.${POLARIZATION}.vcf"

FOCAL_SAMPLES=./${FOCAL_POP}.txt
SISTER_SAMPLES=./${SISTER_POP}.txt
OUTGROUP_SAMPLES=./${OUTGROUP_POP}.txt

OUTDIR=./${VAR_TYPE}/${FOCAL_POP}/${MAF}

MC_REPLICATES=10000

###############################################################################
# Start
###############################################################################

date

echo "[INFO] Variant type: ${VAR_TYPE}"
echo "[INFO] Polarization: ${POLARIZATION}"
echo "[INFO] Focal population: ${FOCAL_POP}"
echo "[INFO] Sister population: ${SISTER_POP}"
echo "[INFO] Outgroup population: ${OUTGROUP_POP}"
echo "[INFO] MAF threshold: ${MAF}"

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Run PBS scan
###############################################################################

echo "[INFO] Running PBS selection scan..."

pbscan \
    -vcf ${INPUT_VCF} \
    -maf ${MAF} \
    -pop1 ${FOCAL_SAMPLES} \
    -pop2 ${SISTER_SAMPLES} \
    -pop3 ${OUTGROUP_SAMPLES} \
    -out ${OUTDIR}/${PREFIX} \
    -mc ${MC_REPLICATES}

###############################################################################
# Done
###############################################################################

echo "[INFO] PBS scan completed successfully."

date
