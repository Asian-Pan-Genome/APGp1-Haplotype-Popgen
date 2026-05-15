#!/usr/bin/env bash

###############################################################################
# Plot multi-population LD decay curves using PopLDdecay
#
# Usage:
#   bash plot_ld_decay_multi_pop.sh
#
# Description:
#   This script generates LD decay comparison plots across multiple
#   populations and datasets (assembly-based vs short-read phased data)
#   using PopLDdecay's Plot_MultiPop.pl utility.
#
# Requirements:
#   - PopLDdecay >= v3.42
#   - Perl
#   - Conda environment with R installed
#
# Input:
#   A tab-delimited list file containing:
#       <LDdecay.stat.gz> <label>
#
# Example:
#   example.stat.gz    AFR
#   another.stat.gz    EUR
#
# Output:
#   assem_full.*
#       ├── assem_full.pdf
#       ├── assem_full.png
#       └── intermediate plotting files
###############################################################################

set -euo pipefail

# ------------------------------ Configuration ------------------------------- #

LIST_FILE="plot.file.full.list"
OUTPUT_PREFIX="assem_full"

PLOT_SCRIPT="${HOME}/biotools/PopLDdecay-3.42/bin/Plot_MultiPop.pl"

# Activate R environment
source "${HOME}/anaconda3/bin/activate" R4.3.1

# ----------------------------- Dependency check ----------------------------- #

[[ -f "${LIST_FILE}" ]] || {
    echo "[ERROR] Input list file not found: ${LIST_FILE}"
    exit 1
}

[[ -f "${PLOT_SCRIPT}" ]] || {
    echo "[ERROR] Plot_MultiPop.pl not found: ${PLOT_SCRIPT}"
    exit 1
}

# ------------------------------ Input check --------------------------------- #

echo "Checking LD decay statistics files..."

while IFS=$'\t' read -r stat_file label; do
    [[ -f "${stat_file}" ]] || {
        echo "[ERROR] Missing file for ${label}: ${stat_file}"
        exit 1
    }
done < "${LIST_FILE}"

echo "All input files detected."

# --------------------------------- Logging ---------------------------------- #

echo "============================================================"
echo "Multi-population LD decay plotting"
echo "Input list  : ${LIST_FILE}"
echo "Output      : ${OUTPUT_PREFIX}"
echo "Start time  : $(date)"
echo "============================================================"

# -------------------------------- Plotting ---------------------------------- #

perl "${PLOT_SCRIPT}" \
    -inList "${LIST_FILE}" \
    --out "${OUTPUT_PREFIX}"

# Optional:
# Add '-measure both' to plot both r^2 and D'

# ------------------------------- Completion --------------------------------- #

echo
echo "============================================================"
echo "LD decay plotting completed successfully"
echo "Output prefix : ${OUTPUT_PREFIX}"
echo "Finish time   : $(date)"
echo "============================================================"
