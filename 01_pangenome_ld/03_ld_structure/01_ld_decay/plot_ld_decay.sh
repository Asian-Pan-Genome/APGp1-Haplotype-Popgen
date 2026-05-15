#!/usr/bin/env bash

###############################################################################
# Plot LD decay curves for multiple variant callsets using PopLDdecay
#
# Usage:
#   bash plot_ld_decay.sh <chromosome>
#
# Example:
#   bash plot_ld_decay.sh chr1
#
# Description:
#   This script generates multi-population LD decay plots from PopLDdecay
#   statistics files using Plot_MultiPop.pl.
#
# Requirements:
#   - PopLDdecay (>= v3.42)
#   - Perl
#   - R environment
#
# Expected directory structure:
#   ├── ngs/
#   ├── pangenome/
#   ├── phased_ngs/
#   └── plot_ld_decay.sh
#
# Input files:
#   ngs/CHM13-APGp1.<chr>...LDdecay.stat.gz
#   pangenome/CHM13-APGp1_MC.<chr>...LDdecay.stat.gz
#   phased_ngs/CHM13-APGp1.<chr>...LDdecay.stat.gz
#
# Output:
#   <chr>.pdf
#   <chr>.png
#   <chr>.list
###############################################################################

set -euo pipefail

# ------------------------------- Arguments ---------------------------------- #

if [[ $# -ne 1 ]]; then
    echo "Usage: bash $0 <chromosome>"
    exit 1
fi

CHR=$1

# ------------------------------ Configuration ------------------------------- #

NGS_DIR="ngs"
PANGENOME_DIR="pangenome"
PHASED_DIR="phased_ngs"

PLOT_SCRIPT="${HOME}/biotools/PopLDdecay-3.42/bin/Plot_MultiPop.pl"

# Activate R environment if needed
source "${HOME}/anaconda3/bin/activate" R4.3.1

# ------------------------------ Input files -------------------------------- #

STAT_NGS="${NGS_DIR}/CHM13-APGp1.${CHR}.recalibrated.addtag.filter.panmasked.biallelic.snp.100k.LDdecay.stat.gz"

STAT_PANGENOME="${PANGENOME_DIR}/CHM13-APGp1_MC.${CHR}.ex_tc.final.SNPs.biallelic.addtag.filter.dip.100k.LDdecay.stat.gz"

STAT_PHASED="${PHASED_DIR}/CHM13-APGp1.${CHR}.recalibrated.addtag.filter.panmasked.biallelic.snp.phased.100k.LDdecay.stat.gz"

# ----------------------------- File checking -------------------------------- #

for file in \
    "${STAT_NGS}" \
    "${STAT_PANGENOME}" \
    "${STAT_PHASED}"
do
    [[ -f "${file}" ]] || {
        echo "[ERROR] Missing input file: ${file}"
        exit 1
    }
done

[[ -f "${PLOT_SCRIPT}" ]] || {
    echo "[ERROR] Plot_MultiPop.pl not found: ${PLOT_SCRIPT}"
    exit 1
}

# ------------------------------ Prepare list -------------------------------- #

LIST_FILE="${CHR}.list"

cat > "${LIST_FILE}" <<EOF
${STAT_NGS}	NGS
${STAT_PANGENOME}	Pangenome
${STAT_PHASED}	Phased_NGS
EOF

# --------------------------------- Logging --------------------------------- #

echo "============================================================"
echo "LD decay plotting started"
echo "Chromosome : ${CHR}"
echo "List file  : ${LIST_FILE}"
echo "Start time : $(date)"
echo "============================================================"

# ------------------------------- Plotting ----------------------------------- #

perl "${PLOT_SCRIPT}" \
    -inList "${LIST_FILE}" \
    --out "${CHR}"

# Optional:
# Add '-measure both' if both r^2 and D' are needed

# -------------------------------- Completion -------------------------------- #

echo
echo "============================================================"
echo "LD decay plotting completed successfully"
echo "Output prefix : ${CHR}"
echo "Finish time   : $(date)"
echo "============================================================"
