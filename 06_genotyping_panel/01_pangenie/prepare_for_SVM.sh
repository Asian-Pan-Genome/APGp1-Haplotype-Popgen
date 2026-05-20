#!/usr/bin/env bash
# ==============================================================================
# PanGenie merged-callset QC and concordance pipeline
# ==============================================================================
#
# Purpose:
#   Run downstream QC/concordance analyses on a complete merged PanGenie biallelic
#   VCF, without chromosome-level splitting and without Slurm submission.
#
# Default input:
#   pangenie_merged_bi_all.vcf.gz
#
# Expected helper scripts:
#   id_to_bubble.py
#   mendelian-consistency.py
#   collect-vcf-stats.py
#   genotype-concordance-variant.py
#   merge-tables.py
#
# Usage:
#   bash run_pangenie_merged_qc_clean.sh
#
# Optional:
#   bash run_pangenie_merged_qc_clean.sh /path/to/workdir
#
# Notes:
#   - Edit the variables in the "User-configurable variables" block as needed.
#   - The script runs directly as a bash workflow. It does not generate per-chromosome
#     work scripts and does not submit jobs to Slurm.
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. User-configurable variables
# ------------------------------------------------------------------------------

# Working directory. Can also be overridden by the first positional argument.
WORKDIR="${1:-/share/home/project/zhanglab/QYTMP/project/pangenome/CHM13-APGp1-HPRCp1-HGSVCp3/pangenie/hgsvc_hprc}"

# Input VCFs.
# MULTI_VCF is used to derive bubble-level statistics.
# INPUT_BI_VCF is the original/reference biallelic VCF used by downstream QC scripts.
# MERGED_BI_VCF is the complete merged PanGenie genotyping result to evaluate.
MULTI_VCF="${MULTI_VCF:-${WORKDIR}/input.multi.vcf}"
INPUT_BI_VCF="${INPUT_BI_VCF:-${WORKDIR}/input.bi.vcf.gz}"
MERGED_BI_VCF="${MERGED_BI_VCF:-${WORKDIR}/pangenie_merged_bi_all.vcf.gz}"

# Pedigree and extra sample list.
# PED_FILE should follow PLINK-style PED columns at minimum:
#   family_id sample_id paternal_id maternal_id ...
# ID_LIST is optional. If present, its second column is used when available;
# otherwise its first column is used.
PED_FILE="${PED_FILE:-${WORKDIR}/1KG/ped.list}"
ID_LIST="${ID_LIST:-${WORKDIR}/id.list}"

# Directory containing helper Python scripts.
# Keep helper script names clean; place all scripts in this directory or export SCRIPT_DIR.
SCRIPT_DIR="${SCRIPT_DIR:-$(pwd)}"

# Output prefix.
PREFIX="${PREFIX:-bi_all}"

# Threads for bcftools/tabix where applicable.
THREADS="${THREADS:-16}"

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

log() {
    echo "[$(date '+%F %T')] $*" >&2
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local f="$1"
    [[ -s "$f" ]] || die "Required file not found or empty: $f"
}

require_script() {
    local s="$1"
    [[ -s "${SCRIPT_DIR}/${s}" ]] || die "Required helper script not found: ${SCRIPT_DIR}/${s}"
}

# ------------------------------------------------------------------------------
# 3. Pre-flight checks
# ------------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

require_file "${MULTI_VCF}"
require_file "${INPUT_BI_VCF}"
require_file "${MERGED_BI_VCF}"
require_file "${PED_FILE}"

require_script "id_to_bubble.py"
require_script "mendelian-consistency.py"
require_script "collect-vcf-stats.py"
require_script "genotype-concordance-variant.py"
require_script "merge-tables.py"

command -v python3 >/dev/null 2>&1 || die "python3 not found in PATH"
command -v bcftools >/dev/null 2>&1 || die "bcftools not found in PATH"
command -v awk >/dev/null 2>&1 || die "awk not found in PATH"

log "Working directory: ${WORKDIR}"
log "MULTI_VCF:       ${MULTI_VCF}"
log "INPUT_BI_VCF:    ${INPUT_BI_VCF}"
log "MERGED_BI_VCF:   ${MERGED_BI_VCF}"
log "PED_FILE:        ${PED_FILE}"
log "ID_LIST:         ${ID_LIST}"
log "SCRIPT_DIR:      ${SCRIPT_DIR}"
log "PREFIX:          ${PREFIX}"

# ------------------------------------------------------------------------------
# 4. Bubble statistics from the multi-allelic input VCF
# ------------------------------------------------------------------------------

BUBBLE_STATS="bubble-statistics_${PREFIX}.tsv"

if [[ ! -s "${BUBBLE_STATS}" ]]; then
    log "Generating ${BUBBLE_STATS}"
    cat "${MULTI_VCF}" \
        | python3 "${SCRIPT_DIR}/id_to_bubble.py" \
        > "${BUBBLE_STATS}"
else
    log "Found existing ${BUBBLE_STATS}; skipping."
fi

# ------------------------------------------------------------------------------
# 5. Mendelian consistency
# ------------------------------------------------------------------------------

PED_SAMPLES="ped.samples.txt"
awk 'NR > 1 && NF >= 2 {print $2}' "${PED_FILE}" > "${PED_SAMPLES}"

log "Running Mendelian consistency"
python3 "${SCRIPT_DIR}/mendelian-consistency.py" statistics \
    -vcf "${MERGED_BI_VCF}" \
    -ped "${PED_FILE}" \
    -samples "${PED_SAMPLES}" \
    -table "mendelian-statistics_${PREFIX}.tsv" \
    -column-prefix pangenie \
    > "trio-statistics_${PREFIX}.tsv"

# ------------------------------------------------------------------------------
# 6. Create unrelated-sample VCF for genotyping statistics
# ------------------------------------------------------------------------------

UNRELATED_SAMPLES="unrelated.samples.txt"

awk 'NR > 1 && NF >= 4 && $3 == 0 && $4 == 0 {print $2}' "${PED_FILE}" > "${UNRELATED_SAMPLES}"

if [[ -s "${ID_LIST}" ]]; then
    awk 'NF >= 2 {print $2; next} NF == 1 {print $1}' "${ID_LIST}" >> "${UNRELATED_SAMPLES}"
else
    log "ID_LIST not found or empty; unrelated-sample list will be based only on PED founders."
fi

sort -u "${UNRELATED_SAMPLES}" -o "${UNRELATED_SAMPLES}"

log "Extracting unrelated samples from ${MERGED_BI_VCF}"
bcftools view \
    --samples-file "${UNRELATED_SAMPLES}" \
    --force-samples \
    "${MERGED_BI_VCF}" \
    --threads "${THREADS}" \
    -Oz \
    -o "unrelated-samples_${PREFIX}.vcf.gz"

tabix -f -p vcf "unrelated-samples_${PREFIX}.vcf.gz"

# ------------------------------------------------------------------------------
# 7. Genotyping statistics
# ------------------------------------------------------------------------------

log "Collecting VCF-level genotyping statistics"
python3 "${SCRIPT_DIR}/collect-vcf-stats.py" \
    "${INPUT_BI_VCF}" \
    "unrelated-samples_${PREFIX}.vcf.gz" \
    "${MERGED_BI_VCF}" \
    > "genotyping-statistics_${PREFIX}.tsv"

# ------------------------------------------------------------------------------
# 8. Self-genotyping concordance
# ------------------------------------------------------------------------------

# Samples are the intersection of PED samples and samples present in MULTI_VCF,
# preserving the order in MULTI_VCF.
VCF_SAMPLES="multi_vcf.samples.txt"
SELF_SAMPLES_FILE="self_concordance.samples.txt"

bcftools query -l "${MULTI_VCF}" > "${VCF_SAMPLES}"

awk 'NR == FNR {keep[$1] = 1; next} ($1 in keep) {print $1}' \
    "${PED_SAMPLES}" \
    "${VCF_SAMPLES}" \
    > "${SELF_SAMPLES_FILE}"

SELF_SAMPLES="$(paste -sd, "${SELF_SAMPLES_FILE}")"

if [[ -z "${SELF_SAMPLES}" ]]; then
    die "No overlapping samples found between ${PED_FILE} and ${MULTI_VCF}; cannot run genotype-concordance-variant.py."
fi

log "Running self-genotyping concordance for $(wc -l < "${SELF_SAMPLES_FILE}") samples"
python3 "${SCRIPT_DIR}/genotype-concordance-variant.py" \
    "${INPUT_BI_VCF}" \
    "${MERGED_BI_VCF}" \
    "self_${PREFIX}" \
    "${SELF_SAMPLES}" \
    pangenie_self-genotyping

# ------------------------------------------------------------------------------
# 9. Variant type annotation
# ------------------------------------------------------------------------------

log "Generating variant-type table"
bcftools query -f '%INFO/ID\t%REF\t%ALT\n' "${MERGED_BI_VCF}" \
    | awk -F '\t' '
        BEGIN {
            OFS = "\t";
            print "variant_id", "variant_type";
        }
        {
            id = $1;
            ref = $2;
            alt = $3;

            if (length(ref) < 50 && length(alt) < 50) {
                if (length(ref) == 1 && length(alt) == 1) {
                    type = "SNP";
                } else if (length(ref) == length(alt)) {
                    type = "MNP";
                } else if (length(ref) > length(alt) && length(alt) == 1) {
                    type = "DEL";
                } else if (length(alt) > 1 && length(ref) == 1) {
                    type = "INS";
                } else {
                    type = "INDEL";
                }
            } else {
                if (length(ref) > length(alt) && length(alt) == 1) {
                    type = "SV_DEL";
                } else if (length(alt) > 1 && length(ref) == 1) {
                    type = "SV_INS";
                } else {
                    type = "SV_COMPLEX";
                }
            }

            print id, type;
        }
    ' > vartype.tsv

# ------------------------------------------------------------------------------
# 10. Merge summary tables
# ------------------------------------------------------------------------------

log "Merging summary tables"
python3 "${SCRIPT_DIR}/merge-tables.py" \
    "genotyping-statistics_${PREFIX}.tsv" \
    "mendelian-statistics_${PREFIX}.tsv" \
    "self_${PREFIX}_variant-stats.tsv" \
    "${BUBBLE_STATS}" \
    vartype.tsv \
    "summary_${PREFIX}.tsv"

log "Done."
log "Main outputs:"
log "  trio-statistics_${PREFIX}.tsv"
log "  mendelian-statistics_${PREFIX}.tsv"
log "  genotyping-statistics_${PREFIX}.tsv"
log "  self_${PREFIX}_variant-stats.tsv"
log "  vartype.tsv"
log "  summary_${PREFIX}.tsv"
