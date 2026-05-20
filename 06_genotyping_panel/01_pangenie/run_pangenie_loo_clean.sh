#!/usr/bin/env bash
# ==============================================================================
# PanGenie leave-one-out genotyping and concordance evaluation
#
# Usage:
#   bash run_pangenie_loo_clean.sh sample.list
#
# The only required command-line input is a one-column sample list.
# Edit the variables below if your input files or references are in different paths.
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# User-editable paths
# ------------------------------------------------------------------------------

# Directory containing input.multi.vcf, input.bi.vcf.gz, biallelic-bubbles.bed,
# and complex-bubbles.bed. By default, use the current directory.
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

# Output directory. One subdirectory will be created for each sample.
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}}"

# Input VCFs. These replace the original long file names.
INPUT_MULTI_VCF="${INPUT_MULTI_VCF:-${PROJECT_DIR}/input.multi.vcf}"
INPUT_BI_VCF_GZ="${INPUT_BI_VCF_GZ:-${PROJECT_DIR}/input.bi.vcf.gz}"

# Bubble annotation BED files.
BIALLELIC_BUBBLES_BED="${BIALLELIC_BUBBLES_BED:-${PROJECT_DIR}/biallelic-bubbles.bed}"
COMPLEX_BUBBLES_BED="${COMPLEX_BUBBLES_BED:-${PROJECT_DIR}/complex-bubbles.bed}"

# Reference files.
CHM13_FASTA="${CHM13_FASTA:-/share/home/zhanglab/user/chenquanyu/rawdata/CHM13/ref/CHM13v2.fasta}"
CRAM_REFERENCE_FASTA="${CRAM_REFERENCE_FASTA:-/share/home/zhanglab/user/chenquanyu/HumanGenomeArchive/GSA/1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# Sample CRAM files are assumed to be: ${CRAM_DIR}/${sample}${CRAM_SUFFIX}
CRAM_DIR="${CRAM_DIR:-/share/home/zhanglab/user/chenquanyu/HumanGenomeArchive/GSA/1KGP/Related}"
CRAM_SUFFIX="${CRAM_SUFFIX:-.final.cram}"

# PanGenie singularity image.
PANGENIE_SIF="${PANGENIE_SIF:-/share/home/zhanglab/user/chenquanyu/softwares/pangenie.sif}"

# Helper Python scripts. Put the required *.py scripts in the same directory as this
# bash file by default, or set HELPER_SCRIPT_DIR to another directory.
PIPELINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_SCRIPT_DIR="${HELPER_SCRIPT_DIR:-${PIPELINE_SCRIPT_DIR}}"

# Runtime settings.
THREADS="${THREADS:-16}"
BGZIP_THREADS="${BGZIP_THREADS:-4}"
LOAD_SINGULARITY_MODULE="${LOAD_SINGULARITY_MODULE:-1}"
SINGULARITY_MODULE="${SINGULARITY_MODULE:-singularity/3.5.3}"

# Extra bind paths for Singularity, for example: /share:/share,/data:/data
SINGULARITY_EXTRA_BIND="${SINGULARITY_EXTRA_BIND:-}"

# Variant and bubble classes to evaluate.
VARIANT_TYPES=(SV_DEL SV_INS SV_COMPLEX)
BUBBLE_TYPES=(biallelic-bubbles complex-bubbles)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

usage() {
    cat >&2 <<USAGE
Usage:
  bash $(basename "$0") sample.list

Required input:
  sample.list: one sample ID per line. Empty lines and lines beginning with # are ignored.

Main editable variables:
  PROJECT_DIR, OUT_DIR, INPUT_MULTI_VCF, INPUT_BI_VCF_GZ,
  BIALLELIC_BUBBLES_BED, COMPLEX_BUBBLES_BED,
  CHM13_FASTA, CRAM_REFERENCE_FASTA, CRAM_DIR, CRAM_SUFFIX,
  PANGENIE_SIF, HELPER_SCRIPT_DIR, THREADS.
USAGE
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local file="$1"
    [[ -f "${file}" ]] || die "File not found: ${file}"
}

require_command() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || die "Command not found: ${cmd}"
}

helper_script() {
    local script_name="$1"
    echo "${HELPER_SCRIPT_DIR}/${script_name}"
}

check_inputs() {
    require_file "${INPUT_MULTI_VCF}"
    require_file "${INPUT_BI_VCF_GZ}"
    require_file "${BIALLELIC_BUBBLES_BED}"
    require_file "${COMPLEX_BUBBLES_BED}"
    require_file "${CHM13_FASTA}"
    require_file "${CRAM_REFERENCE_FASTA}"
    require_file "${PANGENIE_SIF}"

    require_file "$(helper_script remove-missing.py)"
    require_file "$(helper_script convert-to-biallelic.py)"
    require_file "$(helper_script untypable-ids-single.py)"
    require_file "$(helper_script untypable-ids-general.py)"
    require_file "$(helper_script skip-untypable.py)"
    require_file "$(helper_script extract-varianttype.py)"
    require_file "$(helper_script get_ids.py)"
    require_file "$(helper_script genotype-evaluation.py)"

    require_command python3
    require_command bcftools
    require_command bgzip
    require_command tabix
    require_command samtools
    require_command bedtools
    require_command singularity
}

singularity_bind_args() {
    local sample_dir="$1"
    local bind_paths="${sample_dir}:${sample_dir}"
    if [[ -n "${SINGULARITY_EXTRA_BIND}" ]]; then
        bind_paths="${bind_paths},${SINGULARITY_EXTRA_BIND}"
    fi
    echo "${bind_paths}"
}

run_one_sample() {
    local sample="$1"
    local sample_dir="${OUT_DIR}/${sample}"
    local cram="${CRAM_DIR}/${sample}${CRAM_SUFFIX}"
    local bind_paths

    require_file "${cram}"

    mkdir -p "${sample_dir}"
    cd "${sample_dir}"
    bind_paths="$(singularity_bind_args "${sample_dir}")"

    echo "[$(date '+%F %T')] Start sample: ${sample}"

    cat "${INPUT_MULTI_VCF}" \
        | python3 "$(helper_script remove-missing.py)" "${sample}" \
        > "${sample}_multi_no-missing.vcf"

    zcat "${INPUT_BI_VCF_GZ}" \
        | python3 "$(helper_script remove-missing.py)" "${sample}" \
        > "${sample}_bi_no-missing.vcf"

    bcftools view --samples "^${sample}" "${INPUT_MULTI_VCF}" \
        | bcftools view --min-ac 1 -o "panel-${sample}.vcf"

    bcftools view --samples "${sample}" "${sample}_bi_no-missing.vcf" \
        -o "truth-${sample}.vcf"

    singularity run -B "${bind_paths}" "${PANGENIE_SIF}" \
        PanGenie-index \
        -v "panel-${sample}.vcf" \
        -r "${CHM13_FASTA}" \
        -t "${THREADS}" \
        -o index

    singularity run -B "${bind_paths}" "${PANGENIE_SIF}" \
        PanGenie \
        -f "${sample_dir}/index" \
        -i <(samtools view \
                -T "${CRAM_REFERENCE_FASTA}" \
                "${cram}" \
                -Sb \
                -@ "${THREADS}" \
             | samtools fastq -@ "${THREADS}") \
        -s "${sample}" \
        -o pangenie \
        -j "${THREADS}" \
        -t "${THREADS}"

    bgzip -f -@ "${THREADS}" pangenie_genotyping.vcf

    zcat pangenie_genotyping.vcf.gz \
        | python3 "$(helper_script convert-to-biallelic.py)" "${INPUT_BI_VCF_GZ}" "${CHM13_FASTA}" \
        > pangenie_genotyping-biallelic.unsorted.vcf

    {
        grep '^#' pangenie_genotyping-biallelic.unsorted.vcf
        grep -v '^#' pangenie_genotyping-biallelic.unsorted.vcf | sort -k1,1 -k2,2n
    } > pangenie_genotyping-biallelic.vcf
    rm -f pangenie_genotyping-biallelic.unsorted.vcf

    zcat "${INPUT_BI_VCF_GZ}" \
        | python3 "$(helper_script untypable-ids-single.py)" . "${sample}"

    python3 "$(helper_script untypable-ids-general.py)" "${INPUT_BI_VCF_GZ}" "${INPUT_MULTI_VCF}" \
        | cat - "${sample}-untypable.tsv" \
        | sort \
        | uniq \
        > "${sample}-untypable-all.tsv"

    for variant_type in "${VARIANT_TYPES[@]}"; do
        cat pangenie_genotyping-biallelic.vcf \
            | python3 "$(helper_script skip-untypable.py)" "${sample}-untypable-all.tsv" \
            | python3 "$(helper_script extract-varianttype.py)" "${variant_type}" \
            | bgzip -@ "${BGZIP_THREADS}" -c \
            > "pangenie_genotyping-biallelic-typable-${variant_type}.vcf.gz"
        tabix -f -p vcf "pangenie_genotyping-biallelic-typable-${variant_type}.vcf.gz"

        cat "truth-${sample}.vcf" \
            | python3 "$(helper_script skip-untypable.py)" "${sample}-untypable-all.tsv" \
            | python3 "$(helper_script extract-varianttype.py)" "${variant_type}" \
            | bgzip -@ "${BGZIP_THREADS}" -c \
            > "truth-${sample}-typable-${variant_type}.vcf.gz"
        tabix -f -p vcf "truth-${sample}-typable-${variant_type}.vcf.gz"

        for bubble_type in "${BUBBLE_TYPES[@]}"; do
            local bubble_bed
            if [[ "${bubble_type}" == "biallelic-bubbles" ]]; then
                bubble_bed="${BIALLELIC_BUBBLES_BED}"
            else
                bubble_bed="${COMPLEX_BUBBLES_BED}"
            fi

            mkdir -p "concordance/${bubble_type}_${variant_type}"

            cat "${sample}_bi_no-missing.vcf" \
                | python3 "$(helper_script skip-untypable.py)" "${sample}-untypable-all.tsv" \
                | python3 "$(helper_script extract-varianttype.py)" "${variant_type}" \
                | bedtools intersect -header -a - -b "${bubble_bed}" -u -f 0.5 \
                | python3 "$(helper_script get_ids.py)" \
                > "${bubble_type}_${variant_type}.tsv"

            bedtools intersect \
                -header \
                -a "truth-${sample}-typable-${variant_type}.vcf.gz" \
                -b "${bubble_bed}" \
                -u \
                -f 0.5 \
                | bgzip -@ "${BGZIP_THREADS}" -c \
                > "concordance/${bubble_type}_${variant_type}_base.vcf.gz"
            tabix -f -p vcf "concordance/${bubble_type}_${variant_type}_base.vcf.gz"

            bedtools intersect \
                -header \
                -a "pangenie_genotyping-biallelic-typable-${variant_type}.vcf.gz" \
                -b "${bubble_bed}" \
                -u \
                -f 0.5 \
                | bgzip -@ "${BGZIP_THREADS}" -c \
                > "concordance/${bubble_type}_${variant_type}_call.vcf.gz"
            tabix -f -p vcf "concordance/${bubble_type}_${variant_type}_call.vcf.gz"

            python3 "$(helper_script genotype-evaluation.py)" \
                "concordance/${bubble_type}_${variant_type}_base.vcf.gz" \
                "concordance/${bubble_type}_${variant_type}_call.vcf.gz" \
                "${bubble_type}_${variant_type}.tsv" \
                --qual 0 \
                > "concordance/${bubble_type}_${variant_type}/summary.txt"
        done
    done

    echo "[$(date '+%F %T')] Finished sample: ${sample}"
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    usage
    exit 1
fi

SAMPLE_LIST="$1"
require_file "${SAMPLE_LIST}"
mkdir -p "${OUT_DIR}"

if [[ "${LOAD_SINGULARITY_MODULE}" == "1" ]]; then
    if type module >/dev/null 2>&1; then
        module load "${SINGULARITY_MODULE}"
    else
        echo "WARNING: LOAD_SINGULARITY_MODULE=1 but the module command is unavailable; continuing without module load." >&2
    fi
fi

check_inputs

while IFS= read -r sample || [[ -n "${sample}" ]]; do
    sample="${sample%%[[:space:]]*}"
    [[ -z "${sample}" ]] && continue
    [[ "${sample}" =~ ^# ]] && continue
    run_one_sample "${sample}"
done < "${SAMPLE_LIST}"

