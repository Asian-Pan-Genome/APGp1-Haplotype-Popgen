#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Calculate population-specific SV frequencies
#
# This script:
#   1. Indexes phased SV VCF files
#   2. Extracts overlapping SVs using Truvari
#   3. Matches SV coordinates between datasets
#   4. Generates population-specific sample lists
#   5. Calculates SV allele frequencies across 26 populations
#
# Requirements:
#   - bcftools
#   - truvari
#   - awk
#   - grep
#
# Usage:
#   bash get_26_pop_sv_frqs_1kgp_ont.sh \
#       511.candidate.sv.chm13.vcf.gz \
#       shapeit5-phased-callset_final-vcf.phased.vcf.gz \
#       26_pop.txt \
#       1000G_3202_samples_pop.txt
#
###############################################################################

BASE_VCF=$1
PHASED_VCF=$2
POP_LIST=$3
POP_METADATA=$4

###############################################################################
# Step 1. Index phased VCF
###############################################################################

echo "[INFO] Indexing phased VCF..."

bcftools index -t ${PHASED_VCF}

###############################################################################
# Step 2. Extract sample list
###############################################################################

echo "[INFO] Extracting phased sample list..."

bcftools query -l ${PHASED_VCF} > 908.1kgp.sample.list

###############################################################################
# Step 3. Generate population-specific sample lists
###############################################################################

echo "[INFO] Generating population-specific sample lists..."

while read -r pop; do

    grep -w -f 908.1kgp.sample.list \
        <(awk -v pop="${pop}" '$2 == pop' ${POP_METADATA}) \
        > ${pop}.sample.txt

    sample_n=$(wc -l < ${pop}.sample.txt)

    echo -e "${pop}\t${sample_n}" \
        > ${pop}.sample.count.txt

done < ${POP_LIST}

###############################################################################
# Step 4. Benchmark SVs using Truvari
###############################################################################

echo "[INFO] Running Truvari benchmark..."

truvari bench \
    --base ${BASE_VCF} \
    --comp ${PHASED_VCF} \
    --pctovl 0.5 \
    --pctsize 0.7 \
    --refdist 1000 \
    -o truvari_out

###############################################################################
# Step 5. Generate paired SV coordinates
###############################################################################

echo "[INFO] Generating paired SV coordinates..."

paste \
    <(bcftools view -H truvari_out/tp-base.vcf.gz | cut -f1-3) \
    <(bcftools view -H truvari_out/tp-comp.vcf.gz | cut -f1-3) | \
sort -k1,2V > tp.paired.sv.tsv

###############################################################################
# Step 6. Calculate population-specific SV frequencies
###############################################################################

echo "[INFO] Calculating SV allele frequencies..."

mkdir -p frq

cut -f4,5 tp.paired.sv.tsv | \
while read -r chr pos; do

    out="frq/${chr}_${pos}.frq.tsv"

    echo -e "POP\tAC\tAN\tAF" > ${out}

    while read -r pop; do

        bcftools view \
            -r "${chr}:${pos}" \
            -S ${pop}.sample.txt \
            ${PHASED_VCF} | \
        bcftools +fill-tags \
            --threads 1 \
            - -- \
            -t AN,AC,AF,MAC:1=MAC,HWE | \
        bcftools query \
            -f "%AC\t%AN\t%AF\n" | \
        awk -v pop="${pop}" 'BEGIN{OFS="\t"}{
            print pop, $0
        }'

    done < ${POP_LIST} >> ${out}

done

###############################################################################
# Done
###############################################################################

echo "[INFO] Population frequency calculation completed successfully."
