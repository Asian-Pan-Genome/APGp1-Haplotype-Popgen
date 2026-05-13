#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Summarize out-of-panel phasing switch and flip error statistics
#
# This script:
#   1. Aggregates chromosome-level WhatsHap statistics
#   2. Computes autosome-wide switch and flip error rates
#   3. Summarizes out-of-panel phasing accuracy across samples
#
# Requirements:
#   - awk
#   - sed
#   - sort
#
# Usage:
#   bash 02.summarize_switch_flip_error_statistics.sh
#
###############################################################################

SAMPLE_LIST=160.sample.child-01.list

METHOD=assem

OUTDIR=1.stat

###############################################################################
# Start
###############################################################################

date

mkdir -p ${OUTDIR}

###############################################################################
# Step 1. Generate chromosome-level statistics
###############################################################################

echo "[INFO] Generating chromosome-level statistics..."

for chr in chr{1..22}; do

    mkdir -p ${OUTDIR}/${chr}

    {
        echo -e "sample\tall_assessed_pairs\tall_switches\tall_switch_rate\tall_true_switchs\tall_flips\tall_switchflip_rate\tall_true_switch_rate\tall_flip_rate"

        while read -r sample; do

            cut -f9-13 \
                ${chr}/${sample}.${METHOD}.1kgp.out_of_panel.pairwise.tsv | \
            sed '1d;s#/#\t#g' | \
            awk -v sample="${sample}" '
                BEGIN{OFS="\t"}
                {
                    print sample,
                          $0,
                          $(NF-2)/$1,
                          $(NF-1)/$1
                }
            '

        done < ${SAMPLE_LIST}

    } > ${OUTDIR}/${chr}/${chr}.${METHOD}.stat.tsv

done

###############################################################################
# Step 2. Generate autosome-wide statistics
###############################################################################

echo "[INFO] Calculating autosome-wide statistics..."

{
    echo -e "sample\tall_assessed_pairs\tall_switches\tall_switch_rate\tall_switchflips\tall_true_switchs\tall_flips\tall_switchflip_rate\tall_true_switch_rate\tall_flip_rate"

    while read -r sample; do

        for chr in chr{1..22}; do

            sed '1d' \
                ${chr}/${sample}.${METHOD}.1kgp.out_of_panel.pairwise.tsv | \
            cut -f9-13 | \
            sed 's#/#\t#g'

        done | \
        awk -v sample="${sample}" '
            BEGIN{OFS="\t"}

            {
                assessed += $1
                switches += $2
                switchflips += ($4 + $5)
                true_switches += $4
                flips += $5
            }

            END{
                print sample,
                      assessed,
                      switches,
                      switches/assessed,
                      switchflips,
                      true_switches,
                      flips,
                      switchflips/assessed,
                      true_switches/assessed,
                      flips/assessed
            }
        '

    done < ${SAMPLE_LIST}

} > ${OUTDIR}/autosome.${METHOD}.stat.tsv

###############################################################################
# Done
###############################################################################

echo "[INFO] Out-of-panel switch error summary completed successfully."

date
