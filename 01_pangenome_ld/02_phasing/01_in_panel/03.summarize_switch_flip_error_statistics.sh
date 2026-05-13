#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Summarize phasing switch and flip error statistics
#
# This script:
#   1. Aggregates chromosome-level WhatsHap statistics
#   2. Computes autosome-wide switch and flip error rates
#   3. Merges NGS- and assembly-based phasing statistics
#   4. Performs paired statistical comparisons in R
#
# Requirements:
#   - awk
#   - sed
#   - paste
#   - sort
#   - Rscript
#
# Usage:
#   bash 03.summarize_switch_flip_error_statistics.sh
#
###############################################################################

SAMPLE_LIST=136.sample.child-01.list

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

    for method in ngs assem; do

        {
            echo -e "sample\tall_assessed_pairs\tall_switches\tall_switch_rate\tall_true_switchs\tall_flips\tall_switchflip_rate\tall_true_switch_rate\tall_flip_rate"

            while read -r sample; do

                cut -f9-13 \
                    ${chr}/${sample}.overlap.${method}.1kgp.pairwise.tsv | \
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

        } > ${OUTDIR}/${chr}/${chr}.${method}.stat.tsv

    done
done

###############################################################################
# Step 2. Generate autosome-wide statistics
###############################################################################

echo "[INFO] Calculating autosome-wide statistics..."

for method in ngs assem; do

    {
        echo -e "sample\tall_assessed_pairs\tall_switches\tall_switch_rate\tall_switchflips\tall_true_switchs\tall_flips\tall_switchflip_rate\tall_true_switch_rate\tall_flip_rate"

        while read -r sample; do

            for chr in chr{1..22}; do

                sed '1d' \
                    ${chr}/${sample}.overlap.${method}.1kgp.pairwise.tsv | \
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

    } > ${OUTDIR}/autosome.${method}.stat.tsv

done

###############################################################################
# Step 3. Merge NGS and assembly statistics
###############################################################################

echo "[INFO] Merging NGS and assembly statistics..."

paste \
    <(awk '{print $1"\t"$(NF-1)"\t"$NF}' ${OUTDIR}/autosome.ngs.stat.tsv) \
    <(awk '{print $1"\t"$(NF-1)"\t"$NF}' ${OUTDIR}/autosome.assem.stat.tsv) | \
awk 'BEGIN{OFS="\t"}{
    print $1, $2, $5, $3, $6
}' | \
sed '1d' | \
sed '1i sample\tngs_ser\tassem_ser\tngs_fer\tassem_fer' \
> ${OUTDIR}/autosome.merge.stat.tsv

###############################################################################
# Step 4. Sort final summary table
###############################################################################

echo "[INFO] Sorting final summary table..."

cat \
    <(sed -n '1p' ${OUTDIR}/autosome.merge.stat.tsv) \
    <(
        sed '1d' ${OUTDIR}/autosome.merge.stat.tsv | \
        sort -k2,2nr | \
        sed '1,4d' | \
        sort -k1,1V
    ) \
> ${OUTDIR}/136.autosome.merge.stat.tsv

###############################################################################
# Step 5. Statistical testing
###############################################################################

echo "[INFO] Running paired statistical tests..."

Rscript 04.136.paired.calculate.ser_fer.p_value.R

###############################################################################
# Done
###############################################################################

echo "[INFO] Switch error summary completed successfully."

date
