**This pipeline was borrowed from [Ebler](https://github.com/eblerjana/hgsvc3/tree/main/experiments/genotyping), please visit this repo for associated scripts.**


# Requirments
- Pangenie
- singularity
- bcftools
- tabix
- bedtools
- samtools

# Running Pangenie
First, construct the Pangenie index file:
```shell
singularity run -B $(pwd):$(pwd) pangenie.sif PanGenie-index -v input.multi.vcf -r CHM13v2.fasta -t 32 -o index
```
Then, run pangenie individually and convert the resulted mutiallelic genotype file to a biallelic file:
```shell
singularity run -B ${home}:${home} pangenie.sif PanGenie -f index -i $sample.fq -s ${sample} -o pangenie -j 16 -t 16
bgzip -f -@16 pangenie_genotyping.vcf
tabix -f -@16 pangenie_genotyping.vcf.gz
zcat pangenie_genotyping.vcf.gz | python3 convert-to-biallelic.py input.bi.vcf CHM13v2.fasta | bcftools sort -Oz -o pangenie_genotyping_bi_all.vcf.gz
```

After completing pangenie, individual VCFs were collected and merged by bcftools.

# Leave-one-out experiments
First, get biallelic and complex regions from pangenome:
```shell
cat input.multi.vcf | python3 variant-statistics.py alleles-per-bubble.pdf 1 > complex-bubbles.bed
sort -k1,1d -k 2,2n -k 3,3n CHM13v2.fasta.fai > biallelic-bubbles.fai
sort -k1,1d -k 2,2n -k 3,3n complex-bubbles.bed > biallelic-bubbles.tmp
bedtools complement -i biallelic-bubbles.tmp -g biallelic-bubbles.fai > biallelic-bubbles.bed
```

Next, perform the LOO experiments:
```shell
bash run_pangenie_loo_clean.sh sample.list
```

# Genotype filtering and quality assessment
## Data processing
Prepare necessary files before filtering genotypes.
```shell
bash prepare_for_SVM.sh
```

## Running regression model
```shell
python3 analysis-stepwise.py -t summary_bi_all.tsv -o plot_bi_all -n 20 --regression-only
python3 analysis-stepwise.py -t plot_bi_all_regression.tsv -o plot_bi_all -n 20 --plot-only
```
