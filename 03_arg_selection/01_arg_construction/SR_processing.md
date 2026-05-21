# Steps taken for short-read data

Some of these steps may be redundant with previous steps. Each step is run for each autosome.

1. Convert BCF to VCF and remove non-segregating sites

```
bcftools view --min-af 0.0000001:nref --max-af 0.9999999:nref ${INDIR}/${IN}.bcf -o ${OUTDIR}/${IN}.downsamp.vcf.gz -O z
```

2. Remove variants with missingness more than 10%

```
bcftools filter -e "F_MISSING > 0.1" ${INDIR}/${IN}.vcf.gz -o ${OUTDIR}/${IN}.m10p.vcf.gz -O z
```

3. Only SNPs; drop some meta-data

```
bcftools view -v snps ${INDIR}/${IN}.vcf.gz | bcftools annotate -x FORMAT --remove INFO - -o ${OUTDIR}/${IN}.snvs.vcf.gz -O z
```

4. Convert to IGD

```
igdtools ${INDIR}/${IN}.vcf.gz -o ${OUTDIR}/${IN}.igd
```

5. Drop sites with non-SNVs (bcftools doesn't always do this completely...?)

```
igdtools --drop-non-snv-sites ${INDIR}/${IN}.may_have_nonsnv.igd -o ${OUTDIR}/${IN}.igd
```

6. Infer ARGs with tsinfer

```
OUTPUT_DIR=${OUTDIR} DATA_DIR=working/data_dir python make_args.py ${INDIR}/${IN} tsinfer chm13 -j 20 -c ${CHROM}
```

7. Collect statistics

```
OUTPUT_DIR=${OUTDIR} DATA_DIR=working/data_dir python ts_stats.py chr${CHROM} ${TREE} > ${OUTDIR}/cha.chr${CHROM}.stats.json
```
