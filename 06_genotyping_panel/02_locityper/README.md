**Some scripts were borrowed from [Locityper](https://github.com/tprodanov/locityper/tree/main/extra).**
# Requirements
- jellyfish
- locityper

# Preparing panel
```shell
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 --threads 32 --size 3G --output counts.jf CHM13v2.fasta
locityper add -d db -v input.multi.vcf.gz -r CHM13v2.fasta -j counts.jf -L chm13v2.0_RefSeq_Liftoff_v5.1.CMRG.bed --threads 32 -e 300k
```

# Evaluating panel haplotype availability
```shell
python locityper_haplotype_availability_pipeline.py \
  --loci-bed chm13v2.0_RefSeq_Liftoff_v5.1.CMRG.bed \
  --samples samples.list \
  --genotypes hprc_hgsvc.genotypes.txt \
  --db-dir db \
  --output haplotype_availability.summary.tsv \
  --raw-output haplotype_availability.raw.tsv \
  --gt-dist gt_dist.new.py \
  --num-processes 32
```

# Target genotyping
```shell
locityper preproc -i $sample.R1.fq.gz $sample.R2.fq.gz -r CHM13v2.fasta -j counts.jf -o bg --threads 16
locityper genotype -i $sample.R1.fq.gz $sample.R2.fq.gz -d db -p bg -o gts --threads 16
```

# Evaluating genotyping accuracy
```shell
python into_csv.py -I <(cat APGp1.list | while read -r i; do echo -e "$i/gts/\t$i"; done) -o summary.csv
python eval_accuracy.py -i summary.csv -d db/ -a db/loci/{}/haplotypes.paf.gz -o eval.csv
```

# Classifing 268 CMRG locus
After evaluating genotyping accuracy, you should have two files corresponding to `APGp1+` and `APGp1-` panels.
```shell
python locityper_locus_performance.py \
  --panel-plus-eval eval.ALL_260.all.csv \
  --panel-minus-eval eval.ALL_260.sub.csv \
  --sample-pop samples.id.population \
  --outdir locityper_performance \
  --prefix CMRG \
  --panel-minus-use-call-prefixes C,K \
  --make-plots \
  --make-qv-stack-plot
```
