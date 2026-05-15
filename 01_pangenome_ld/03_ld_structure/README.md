# Linkage Disequilibrium and Haplotype Structure Analysis

This repository contains scripts used for linkage disequilibrium (LD), haplotype block, and SV-centric linkage analyses in the APGp1 pangenome study.

---

## Overview

The workflow includes:

1. Population-specific LD decay analysis
2. Haplotype block inference
3. SV-centric LD analysis
4. GWAS/QTL tagging analysis
5. LD density visualization

Analyses were performed using assembly-derived and statistically phased short-read callsets from unrelated individuals.

---

## Requirements

### Software

- PLINK v1.90b7.2
- PopLDdecay v3.43
- bcftools
- Python ≥ 3.8
- R ≥ 4.3
- Perl

### Python packages

```bash
pip install pandas numpy matplotlib seaborn
```

---

## Directory Structure

```text
.
├── 01_ld_decay
│   ├── popLDdecay.subset.sh
│   ├── plot_ld_decay.sh
│   ├── plot_ld_decay_multi_pop.sh
│   ├── plot_kd.py
│   ├── classify.py
│   └── plot.file.full.list
│
├── 02_haploblock
│   ├── 00.plink_preprocessing.sh
│   ├── 01.plink_hapblock.sh
│   └── 155.unrelated.APGp1.sample.list
│
└── 03_sv_centric_linkage
    ├── 03.cal_LD_site.sh
    ├── 04.get_orphan.linked.sv.sh
    ├── 05.sv.linkage.stat.sh
    ├── 06.orphan.sv.stat.sh
    ├── 06.non_orphan.sv.stat.sh
    ├── 07.high_linked.sv.stat.sh
    ├── 08.sv.linkage.short.sh
    ├── 09.sv.tagged.short.var.gwas.sh
    ├── 09.sv.tagged.short.var.qtl.sh
    └── 10.sv.tagged.short.var.summary.sh
```

---

# 1. LD Decay Analysis

Scripts in `01_ld_decay/` were used to calculate and visualize LD decay across populations.

## Run PopLDdecay

```bash
bash popLDdecay.subset.sh <vcf.gz> <output_dir> <sample_list> <max_dist_kb>
```

Example:

```bash
bash popLDdecay.subset.sh \
    input.vcf.gz \
    results \
    EAS.samples.txt \
    100
```

## Plot LD decay

Single chromosome:

```bash
bash plot_ld_decay.sh chr1
```

Multi-population comparison:

```bash
bash plot_ld_decay_multi_pop.sh
```

## KDE visualization

```bash
python plot_kd.py \
    --input sv.100k.ld.dis.maxr2.tsv \
    --output sv.100k.ld.dis.maxr2.pdf
```

## SV classification

```bash
python classify.py \
    --blocks autosome.block \
    --svs input.sv.txt \
    --output-prefix classified
```

---

# 2. Haplotype Block Inference

Scripts in `02_haploblock/` were used to infer haplotype blocks from SNP-only, SNP-indel, and SNP-indel-SV callsets.

## PLINK preprocessing

```bash
bash 00.plink_preprocessing.sh <chromosome> <variant_type>
```

Example:

```bash
bash 00.plink_preprocessing.sh chr1 SNPs
```

## Haplotype block calling

```bash
bash 01.plink_hapblock.sh <chromosome> <variant_type> <group>
```

Example:

```bash
bash 01.plink_hapblock.sh chr1 SNPs EAS
```

Haplotype blocks were inferred using PLINK `--blocks` following the definition of Gabriel et al. (2002).

---

# 3. SV-Centric Linkage Analysis

Scripts in `03_sv_centric_linkage/` were used for SV–short variant LD analyses.

## Calculate LD around SVs

```bash
bash 03.cal_LD_site.sh <chromosome> <window_kb>
```

Example:

```bash
bash 03.cal_LD_site.sh chr1 1000
```

## Extract orphan SVs

```bash
bash 04.get_orphan.linked.sv.sh
```

## Generate lightweight LD tables

```bash
bash 05.sv.linkage.stat.sh chr1 1000
```

## Region-stratified statistics

```bash
bash 06.orphan.sv.stat.sh
bash 06.non_orphan.sv.stat.sh
```

## Top linked short variants

```bash
bash 07.high_linked.sv.stat.sh chr1 1000
```

## Count tagged short variants

```bash
bash 08.sv.linkage.short.sh chr1 0.8
```

## GWAS/QTL intersection

GWAS:

```bash
bash 09.sv.tagged.short.var.gwas.sh 0.8
```

QTL:

```bash
bash 09.sv.tagged.short.var.qtl.sh chr1 eQTL
```

## Summarize GWAS/QTL annotations

```bash
bash 10.sv.tagged.short.var.summary.sh chr1 0.8
```

---

# Methods Summary

## LD decay

LD decay was calculated using PopLDdecay v3.43 with phased autosomal SNPs (MAF ≥ 0.01).

## Haplotype blocks

Haplotype blocks were inferred using PLINK v1.90b7.2 with:

```text
--blocks-min-maf 0.05
--blocks-max-kb 10000
--blocks 'no-pheno-req'
```

## SV-centric LD

Pairwise LD (r²) between common SVs and nearby short variants (<50 bp) was calculated within 50 kb–1 Mb windows.

## GWAS/QTL annotation

GWAS Catalog and GTEx v10 QTL sites were intersected using:

```bash
bcftools isec -c none
```

Coordinates were lifted from GRCh38 to CHM13v2.0 before analysis.

---

# Citation

If you use this repository, please cite:

```text
[Manuscript information to be added after publication]
```
