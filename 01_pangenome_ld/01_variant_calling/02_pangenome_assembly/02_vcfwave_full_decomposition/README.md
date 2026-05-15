# Pangenome Variant Processing Pipeline

This repository contains scripts used for pangenome-based variant processing, filtering, and summary statistics generation in the Asian Pan-Genome Project (APGp1).

## Overview

A total of 539 haplotype-resolved assemblies were used throughout the study, including:

- T2T-CHM13v2.0
- T2T-CN1v1.0
- GRCh38
- 536 haplotypes from 268 diploid samples

Global pangenome graphs were constructed using Minigraph-Cactus with T2T-CHM13v2.0 as the reference backbone.

The downstream workflow includes:

1. Graph-based variant decomposition
2. Variant normalization and deduplication
3. Diploid VCF conversion
4. Sample filtering
5. Biallelic variant extraction
6. INFO tag recalculation
7. Population-level variant statistics

---

# Repository Structure

```text
.
├── scripts/
│   ├── sv_filter_pipeline.sh
│   ├── subset_unrelated.sh
│   └── variant_stats.sh
│
├── sample_lists/
│   ├── 536_hap.txt
│   ├── 258.unrelated.sample.txt
│   └── 99.unrelated.CHA.sample.txt
│
├── resources/
│   ├── chm13.pm151a-v2.easy.bed.gz
│   ├── CHM13_notinalldifficultregions.bed.gz
│   └── chm13v2.0.fa
│
└── README.md
```

---

# Requirements

## Software

- bcftools >= 1.16
- Python >= 3.8
- Minigraph-Cactus >= 2.7.2
- Minimap2 >= 2.26
- vcfwave >= 1.0.13

## External workflows

This repository relies on the following published workflows:

- VCF preparation pipeline  
  https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC

- SV collapsing workflow  
  https://github.com/Han-Cao/collapse-bubble/tree/pangenie-pipeline

- SVAN annotation  
  https://github.com/REPBIO-LAB/SVAN

---

# Workflow

## 1. Pangenome graph construction

Global pangenome graphs were generated using Minigraph-Cactus with T2T-CHM13v2.0 as the backbone reference.

Unanchored contigs and mitochondrial assemblies were excluded.

For PAR regions, fragmented assemblies were rescued using Minimap2 alignments against CHM13v2.0 PAR intervals.

---

## 2. Variant calling and decomposition

Variants were extracted from the pangenome graph using:

```bash
vg deconstruct
```

Large and nested indels (>100 kb) were removed using:

```bash
vcfbub
```

Subsequent processing included:

- decomposition of nested alleles
- full normalization using vcfwave
- SV collapsing and deduplication
- genotype-aware variant merging
- phased VCF conversion

---

## 3. Variant filtering

Final filtering steps implemented in this repository include:

- telomere/centromere masking
- diploid genotype conversion
- normalization and sorting
- biallelic extraction
- missingness filtering
- HWE filtering
- allele frequency annotation

Example:

```bash
bash scripts/sv_filter_pipeline.sh SVs chr1
```

---

## 4. Subset unrelated samples

Extract unrelated sample subsets and recalculate INFO statistics.

### All unrelated samples

```bash
bash scripts/subset_unrelated.sh SVs chr1 unrelated
```

### CHA unrelated subset

```bash
bash scripts/subset_unrelated.sh SVs chr1 unrelated_CHA
```

---

## 5. Variant statistics

Generate chromosome-level variant summary statistics.

### MAF-based statistics

```bash
bash scripts/variant_stats.sh maf chr1
```

### AF-based statistics

```bash
bash scripts/variant_stats.sh af chr1
```

### Primary site extraction

```bash
bash scripts/variant_stats.sh primary SVs
```

---

# Variant Stratification

Genome-wide stratifications were applied using:

## Panmask easy regions

- pm151a-v2 (lenient SR-accessible regions)

Used for genome-wide easy/difficult region analyses.

## GIAB stratifications

- CHM13 difficult-region mask
- CHM13 short-read accessibility mask

Used for supplementary benchmarking analyses.

---

# Output Files

Main outputs include:

| File | Description |
|---|---|
| `*.biallelic.vcf.gz` | Biallelic variant callset |
| `*.addtag.vcf.gz` | Variants with recalculated INFO tags |
| `*.filter.vcf.gz` | Final filtered callset |
| `*.unrelated.vcf.gz` | Unrelated sample subset |
| `*.stat.MAF.txt` | MAF-based statistics |
| `*.stat.AF.txt` | AF-based statistics |

---

# Citation

If you use this repository, please cite:

> Asian Pan-Genome Project (APGp1) study manuscript.

---

# Notes

- All analyses were performed on the T2T-CHM13v2.0 coordinate system.
- Telomeric and centromeric regions were masked prior to downstream analyses.
- Structural variants were sequence-resolved and fully normalized before filtering.
