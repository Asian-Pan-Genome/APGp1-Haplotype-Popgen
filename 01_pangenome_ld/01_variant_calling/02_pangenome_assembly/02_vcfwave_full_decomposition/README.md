# Full Decomposition of Pangenome Variant

This repository contains scripts used for fully decomposing a pangenome-based variant callset.

## Overview

Phased haploid/diploid callsets can be used for the decomposition workflow.

- T2T-CHM13v2.0 (reference)

The downstream workflow includes:

1. Graph-based variant decomposition
2. Variant normalization and deduplication
3. Diploid VCF conversion
4. Sample filtering
5. Biallelic variant extraction
6. INFO tag recalculation

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

---

Subsequent processing included:

- decomposition of nested alleles
- full normalization using vcfwave
- SV collapsing and deduplication
- genotype-aware variant merging
- phased VCF conversion

---


# Citation

If you use this repository, please cite:

> Asian Pan-Genome Project (APGp1) study manuscript.

---

# Notes

- All analyses were performed on the T2T-CHM13v2.0 coordinate system.
- Telomeric and centromeric regions were masked prior to downstream analyses.
- Structural variants were sequence-resolved and fully normalized before filtering.
