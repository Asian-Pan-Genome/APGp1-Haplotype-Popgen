# Positive Selection Scan of Structural Variants

This repository contains scripts used for genome-wide positive selection analyses of polarized biallelic structural variants (SVs) based on phased assembly-derived pangenome variation datasets.

## Overview

We polarized biallelic SNPs and SVs (MAF ≥ 0.05) into ancestral and derived states and applied five complementary selection statistics:

- PBS
- nSL
- iHH12
- XP-nSL
- XP-EHH

Haplotype-based statistics were computed using [selscan](https://github.com/szpiech/selscan) and normalized using the companion `norm` utility.
Population-averaged recombination maps were projected onto variant positions using [predictGMAP](https://github.com/szpiech/predictGMAP).

We further assessed:
- overlap with public short-read SV resources [(gnomAD SV v4.1)](https://gnomad.broadinstitute.org/data#v4-structural-variants)
- overlap with 1KCP SV and eQTL resources [(Wang et al. 2026)](https://www.nature.com/articles/s41586-026-10315-y)
- alternatively represented adaptive SV candidates in CHM13V2.0 and GRCh38
- population-specific SV frequencies based on 1019-sample 1kGP ONT long-read cohort [(Schloissnig et al. 2025)](https://www.nature.com/articles/s41586-025-09290-7)
- functional relevance and associations of adaptive SV candidates in UKB WGS sv callset

---

## Repository Structure

```text
.
├── 00_SV_ancestral_inference
├── 01_selection_scan
│   ├── 01_PBS
│   ├── 02_nSL
│   ├── 03_iHH12
│   ├── 04_XP-nSL
│   └── 05_XP-EHH
└── 02_functional_relevance
    ├── 01_GnomAD
    ├── 02_1KCP
    ├── 03_1KGP_1019_ONT
    └── 04_UKB_WAS
```

---

# Selection Scan Workflow

## 1. PBS

Population branch statistics were computed using `pbscan`.

### Site-based PBS
```bash
01_selection_scan/01_PBS/01.run_pbs_selection_scan_site.sh
```

### Window-based PBS
```bash
01_selection_scan/01_PBS/02.run_pbs_selection_scan_window.sh
```

---

## 2. nSL

### Prepare population-specific VCFs
```bash
01_selection_scan/02_nSL/00.prep_input_vcf.sh
```

### Run nSL
```bash
01_selection_scan/02_nSL/01.run_nsl_selection_scan.sh
```

### Normalize nSL scores
```bash
01_selection_scan/02_nSL/02.normalize_nsl_scores.sh
```

---

## 3. iHH12

### Predict genetic map positions
```bash
01_selection_scan/03_iHH12/00.predict_genetic_map_positions.sh
```

### Run iHH12
```bash
01_selection_scan/03_iHH12/01.run_ihh12_selection_scan.sh
```

### Normalize iHH12 scores
```bash
01_selection_scan/03_iHH12/02.normalize_ihh12_scores.sh
```

---

## 4. XP-nSL

### Prepare population-specific VCFs
```bash
01_selection_scan/04_XP-nSL/00.prep_input_vcf.sh
```

### Run XP-nSL
```bash
01_selection_scan/04_XP-nSL/01.run_xpnsl_selection_scan.sh
```

### Normalize XP-nSL scores
```bash
01_selection_scan/04_XP-nSL/02.normalize_xpnsl_scores.sh
```

---

## 5. XP-EHH

### Predict genetic map positions
```bash
01_selection_scan/05_XP-EHH/00.predict_genetic_map_positions.sh
```

### Run XP-EHH
```bash
01_selection_scan/05_XP-EHH/01.run_xpehh_selection_scan.sh
```

### Normalize XP-EHH scores
```bash
01_selection_scan/05_XP-EHH/02.normalize_xpehh_scores.sh
```

---

# Functional Relevance Analyses

## 1. Benchmark against gnomAD SV

Candidate adaptive SVs were lifted from CHM13v2.0 to GRCh38 using CrossMap and benchmarked against gnomAD SV v4.1 using Truvari.

```bash
02_functional_relevance/01_GnomAD/sv_crossmap_and_benchmark_GnomAD.sh
```

---

## 2. 1KCP SV and eQTL overlap

### Benchmark against 1KCP SVs
```bash
02_functional_relevance/02_1KCP/01.benchmark_lifted_svs_against_1kcp.sh
```

### Identify alternatively represented SVs
```bash
02_functional_relevance/02_1KCP/02.get_alternatively_represented_svs.sh
```

### Extract overlapping SV-eQTLs
```bash
02_functional_relevance/02_1KCP/03.extract_overlapping_sv_eqtls.sh
```

---

## 3. Population SV frequencies

Population-specific SV frequencies in ONT-phased 1KGP samples:

```bash
02_functional_relevance/03_1KGP_1019_ONT/get_26_pop_sv_frqs_1kgp_ont.sh
```

---

## 4. UK Biobank association analyses

Genome-wide and phenome-wide association analyses were performed using UK Biobank SV resources.

```text
02_functional_relevance/04_UKB_WAS
```

---

# Main Software

- [selscan](https://github.com/szpiech/selscan)
- [predictGMAP](https://github.com/szpiech/predictGMAP)
- [CrossMap](https://crossmap.sourceforge.net)
- [Truvari](https://github.com/ACEnglish/truvari)
- [bcftools](https://samtools.github.io/bcftools)
- [bedtools](https://bedtools.readthedocs.io/en/latest)

---

# Data Resources

- [gnomAD SV v4.1](https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz)
- [1KCP SV resource](https://yanglab.westlake.edu.cn/resources/1kcp/variant/1kcp.sv.vcf.gz)
- [1KCP eQTL resource](https://yanglab.westlake.edu.cn/resources/1kcp/eqtl/merge/1kcp.eqtl.summary.tar.gz)

---

# Citation

If you use this repository, please cite the corresponding study and the original software packages referenced above.
