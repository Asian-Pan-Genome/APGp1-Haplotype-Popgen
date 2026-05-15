# APGp1 Short-Read Variant Calling Pipeline

This repository contains the short-read germline variant calling workflow used in the APGp1 study, supporting both :contentReference[oaicite:0]{index=0} (CRGh38) and :contentReference[oaicite:1]{index=1} (T2T-CHM13v2.0+HG002Y) references.

---

## Overview

For 160 APGp1 samples, paired-end 150 bp short reads were quality filtered, uniformly downsampled to ~33× coverage, aligned to both reference assemblies, and processed following :contentReference[oaicite:2]{index=2} Best Practices. Small variants were jointly genotyped, recalibrated using VQSR, and further stratified into short-read accessible and inaccessible genomic regions.

---

## Repository Structure

```bash
project/
├── config/
├── reference/
├── resources/
├── scripts/
│   ├── 00.prepare_reference.sh
│   ├── 01.mapping.sh
│   ├── 02.postprocess_bam.sh
│   ├── 03.bqsr_haplotypecaller.sh
│   ├── 04.combine_gvcf.sh
│   ├── 05.genotype_gvcf.sh
│   ├── 06.vqsr_snp.sh
│   └── 07.vqsr_indel.sh
├── bam/
├── gvcf/
├── joint_calling/
└── vqsr/
```

---

## Software

- :contentReference[oaicite:3]{index=3} v0.23.4
- :contentReference[oaicite:4]{index=4}
- :contentReference[oaicite:5]{index=5}
- :contentReference[oaicite:6]{index=6}
- :contentReference[oaicite:7]{index=7} v4.1.8.1
- :contentReference[oaicite:8]{index=8} v1.16

---

## Workflow

### 1. Read preprocessing

Paired-end short reads were processed using :contentReference[oaicite:9]{index=9} for adapter trimming and quality filtering (`-u 30 -q 20`), followed by uniform downsampling to ~33× coverage.

---

### 2. Reference preparation

Reference indexing and sequence dictionary generation were performed for both CRGh38 and T2T-CHM13v2.0+HG002Y.

Script:

- `scripts/00.prepare_reference.sh`

---

### 3. Read mapping and BAM preprocessing

Reads were aligned using :contentReference[oaicite:10]{index=10} and filtered to remove low-quality, secondary, and supplementary alignments.

Duplicate removal, coordinate sorting, and BAM indexing were subsequently performed.

Scripts:

- `scripts/01.mapping.sh`
- `scripts/02.postprocess_bam.sh`

---

### 4. BQSR and per-sample variant calling

Base Quality Score Recalibration (BQSR) and per-sample GVCF generation were performed using GATK `BaseRecalibrator`, `ApplyBQSR`, and `HaplotypeCaller`.

Script:

- `scripts/03.bqsr_haplotypecaller.sh`

---

### 5. Joint genotyping

Per-sample GVCFs were merged and jointly genotyped across all samples.

Scripts:

- `scripts/04.combine_gvcf.sh`
- `scripts/05.genotype_gvcf.sh`

Merged VCFs were additionally concatenated using :contentReference[oaicite:11]{index=11} where appropriate.

---

### 6. Variant Quality Score Recalibration (VQSR)

Variant recalibration was performed separately for SNPs and INDELs using lifted public truth resources compatible with T2T-CHM13.

Scripts:

- `scripts/06.vqsr_snp.sh`
- `scripts/07.vqsr_indel.sh`

Public CHM13 VQSR resources were downloaded from:

- [T2T CHM13 GATK Resource Bundle](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T%2FCHM13%2Fassemblies%2Fvariants%2FGATK_CHM13v2.0_Resource_Bundle%2F&utm_source=chatgpt.com)

Including:

- dbSNP
- HapMap
- Omni
- 1000G Phase1
- known_indels
- Mills_and_1000G_gold_standard.indels

---

## Variant Filtering

For variants in autosomes and PAR1/2, the following filters were applied:

- genotype missingness < 5%
- minor allele count (MAC) ≥ 1
- Hardy-Weinberg equilibrium P-value > 1 × 10⁻¹⁰
- `FILTER == PASS`
- ALT allele ≠ `"*"`
- Mendelian error rate < 5%
- `VQSLOD > 0`

For non-PAR chrX and chrY, variants were retained in haploid representation and HWE filtering was not applied due to sex imbalance and reduced sample size.

---

## Allele Frequency Comparison

Assembly-derived and short-read-derived variants on autosomes were intersected using `bcftools isec -c none`, and allele frequencies were compared across shared variants.

---

## Variant Masks and Genome Stratification

Telomeric and centromeric regions were masked using:

- [CHM13 annotations](https://github.com/marbl/CHM13?utm_source=chatgpt.com)

Genome-wide stratification into short-read accessible and inaccessible regions was performed using the lenient Panmask file:

- [Panmask pm151a-v2](https://zenodo.org/records/15683328?utm_source=chatgpt.com)

For selection analyses, a more stringent mask (`pm151b-v3`) was additionally applied to exclude low-complexity regions identified by SDUST.

Additional :contentReference[oaicite:15]{index=15} stratifications included:

- [GIAB CHM13 difficult regions v3.6](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/Union/CHM13_alldifficultregions.bed.gz?utm_source=chatgpt.com)
- [CHM13 short-read accessibility mask](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/combined_mask.bed.gz?utm_source=chatgpt.com)

---

## Citation

If you use this workflow, please cite:

- :contentReference[oaicite:18]{index=18}
- :contentReference[oaicite:19]{index=19}
- :contentReference[oaicite:20]{index=20}
- :contentReference[oaicite:21]{index=21}
