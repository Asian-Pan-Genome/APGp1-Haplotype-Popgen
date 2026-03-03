# APGp1-Haplotype-Popgen
This repository hosts a comprehensive data and workflow resource for **assembly-based population genetic analysis** in the **APGp1** genome assembly.

---

## 1. Overview

Main contents:

* 1. variants and phasing comparison
* 2. recombination landscape and PRDM9 genotyping
* 3. ARG reconstruction and selection scan
* 4. LD and recombination feature in genomic complex regions
* 5. pangenome-based genotyping for short-read data

---

## 2. Data Access

The `Annotation` directory contains the final, high-quality annotation data.

- Complete centromere coordinates in APGp1 : [APGp1_complete_centromere](https://github.com/Asian-Pan-Genome/)

- Satellite tracks for each phased assembly: [native_bed_format](https://github.com/Asian-Pan-Genome/)

- HOR identified by HORmon: [HOR_HORmon](https://github.com/Asian-Pan-Genome/)

- HOR identified by HiCAT: [HOR_HiCAT](https://github.com/Asian-Pan-Genome/)

- CENP-A enrichment boundaries: [native_bed_format](https://github.com/Asian-Pan-Genome/)

## 3. Detailed Workflows

We provide and highly recommend the following two pipelines for future T2T human assemblies and related popgen studies.

### Varinat calling & haplotype phasing

**Directory:** `SatelliteAnnotationWorkflow`

This workflow was developed for the comprehensive annotation of centromeric satellites (including $\alpha$, $\beta$, $\gamma$, HSat1, HSat2, and HSat3) in human assemblies. The pipeline outputs precise centromeric coordinates for each genome and provides scripts for result visualization.

* [**Access Workflow**](https://github.com/Asian-Pan-Genome/Centromere/tree/main/SatelliteAnnotationWorkflow)

### HORmining Pipeline

**Directory:** `HORmining`

HORmining is a bioinformatics pipeline designed for robust identification of Higher-Order Repeat (HOR) structures within alpha satellite DNA. It integrates two complementary computational approaches: the graph-based **HORmon** algorithm and the hierarchical tandem repeat mining (HTRM)-based **HiCAT** algorithm. 

* [**Access Workflow**](https://github.com/Asian-Pan-Genome/Centromere/tree/main/HORmining)

---

## 4. Integrated Analysis (Code Archive)

The `Analysis`
