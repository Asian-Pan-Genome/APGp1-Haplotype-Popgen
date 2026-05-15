# Variant Calling

This directory contains workflows and helper scripts for short-read and pangenome assembly-based variant discovery, filtering, and processing.

---

## Directory Overview

```text
01_variant_calling/
├── 01_short_read/
├── 02_pangenome_assembly/
```

---

# 01_short_read

Short-read variant calling workflow based on reference-guided alignment and GATK best practices.

## Structure

```text
01_short_read/
├── config/        # Reference-specific configuration files
├── resources/     # External resources and known variant databases
├── reference/     # Reference genome files and indexes
├── scripts/       # Main analysis scripts
└── README.md
```

## Configurations

```text
config/chm13.config    # Parameters for CHM13 reference
config/hg38.config     # Parameters for GRCh38 reference
```

## Scripts

```text
scripts/00.prepare_reference.sh    # Prepare and index reference genome
scripts/01.mapping.sh              # Read alignment
scripts/02.postprocess_bam.sh      # BAM sorting and duplicate processing
scripts/03.bqsr_haplotypecaller.sh # BQSR and GATK HaplotypeCaller
scripts/04.combine_gvcf.sh         # Merge per-sample GVCFs
scripts/05.genotype_gvcf.sh        # Joint genotyping
scripts/06.vqsr_snp.sh             # SNP VQSR filtering
scripts/07.vqsr_indel.sh           # INDEL VQSR filtering
```

---

# 02_pangenome_assembly

Assembly/pangenome-based variant processing and downstream filtering workflows.

## Structure

```text
02_pangenome_assembly/
├── 01_primary/                         # Primary variant processing workflow
├── 02_vcfwave_full_decomposition/      # Variant decomposition/normalization
├── 01.variant_filter.sh                # Variant filtering
├── subset_unrelated.sh                 # Extract unrelated samples
├── variant_stats.sh                    # Variant statistics summary
├── hap2dip.py                          # Convert haplotype calls to diploid format
└── README.md
```

---
