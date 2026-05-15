# Variant-Phasing-Linkage Analysis Pipeline

This repository contains workflows for variant discovery, haplotype phasing, and linkage disequilibrium (LD) structure analysis using short-read and pangenome-based genomic data.

---

## Directory Overview

```text
.
├── 01_variant_calling/   # Variant discovery and filtering workflows
├── 02_phasing/           # Haplotype phasing and comparison analyses
└── 03_ld_structure/      # LD, haplotype block, and recombination analyses
```

---

# 01_variant_calling

Variant calling workflows from short-read sequencing and pangenome assemblies.

## Contents

```text
01_variant_calling/
├── 01_short_read/            # Short-read alignment and GATK-based variant calling
└── 02_pangenome_assembly/    # Assembly/pangenome-based variant processing
```

---

# 02_phasing

Haplotype phasing workflows and reference panel analyses.

## Includes

- Reference panel preparation
- Haplotype phasing
- Phasing comparison and evaluation

---

# 03_ld_structure

Analyses of linkage disequilibrium and haplotype structure across genomic regions.

## Includes

- LD calculation
- Haplotype block inference
- Structural variant effects on LD structure


---
