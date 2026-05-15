# Genotyping Panels and Targeted Genotyping

This directory contains scripts and workflows for graph-based SV genotyping using `PanGenie` and locus-specific targeted genotyping using `Locityper`.

## Directory Structure

```text
06_genotyping_panel/
├── 01_pangenie/    # PanGenie genotyping and evaluation
└── 02_locityper/   # Locityper targeted genotyping analyses
```

---

## Overview

### 1. PanGenie Genotyping

We used `PanGenie` (v4.2.1) to genotype:

- 449 1kGP samples
- 824 HGDP samples

using globally diverse and balanced sample sets.

Two graph panels were evaluated:

- `APGp1+` = APGp1 + HGSVC3 + HPRCy1
- `APGp1−` = HGSVC3 + HPRCy1

#### Data Preprocessing

The original multiallelic MC graph bubble VCFs were decomposed into biallelic VCFs using the preprocessing pipeline described in the pangenome construction workflow.

The APGp1− bubble VCF was downloaded from:

```text
https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.vcf.gz
```

Bubble VCFs were used to construct `PanGenie` index files.

After genotyping, multiallelic genotype outputs were converted into biallelic VCFs using:

```text
convert-to-biallelic.py
```

from:

```text
https://bitbucket.org/jana_ebler/hprc-experiments/src/master/genotyping-experiments/workflow/scripts/convert-to-biallelic.py
```

Scripts are available in:

```text
01_pangenie/
```

#### Leave-One-Out Experiments

Leave-one-out analyses were performed by:

1. Removing one sample from the panel
2. Re-genotyping the excluded sample using the remaining panel haplotypes
3. Comparing inferred genotypes against the original assembly-derived genotypes

Private alleles unique to the left-out sample were excluded from evaluation because `PanGenie` only genotypes alleles present in the panel VCF.

Performance was evaluated using weighted genotype concordance for:

- `SV-INS`
- `SV-DEL`
- `SV-Others`

in both biallelic and multiallelic bubbles.

#### Genotype Filtering and Quality Assessment

Genotypes were filtered using a support vector regression framework adapted from previous studies.

The `gq_fail` threshold was adjusted from 50 to 20 to match the sample size:

```text
n = 449 + 824
```

Allele frequencies observed in assembly-based decomposed VCFs were compared against filtered PanGenie genotypes.

Pearson correlation coefficients between assembly-derived and genotyped allele frequencies were:

| Variant Type | Pearson r |
|---|---|
| SV-DEL | 0.980 |
| SV-INS | 0.977 |
| SV-Others | 0.977 |

Observed heterozygosity was additionally compared with Hardy-Weinberg expectations across allele frequency bins.

---

### 2. Locityper Genotyping for Target Loci

We used `Locityper` (v1.0.0) to genotype:

- 160 APGp1 samples
- 100 HGSVC3 + HPRCy1 samples

across 268 challenging medically relevant genes (CMRGs).

Five of the original 273 CMRGs absent from the T2T-CHM13 reference were excluded.

#### Genotyping Workflow

`Locityper` recruits and aligns reads to locus-specific haplotypes and infers the most likely haplotype pair using:

- Read alignment
- Insert size
- Read depth profiles

MC graph bubble VCFs were used to construct locus-specific genotype panels.

Scripts are available in:

```text
02_locityper/
```

#### Haplotype Availability Evaluation

For each sample:

1. Reads were aligned to unrelated panel haplotypes
2. Sequence divergence was calculated as:
   - edit distance / alignment length
3. Divergence values were transformed into Phred-like quality values (QV)

The haplotype pair with the highest QV was retained for comparison.

#### Genotyping Accuracy Assessment

Predicted and true haplotype pairs were permuted to minimize sequence divergence.

The minimum divergence between predicted and true haplotypes was used as the genotyping accuracy metric.

Genotypes were discarded if either:

1. Weighted distance to other probable genotypes > 30
2. Unexplained reads exceeded:
   - 1000 reads
   - and 20% of total reads

All analyses were performed for:

- APGp1+
- APGp1+ leave-one-out
- APGp1−
- APGp1− leave-one-out

#### Locus-Level Categorization

CMRG loci were categorized according to leave-one-out genotyping performance differences between APGp1+ and APGp1− panels.

For ΔQV analyses:

- QV values were capped at 45

Improvement was defined as:

```text
ΔQV = QV(APGp1+) − QV(APGp1−)
```

Fold error reduction (FER) was calculated as:

```text
FER = divergence(APGp1−) / divergence(APGp1+)
```

where:

```text
FER > 1
```

indicates improved genotyping accuracy after APGp1 inclusion.

Loci were assigned into three mutually exclusive categories:

| Category | Definition |
|---|---|
| SR-infeasible | ≥5% no-calls in both leave-one-out panels |
| Ceiling-effected | median EAS ΔQV ≤ 1 |
| Panel/ancestry-dependent | median EAS ΔQV > 1 |

Final classification counts:

| Category | Number of Loci |
|---|---|
| SR-infeasible | 11 |
| Ceiling-effected | 54 |
| Panel/ancestry-dependent | 203 |

---

## Software

Main tools used in this workflow:

- `PanGenie`
- `Locityper`
- `minigraph-cactus`
- `bcftools`
- `SVR-based genotype filtering`
