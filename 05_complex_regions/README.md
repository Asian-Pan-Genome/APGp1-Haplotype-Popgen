# Structurally Complex Regions and Structural Haplotypes

This directory contains workflows for graph-based stratification of structurally complex regions (SCRs), structural haplotype decomposition, linkage characterization, population differentiation analyses, and great ape comparative analyses.

## Directory Structure

```text
05_complex_regions/
├── 01_SCRs/                  # SCR identification and linkage analyses
└── 02_structural_haplotype/  # Structural haplotype construction and refinement
```

---

## Overview

### 1. Graph-Based Stratification of Structurally Complex Regions

In this study, “structural complexity” broadly refers to genomic features that limit conventional short-read population analyses, including:

- Repeat arrays
- Segmental duplications (SDs)
- Inversions
- Copy-number variants (CNVs)
- Complex structural variants (SVs)

Structurally complex regions (SCRs) were defined as autosomal loci:

- larger than 10 kb
- containing multiple co-occurring SVs

using pangenome graph topology and assembly-derived SV calls.

#### SCR Identification

SCRs were identified from 1-kb windows with 500-bp steps using graph-derived complexity annotations [(Han et al. 2026)](under review).

Windows were classified into:

- `simple`
- `CC1–CC7`

where:

- `CC1` = inter-segmental homology
- `CC2` = intra-segmental homology
- `CC3–CC7` = increasing population-scale structural diversity

SCR construction procedure:

1. Retained windows annotated as `CC2–CC7`
2. Merged nearby windows using:

```text
bedtools merge -d 10000
```

3. Adjacent regions were merged if:
   - each region length exceeded the separating distance
   - combined region occupied >80% of merged interval
4. Regions were iteratively refined until contiguous `simple` or `CC1` segments accounted for ≤40% of total length
5. Regions composed entirely of `simple`, `CC1`, or `CC2` annotations were removed

Results:

- 47.2% of Panmask SR-inaccessible regions overlapped SCRs
- 92.5% of PAV-derived complex variants overlapped SCRs

Scripts are available in:

```text
01_SCRs/
```

---

### 2. PAV Structural Variant Calling

To characterize structural diversity within SCRs, we used `PAV` (v3.0.0) with T2T-CHM13v2.0 as reference across 536 haplotypes.

Extracted variant categories included:

- `svindel_ins`
- `svindel_del`
- `sv_inv`

Simple SVs and complex SVs were merged into the final PAV-derived SV callset.

---

### 3. Inter-Flank Linkage Analyses

#### SCR Inter-Flank Linkage

For each SCR:

1. Flanking `CC2` intervals adjacent to the SCR core were removed
2. The core region was extended by ±100 kb
3. SNP LD (`r²`) was calculated using `PLINK`

To quantify inter-flank LD:

- SNPs inside SCR cores and immediate ±100 kb regions were excluded
- Maximum pairwise `r²` among remaining SNPs was recorded

SCR length was defined as the core-region length.

#### Random Easy Regions

To establish background linkage patterns:

- 200,000 random easy regions were generated
- Regions overlapping any non-simple or non-CC1 windows were removed

This yielded:

```text
52,016 easy regions
```

The same LD analyses were applied to these regions.

---

## Structural Haplotype Construction

### 4. Structural Haplotype Construction Using PGR-TK

Structural haplotypes were constructed across 2,111 autosomal SCRs using `PGR-TK` (v0.6.0).

For each SCR:

1. SCR ±5 kb flanking regions were lifted to 536 haplotypes using `halLiftover`
2. BED intervals were merged using:

```text
bedtools merge -d 10000000
```

3. Sequences were extracted
4. Haplotypes containing assembly gaps were excluded

#### Repeat Masking

Simple repeats were annotated using `TRF` with parameters:

```text
2 7 7 80 10 50 2000 -d -h
```

Repeats satisfying:

- ≥80% identity
- >600 bp length

were masked before decomposition.

#### PGR-TK Decomposition

PGR-TK was evaluated using nine parameter combinations for:

- `-k`
- `-w`
- `min-span`
- `--min-branch-size`
- `-bundle-length-cutoff`

The `-r` parameter was adjusted according to average sequence length.

For each SCR:

- the optimal parameter set maximized mean principal bundle coverage

The ordered principal bundle composition of each haplotype was defined as its structural haplotype.

Scripts are available in:

```text
02_structural_haplotype/
```

---

### 5. Fst Calculation

Population differentiation between EAS and AFR populations was estimated using Hudson Fst.

For each population:

```text
H = 1 − Σ(x_i^2 / X^2)
```

where:

- `x_i` = structural haplotype count
- `X` = total haplotypes

Mean within-population dissimilarity:

```text
Hw = (H_EAS + H_AFR) / 2
```

Between-population distance:

```text
D = 1 − Σ(x_i^EAS x_i^AFR / X_EAS X_AFR)
```

Final Hudson Fst:

```text
Fst = 1 − Hw / D
```

---

### 6. Manual Refinement of Structural Haplotypes

Manual refinement was performed on the top 200 high-Fst SCRs.

Adjacent SCRs were merged when:

- >5 haplotypes carried inversions spanning intervening intervals

Merged regions included:

- `(SCR_287, SCR_799, SCR_909)`
- `(SCR_54, SCR_108)`

#### Refinement Workflow

1. Duplicated principal bundles were identified from T2T-CHM13
2. Self-dotplots were generated using `Gepard`
3. Bundle boundaries were manually refined
4. Representative bundles were mapped back using `minimap2`

Parameters:

```text
-x asm20 -p 0.2 -N 10 -g 10000 -c
```

Bundles showing truncated alignments (<80% aligned length) in >5 haplotypes were iteratively split and remapped.

Additional segments identified from bundle coverage analyses were incorporated until:

- no bundle was truncated in >5 haplotypes
- <5 haplotypes showed <95% bundle coverage

Fst values were subsequently recalculated.

---

## Great Ape Comparative Analyses

### 7. Haplotype Characterization in Great Apes

Great ape genome assemblies were obtained from:

```text
https://github.com/marbl/Primates
```

SCRs were localized in great ape genomes using `minimap2` with:

```text
-x asm20 -p 0.6 -N 10 -g 10000
```

#### REXO1L VNTR Array

For the 8q21.2 REXO1L VNTR:

- repeat motifs were defined from recurrent transposable elements
- annotations were obtained from GenomeArk
- motifs were remapped to assess tandem repeat organization

#### ACOT Loci

Manually resolved structural segments were aligned to syntenic great ape regions using the same refinement strategy described above.

---

## α-Globin Structural and Phylogenetic Analyses

### 8. α-Globin Gene Cluster Analyses

To investigate chimeric haplotypes:

- BED coordinates for chimeric regions and flanking segments were extracted
- sequences were aligned to masked T2T-CHM13 chr16 using `Winnowmap`

Parameters:

```text
-k 9 -w 2 -ax map-ont
```

Alignments were visualized using `IGV`.

#### Phylogenetic Analyses

Variants satisfying:

```text
AF ≥ 0.01
```

within ±100-kb flanking regions were extracted from pangenome-derived VCFs.

Maximum-likelihood phylogenies were reconstructed using `IQ-TREE` under:

```text
HKY+F+I
```

with:

```text
1000 bootstrap replicates
```

---

## Software

Main tools used in this workflow:

- `PAV`
- `PGR-TK`
- `halLiftover`
- `BEDTools`
- `PLINK`
- `TRF`
- `Gepard`
- `minimap2`
- `Winnowmap`
- `IQ-TREE`
- `IGV`
