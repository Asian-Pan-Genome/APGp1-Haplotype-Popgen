# UK Biobank Structural Variant Association Analysis

This repository contains the scripts and workflows used for structural variant (SV) genotyping and association analyses in the UK Biobank (UKB), including proteomic, metabolomic, and phenome-wide association studies (PheWAS).

---

# Overview

A total of 511 candidate SVs identified on the T2T-CHM13v2.0 reference were lifted over to GRCh38 and genotyped across 490,086 UK Biobank individuals using DRAGEN SV callsets. After merging and collapsing overlapping calls, a final set of 970 SV loci was obtained for downstream association analyses.

Association testing was performed for:

- Plasma proteomic biomarkers
- Metabolomic traits
- Disease phenotypes (PheWAS)

All analyses used additive SV genotype dosage encoding (0, 1, 2).

---

# Workflow

## 1. SV LiftOver and Genotyping

### LiftOver to GRCh38

Candidate SVs were converted from T2T-CHM13v2.0 coordinates to GRCh38 using:

- CrossMap v0.7.3

### SV Matching Against UKB DRAGEN Callsets

SV matching between candidate SVs and UKB individual SV VCFs was performed using:

- Truvari v5.3.0
  - `truvari bench`
  - `truvari refine`
  - `truvari collapse`

To improve matching robustness across different SV representations:

- All SV genotypes were fixed to `1/1`
- Matching was performed at the variant level
- Breakpoint uncertainty and allele decomposition were tolerated

### Example Commands

```bash
truvari bench \
    -b sv.hg38.vcf.gz \
    -c $sample.fix_gt.vcf.gz \
    -o $sample \
    -r 2000 \
    --passonly \
    --pick ac \
    -C 5000 \
    -f hg38.fa

truvari refine \
    --use-original-vcfs \
    -a mafft \
    -t $threads
```

```bash
truvari collapse \
    -k common \
    -r 2000 \
    -p 0.7 \
    -P 0.7 \
    -D \
    --passonly \
    -f hg38.fa \
    -i $sample.tp.vcf.gz \
    -o $sample.tp.collapsed.vcf.gz \
    -c $sample.tp.removed.vcf.gz
```

After processing all individuals:

- Sample IDs were renamed to UKB EIDs
- Individual VCFs were merged and collapsed
- Final SV callset contained 970 SVs

---

# 2. Sample Quality Control

Individuals were retained using stringent genetic QC criteria.

## Inclusion Criteria

- Concordant self-reported and genetically inferred sex
- European/Caucasian ancestry
- Included in UKB principal-component calculation set
- Complete covariate information:
  - age
  - sex
  - PC1–PC10

## Exclusion Criteria

- Sex chromosome aneuploidy
- ≥10 third-degree relatives
- Heterozygosity outliers
- Missingness outliers

After QC filtering:

- Final analysis cohort: 328,896 individuals

For each association test, individuals with missing:

- SV genotype
- phenotype value
- covariates

were removed.

---

# 3. Proteomic and Metabolomic Association Analysis

## Data Processing

Quantitative biomarker measurements were merged with:

- SV genotype matrix
- covariate table

using participant IDs.

All biomarkers were transformed using:

- Rank-based inverse normal transformation (RINT)

## Association Model

Ordinary least-squares regression was fitted under an additive genetic model:

```text
RINT(biomarker) ~ SV genotype + age + sex + PC1 + ... + PC10
```

where the biomarker represents either:

- plasma protein abundance
- metabolomic trait

## Filtering

Tests were excluded if:

- fewer than 100 complete samples
- non-varying SV genotype

## Output

For each SV-biomarker pair, the following were recorded:

- regression coefficient
- standard error
- P-value

Bonferroni correction was applied separately for:

- proteomic associations
- metabolomic associations

---

# 4. Phenome-Wide Association Analysis (PheWAS)

## ICD-10 to PheCode Conversion

ICD-10 diagnosis records were converted to PheCodes using:

- R package: PheWAS

PheCode processing included:

- exclusion rules
- sex-specific filtering
- minimum code count = 2

## Association Model

Logistic regression was performed under an additive genetic model:

```text
PheCode status ~ SV genotype + age + sex + PC1 + ... + PC10
```

## Filtering

Models were excluded if they contained:

- insufficient case/control counts
- non-varying phenotypes
- non-varying SV genotypes

## Output

For each SV-PheCode pair, the following statistics were recorded:

- odds ratio
- regression coefficient
- standard error
- P-value

Bonferroni correction was applied across PheWAS tests.

---

# Software

| Software | Version |
|---|---|
| CrossMap | v0.7.3 |
| Truvari | v5.3.0 |
| MAFFT | integrated in Truvari refine |
| R PheWAS package | latest |

---

# Input Data

## Required Inputs

- Candidate SV VCF (T2T-CHM13v2.0)
- LiftOver chain files
- UKB DRAGEN SV VCFs
- UKB covariates
- Proteomic measurements
- Metabolomic traits
- ICD-10 records
- Reference genome (GRCh38 / hg38)

---

# Output Files

## SV Genotyping

- Individual matched SV VCFs
- Collapsed TP SV VCFs
- Merged cohort-level SV matrix

## Association Results

- Proteomic association summary statistics
- Metabolomic association summary statistics
- PheWAS summary statistics
- Multiple-testing corrected results

---

# Citation

If you use this workflow or codebase, please cite the corresponding study and the following software:

- CrossMap
- Truvari
- PheWAS R package

---

# References

- Zhao et al. 2014. CrossMap
- English et al. 2022. Truvari
- Carss et al. 2025. UKB DRAGEN SV callsets
- PheWAS R package: https://github.com/PheWAS/PheWAS
