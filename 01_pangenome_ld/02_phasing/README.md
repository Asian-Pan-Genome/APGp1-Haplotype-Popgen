\# Haplotype Phasing and Comparison



This repository contains scripts used for haplotype phasing, scaffold-based statistical phasing, and phasing accuracy benchmarking for assembly-derived and short-read (SR)-derived variant callsets.



\## Overview



We performed chromosome-wise statistical phasing using the T2T-CHM13-native 1000 Genomes Project (1kGP) reference panel and corresponding recombination maps.



Phasing analyses included:



\- trio-based phasing for benchmarking haplotypes

\- panel-based statistical phasing of SR-derived callsets

\- scaffold-guided phasing of out-of-panel variants

\- comparison between statistical and assembly-derived haplotypes

\- evaluation of switch error rate (SER) and flip error rate (FER)



Downstream analyses were restricted to autosomes.



\---



\## Repository Structure



```text

02\_phasing/

├── 01\_in\_panel

│   ├── 00.intersect\_phased\_sr\_and\_assembly\_variants.sh

│   ├── 01.sr\_error\_rate\_whatshap\_evaluation.sh

│   ├── 02.assembly\_error\_rate\_whatshap\_evaluation.sh

│   ├── 03.summarize\_switch\_flip\_error\_statistics.sh

│   ├── 04.136\_paired.calculate.ser\_fer.p\_value.R

│   └── 136.sample.child-01.list

└── 02\_out\_of\_panel

&#x20;   ├── 00.phase\_out\_of\_panel\_variants.sh

&#x20;   ├── 01.evaluate\_out\_of\_panel\_phasing\_accuracy.sh

&#x20;   ├── 02.summarize\_switch\_flip\_error\_statistics.sh

&#x20;   └── 160.sample.child-01.list

```



\---



\# Reference Resources



\## 1kGP CHM13-native reference panel



The T2T-CHM13-native phased reference panel for 3,202 samples:



\[1kGP CHM13-native phased panel](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T%2FCHM13%2Fassemblies%2Fvariants%2F1000\_Genomes\_Project%2Fchm13v2.0%2FPhased\_SHAPEIT5\_v1.1%2F\&utm\_source=chatgpt.com)



\## Recombination maps



Population-averaged and scaled genetic maps:



\[CHM13 recombination maps](https://zenodo.org/records/14891074?utm\_source=chatgpt.com)



\---



\# In-panel Phasing Evaluation



\## 1. Intersect phased SR and assembly-derived variants



Shared variants between statistically phased SR callsets and assembly-derived haplotypes were identified using `bcftools isec`.



```bash

01\_in\_panel/00.intersect\_phased\_sr\_and\_assembly\_variants.sh

```



\---



\## 2. Evaluate SR-based phasing accuracy



Trio-phased haplotypes were used as benchmark truth sets for evaluating SR-derived statistical phasing.



```bash

01\_in\_panel/01.sr\_error\_rate\_whatshap\_evaluation.sh

```



\---



\## 3. Evaluate assembly-derived phasing accuracy



Assembly-derived haplotypes were compared against trio-phased benchmark haplotypes.



```bash

01\_in\_panel/02.assembly\_error\_rate\_whatshap\_evaluation.sh

```



\---



\## 4. Summarize SER and FER statistics



Switch error rate (SER) and flip error rate (FER) were summarized across autosomes.



```bash

01\_in\_panel/03.summarize\_switch\_flip\_error\_statistics.sh

```



Paired statistical testing:



```bash

01\_in\_panel/04.136\_paired.calculate.ser\_fer.p\_value.R

```



\---



\# Out-of-panel Phasing Evaluation



\## 1. Scaffold-based statistical phasing



Variants absent from the 1kGP reference panel were phased using SHAPEIT5 scaffold-based phasing.



```bash

02\_out\_of\_panel/00.phase\_out\_of\_panel\_variants.sh

```



Assembly-derived haplotypes without missingness at overlapping in-panel sites were used as scaffolds.



\---



\## 2. Evaluate out-of-panel phasing accuracy



Scaffold-phased variants were benchmarked against assembly-derived haplotypes using WhatsHap.



```bash

02\_out\_of\_panel/01.evaluate\_out\_of\_panel\_phasing\_accuracy.sh

```



\---



\## 3. Summarize out-of-panel SER and FER



```bash

02\_out\_of\_panel/02.summarize\_switch\_flip\_error\_statistics.sh

```



\---



\# Main Software



\- \[SHAPEIT5](https://odelaneau.github.io/shapeit5/?utm\_source=chatgpt.com)

\- \[WhatsHap](https://whatshap.readthedocs.io/en/latest/?utm\_source=chatgpt.com)

\- \[bcftools](https://samtools.github.io/bcftools/?utm\_source=chatgpt.com)



\---



\# Metrics



\## Switch Error Rate (SER)



SER measures the proportion of heterozygous variants incorrectly phased relative to the preceding heterozygous site.



\## Flip Error Rate (FER)



FER measures the proportion of consecutive heterozygous variants exhibiting back-to-back switch errors.



\---



\# Citation



If you use this repository, please cite the corresponding study and the original software packages referenced above.

