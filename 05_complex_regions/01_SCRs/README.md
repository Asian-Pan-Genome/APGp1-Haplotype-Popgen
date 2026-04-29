# Complex region identification

This repository contains a pipeline for identifying complex regions and lifting them to other genomes.

First, the complex interval annotation ```${complex_interval_annotation}.bed``` is obtained from the pipeline at https://github.com/Hanjunmin/complex_region/tree/main/01_complexregion_decomposition. The ```find_complex.sh``` iscript is then applied to extract the complex regions.

The resulting complex regions file, ```final.CC2_1.10k.bed```, is extended by 5 kb on both sides using ```bedtools slop```:
```bash 
bedtools slop -i final.CC2_1.10k.bed -g ${CHM13v2}.fasta.fai -b 5000 | awk 'OFS="\t" {print $0,"COMPLEX_"NR}' > chm13.complex_et5k.bed
```
For each extended regions in ```chm13.complex_et5k.bed```, the region is lifted to each sample using ```halLiftover``` from the [hal](https://github.com/ComparativeGenomicsToolkit/hal) repository. The lifted region is then merged with  ```bedtools merge```:
```bash
halLiftover CHM13v2_APGp1-HPRCp1-HGSVCp3_MC_${chrom}.full.hal CHM13_Hap_${chrom} CHM13.region.bed ${sample}.lift_regions.bed
sort -k1,1 -k2,2n ${sample}.lift_regions.bed | bedtools merge -i - -d 10000000 > ${sample}.region.bed
```
