# Complex region determination

This repository contains pipeline for acquiring complex regions and lift the regions to other genome.

First, the complex interval annotation ```${complex_interval_annotation}.bed``` is derived from the pipeline (https://github.com/Hanjunmin/complex_region/tree/main/01_complexregion_decomposition) and ```find_complex.sh``` is applied to get the complex regions. 

The resulting complex regions file ```final.CC2_1.10k.bed``` is extended 5 Kbp both sides with ```bedtools slop``` as
```bash 
bedtools slop -i final.CC2_1.10k.bed -g ${CHM13v2}.fasta.fai -b 5000 | awk 'OFS="\t" {print $0,"COMPLEX_"NR}' > chm13.complex_et5k.bed
```
For each extended regions in ```chm13.complex_et5k.bed```, the region is lifted to every sample via halLiftover in the repository (https://github.com/ComparativeGenomicsToolkit/hal) and merged via ```bedtools merge``` as
```bash
halLiftover CHM13v2_APGp1-HPRCp1-HGSVCp3_MC_${chrom}.full.hal CHM13_Hap_${chrom} CHM13.region.bed ${sample}.lift_regions.bed
sort -k1,1 -k2,2n ${sample}.lift_regions.bed | bedtools merge -i - -d 10000000 > ${sample}.region.bed
```
