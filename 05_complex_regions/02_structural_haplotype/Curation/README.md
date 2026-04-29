## Manually curation pipeline

This repository provides key steps for manually curating structural haplotypes.

### Prerequisites
- [gepard](https://github.com/univieCUBE/gepard)
- [minimap2](https://github.com/lh3/minimap2)

### Workflow

First, the PGR‑TK decomposition visualization and self‑dot plot of each complex region in CHM13 are generated:
```
bash
pgr-pbundle-bed2svg ${CHM13_complex}.pgrtk.bed chm13_pgrtk
java -cp ${gepard}/dist/Gepard-2.1.jar org.gepard.client.cmdline.CommandLine -seq ${CHM13_complex}.fa ${CHM13_complex}.fa -matrix ${gepard}/resources/matrices/edna.mat -outfile ${CHM13_complex}.bmp -format bmp -word 25
```

A self‑alignment is also generated to help curate the PGR‑TK decomposition:
```
bash
minimap2 -xasm20 -X -p 0.2 -N 10 -g 10000 -c ${CHM13_complex}.fa ${CHM13_complex}.fa > seg.paf
```

For each complex region, we then adjust the bundle boundaries to match the duplicated segments revealed by the self‑dot plot and document them in a BED file ```chm13.seg.bed```. The decompositions based on the curated bundle of other samples are generated as:
```
bash
minimap2 -xasm20 -p 0.2 -N 10 -g 10000 -c chm13.seg.fa ${complex}.fa > decom.sim.paf
awk 'OFS="\t" {if(($9-$8)/$7 > 0.6) print $1,$3,$4,$6,$10,$10/$11,$5}' decom.sim.paf | sort -k1,1 -k2,2n | bedmap --max-element --fraction-either 0.8 - | uniq | awk 'OFS="\t" {print $1,$2,$3,$4,$6,$7}' | bedmap --max-element --fraction-either 0.8 - | uniq > decom.notr.bed
sed -e 's/dup1_1/dup1/g' -e 's/dup1_2/dup1/g' decom.notr.bed | sort -k4,4 | join -1 4 -2 1 - <(cut -f4 decom.notr.bed |sed -e 's/dup1_1/dup1/g' -e 's/dup1_2/dup1/g' | sort | uniq -c | awk 'OFS=":" {print $2"\t"NR,"0\t0",$1,"U"}') | awk 'OFS="\t" {if($6 == "+") {print $2,$3,$4,$7":0:"$8}else {print $2,$3,$4,$7":1:"$8}}' | sort -k1,1 -k2,2n > decom.pgrtklike.bed 
bedtools coverage -a <(awk 'OFS="\t" {print $1,0,$2}' /share/home/zhanglab/user/nielei/project/APG/haplotype/PGRTK/${chrom}/${complex}/${complex}.rename.fa.fai) -b decom.pgrtklike.bed | cut -f1,7 > sample.cov
```

We check the ```sample.cov``` file for the bundle coverage of each sample and ```decom.sim.paf``` for truncated bundles. We repeatedly generate the CHM13 decomposition file ```chm13.seg.bed``` and perform the above process until:
- (i) no more than five haplotypes showed a truncated alignment for any bundle, and 
- (ii) no more than five haplotypes exhibited <95% coverage
