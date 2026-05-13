## PGR-TK based haplotype structrue identification
### Prerequisites
- [trf](https://github.com/Benson-Genomics-Lab/TRF)
- [PGR-TK](https://github.com/GeneDx/pgr-tk)

### Workflow
For each collection of sequence in complex region among samples, the simple repeats with length over 600bp and identity over 80% are first masked:
```
bash
trf ${complex}.fa 2 7 7 80 10 50 2000 -d -h
awk 'BEGIN{chr=""; repeat_counter=1; OFS="\t"} {if($1 == "Sequence:") {chr=$3; repeat_counter=1} else if(NF==15){print chr,$1-1,$2,chr"_repeat_"repeat_counter++,$3,$4,$6,$14,$15}}' ${complex}.fa.2.7.7.80.10.50.2000.dat | bedtools sort -i - | awk '$3 - $2 > 600 && $7 >= 80' > all.fa.trf.bed
python3 01.parse_trf_ouput.py all.fa.trf.bed 0.3 0.1 all.fa.trf.refined.bed all.fa.trf.refined.motifs.fa 
sed 's/:/\t/g' all.fa.trf.refined.bed | cut -f1,3,4 > mask.bed
python3 getInterval.py -mask -all --cutfile mask.bed -fasta ${complex}.fa > ${complex}.mask.fa 
samtools faidx ${complex}.mask.fa 
```

Next, [PGR-TK](https://github.com/GeneDx/pgr-tk) is applied to the resulting masked FASTA using various parameter combinations:
```bash
pgrtk.sh ${complex}.mask 56 48 12 32 
pgrtk.sh ${complex}.mask 31 20 12 32
pgrtk.sh ${complex}.mask 17 10 8 32
pgrtk.sh ${complex}.mask 56 48 12 16 
pgrtk.sh ${complex}.mask 31 20 12 16
pgrtk.sh ${complex}.mask 17 10 8 16
pgrtk.sh ${complex}.mask 56 48 12 8
pgrtk.sh ${complex}.mask 31 20 12 8
pgrtk.sh ${complex}.mask 17 10 8 8
```
Finally, ```best_hap_getter.sh``` is used to obtain the final haplotype structure of this complex region.