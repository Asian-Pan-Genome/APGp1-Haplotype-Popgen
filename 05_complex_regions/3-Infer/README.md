# SV acienstro sequence analysis

This repository contains a pipeline for evauate ancestro state of SVs in human.

## Requirements
- Apes genome
- muscle3

## Pipeline
First, the flanking 5Kb cordinates and SV refrenece and alternative sequences with its 5Kb flanking sequence are prepared:
```
bash

for i in `seq 1 22` X Y; do
chrom="chr${i}"
bcftools view -e "INFO/MAF <= 0.05" ${chrom}.vcf.gz | bcftools query -f '%CHROM\t%POS\t%END\t%TYPE\t%REF\t%ALT\n' | awk 'OFS="\t"{print $1, $2-1,$3,$4,$5,$6}' 
done | awk 'OFS="\t" {print $1,$2,$3,$4"_"NR,$5,$6}' > SV.maf5.bed
python3 alelle_getter.py -b SV.maf5.bed -fa {CHM13}.fasta
```
Then, the flanking 5Kb is lifted to the Gorilla, Bonobo and chimpanzee and only the SV with which left and right flankings can be lifted to uniq position in a chromosome are kepted:
```
bash

## Lift over the SV flanking regions to the ape genomes and filter the lifted regions. Then, extract the sequences of the lifted regions for further analysis.
awk 'OFS="\t" {print $1,$2-5000,$2,$4"|left"; print $1,$3, $3+5000,$4"|right"}' SV.maf5.bed > SV.maf5.flank.info

##The idcheck.py script is used to check if the left and right flanking regions of the same SV are located on the same chromosome. 
##Chimpanzee
liftOver SV.maf5.flank.info {hs1_vs_pantro}.chain.gz SV.maf5.flank.pantro.info SV.maf5.flank.pantro.missing -minMatch=0.80 -multiple
cut -f1,4 SV.maf5.flank.pantro.info | sed "s/|/ /g" | python3 idcheck.py > pantro.finelift.ids
sed "s/|/\t/g" SV.maf5.flank.pantro.info | awk '{print $4"|"$1"\t"$0}' | sort -k1,1 | join -1 1 -2 1 - <(awk '{print $1"|"$2}' pantro.finelift.ids | sort -k1,1) | sed -f <(awk '{print "s/"$1"/"$3"/g"}' pantro.txt) | awk 'OFS="\t" {print $2,$3,$4,$5,$6}' | sort -k1,1V -k4,4V -k5,5 > SV.maf5.flank.pantro.finelift.info
awk 'OFS="\t" {if (NR % 2 == 1) {start = $2; end = $3} else {if(start < $2) {print $1,start,$3,$4,"0","+"} else {print $1,$2,end,$4,"0","-"}}}' SV.maf5.flank.pantro.finelift.info | awk '$3-$2 < 200000'> SV.maf5.flank.pantro.finelift.bed
python3 pyscript/cutbyBed.py -bed SV.maf5.flank.pantro.finelift.bed -fasta {mPanTro3.analysis-dip.20231122}.fasta -infoonly -dic > SV.maf5.flank.pantro.finelift.fa 

##Bonobo
liftOver SV.maf5.flank.info {hs1_vs_panpan}.chain.gz SV.maf5.flank.panpan.info SV.maf5.flank.panpan.missing -minMatch=0.80 -multiple
cut -f1,4 SV.maf5.flank.panpan.info | sed "s/|/ /g" | python3 idcheck.py > panpan.finelift.ids
sed "s/|/\t/g" SV.maf5.flank.panpan.info | awk '{print $4"|"$1"\t"$0}' | sort -k1,1 | join -1 1 -2 1 - <(awk '{print $1"|"$2}' panpan.finelift.ids | sort -k1,1) | sed -f <(awk '{print "s/"$1"/"$3"/g"}' panpan.txt) | awk 'OFS="\t" {print $2,$3,$4,$5,$6}' | sort -k1,1V -k4,4V -k5,5 > SV.maf5.flank.panpan.finelift.info
awk 'OFS="\t" {if (NR % 2 == 1) {start = $2; end = $3} else {if(start < $2) {print $1,start,$3,$4,"0","+"} else {print $1,$2,end,$4,"0","-"}}}' SV.maf5.flank.panpan.finelift.info | awk '$3-$2 < 200000'> SV.maf5.flank.panpan.finelift.bed
python3 pyscript/cutbyBed.py -bed SV.maf5.flank.panpan.finelift.bed -fasta {mPanPan1.analysis-dip.20231122}.fasta -infoonly -dic > SV.maf5.flank.panpan.finelift.fa 

##Gorilla
liftOver SV.maf5.flank.info {hs1_vs_gorgor}.chain.gz SV.maf5.flank.gorgor.info SV.maf5.flank.gorgor.missing -minMatch=0.80 -multiple
cut -f1,4 SV.maf5.flank.gorgor.info |  sed "s/|/ /g" | python3 idcheck.py> gorgor.finelift.ids
sed "s/|/\t/g" SV.maf5.flank.gorgor.info | awk '{print $4"|"$1"\t"$0}' | sort -k1,1 | join -1 1 -2 1 - <(awk '{print $1"|"$2}' gorgor.finelift.ids | sort -k1,1) | sed -f <(awk '{print "s/"$1"/"$3"/g"}' gorgor.txt) | awk 'OFS="\t" {print $2,$3,$4,$5,$6}' | sort -k1,1V -k4,4V -k5,5 > SV.maf5.flank.gorgor.finelift.info
awk 'OFS="\t" {if (NR % 2 == 1) {start = $2; end = $3} else {if(start < $2) {print $1,start,$3,$4,"0","+"} else {print $1,$2,end,$4,"0","-"}}}' SV.maf5.flank.gorgor.finelift.info | awk '$3-$2 < 200000'> SV.maf5.flank.gorgor.finelift.bed
python3 pyscript/cutbyBed.py -bed SV.maf5.flank.gorgor.finelift.bed -fasta {mGorGor1.analysis-dip.20231122}.fasta -infoonly -dic > SV.maf5.flank.gorgor.finelift.fa 
``` 
Finally, the reference and alternatinve SV sequences with the 5Kb flanking sequence and the corresponding 3 apes sequence is aligned via muscle.
```
bash
cat *.ids | cut -f1 | sort | uniq | while read a; do

   v_pos=`awk -v idx=${a} '{if ($4 == idx) {print $2"-"$3}}' SV.maf5.bed`
   f_pos=`sed "s/|/\t/g" SV.maf5.flank.info | awk -v idx=${a} '{if ($4 == idx) {if($5 == "left") {st = $2} else {ed = $3}}} END{print st"-"ed}'`

   for sp in gorgor panpan pantro; do
      seqtk subseq SV.maf5.flank.${sp}.finelift.fa <(printf "${a}") > temp/${sp}.fa
      if [ ! -s "./temp/${sp}.fa" ]; then
      rm temp/${sp}.fa
      continue
   fi
   sed "s/${a}/${sp}/g" temp/${sp}.fa
   rm temp/${sp}.fa
done | cut -f1 -d " " > temp/apes.fa

cat <(sed "s/${a}/REF/g" temp/ref.fa) <(sed "s/${a}/ALT/g" temp/alt.fa) > temp/human.fa
muscle3 -seqtype dna -quiet -in temp/human.fa -out temp/human.aln

if [ ! -s "temp/apes.fa" ]; then
cp temp/human.aln temp/all.aln
else

muscle3 -seqtype dna -quiet -in temp/apes.fa -out temp/apes.aln
muscle3 -profile -in1 temp/human.aln -in2 temp/apes.aln -out temp/all.aln -seqtype dna

fi

python3 ../stat_via_muscle.py -v ${v_pos} -f ${f_pos} -idx "${a}" -n temp/${tag}.aln

rm temp/all.aln
rm temp/human.aln
rm temp/apes.aln
done > SVid.info
```