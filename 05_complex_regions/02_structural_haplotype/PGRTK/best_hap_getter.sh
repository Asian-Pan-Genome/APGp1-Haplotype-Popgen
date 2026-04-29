#!/bin/bash

p1=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.17.10.8.32.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p2=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.31.20.12.32.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p3=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.56.48.12.32.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p4=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.17.10.8.16.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p5=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.31.20.12.16.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p6=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.56.48.12.16.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p7=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.17.10.8.8.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p8=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.31.20.12.8.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`
p9=`bedtools coverage -a <(awk 'OFS="\t" {print $1,"0",$2}' ${complex}.rename.fa.fai) -b pgrtkfilter_maskkeep.56.48.12.8.bed |  awk '{sum = sum + $NF; l = l + 1}END{print sum/l}'`

pp8=$(echo "$p8 > $p9" | bc)
if [ $pp8 -ne 1 ]; then p8=${p9}; fi
pp7=$(echo "$p7 > $p8" | bc)
if [ $pp7 -ne 1 ]; then p7=${p8}; fi
pp6=$(echo "$p6 > $p7" | bc)
if [ $pp6 -ne 1 ]; then p6=${p7}; fi
pp5=$(echo "$p5 > $p6" | bc)
if [ $pp5 -ne 1 ]; then p5=${p6}; fi
pp4=$(echo "$p4 > $p5" | bc)
if [ $pp4 -ne 1 ]; then p4=${p5}; fi
pp3=$(echo "$p3 > $p4" | bc)
if [ $pp3 -ne 1 ]; then p3=${p4}; fi
pp2=$(echo "$p2 > $p3" | bc)
if [ $pp2 -ne 1 ]; then p2=${p3}; fi
pp1=$(echo "$p1 > $p2" | bc)

if [ $pp1 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.17.10.8.32.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp2 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.31.20.12.32.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp3 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.56.48.12.32.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp4 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.17.10.8.16.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp5 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.31.20.12.16.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp6 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.56.48.12.16.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp7 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.17.10.8.8.bed -chr ${chrom} -st ${start} -ed ${end}
elif [ $pp8 -eq 1 ]; then
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.31.20.12.8.bed-chr ${chrom} -st ${start} -ed ${end}
else
    python3 bedtohap.py -bed pgrtkfilter_maskkeep.56.48.12.8.bed -chr ${chrom} -st ${st} -ed ${ed}
fi
