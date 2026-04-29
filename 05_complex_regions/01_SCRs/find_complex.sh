#!/bin/bash
#Merge the complex intervals
grep -Ev "simple|CC1" ${complex_interval_annotation}.bed | bedtools merge -i - -d 10000 > complex.r1.v1.bed
awk '$3 -$2 >= 10000' complex.r1.v1.bed > complex.r1.v1.long.bed
awk '$3 -$2 < 10000' complex.r1.v1.bed > complex.r1.v1.short.bed
python3 ~/pyscript/Merge_bed.py -bed complex.r1.v1.long.bed -tr 0.8 -mt 2 > complex.r1merge.bed
cat complex.r1merge.bed complex.r1.v1.short.bed | cut -f 1-3 | sort -k1,1V -k2,2n | bedtools merge -i -  > temp 
python3 ~/pyscript/Merge_bed.py -bed temp -tr 0.6 -mt 2 | awk '$3-$2 >= 10000' > complex.r2merge.bed 

#Break the draft complex regions if it contain consecutive simple or CC1 intervals with length larger 40% than its length 
bedtools intersect -a <(grep -E "simple|CC1" ${complex_interval_annotation}.bed | bedtools merge -i - -d 1000) -b complex.r2merge.bed | awk '$3 - $2 > 500'> simple.bed
cp complex.r2merge.bed init.bed
numm=`cat complex.r2merge.bed | wc -l`
while true
do
sleep 1
bedtools subtract -a init.bed -b simple.bed -f 0.4 | awk '$3-$2 >= 10000' > temp.bed
num=`cat temp.bed | wc -l`
if [ $num == $numm ]; then
  break
fi
numm=${num}
cp temp.bed init.bed
done

#Refine the complex regions containing CC2 intervals with length larger 40% than its length and remove the complex regions containing only simple, CC1 and CC2 intervals
mv init.bed complex.r3merge.bed
bedtools coverage -a complex.r3merge.bed -b <(grep -E "simple|CC1" ~/project/APG/data/complex.bed | bedtools merge -i -) | bedtools coverage -a - -b <(grep -E "CC2" ~/project/APG/data/complex.bed | bedtools merge -i -) | sort -k8,8nr | awk '$8 <= 0.4' > complex.r1certain.bed 
bedtools coverage -a complex.r3merge.bed -b <(grep -E "simple|CC1" ~/project/APG/data/complex.bed | bedtools merge -i -) | bedtools coverage -a - -b <(grep -E "CC2" ~/project/APG/data/complex.bed | bedtools merge -i -) | sort -k8,8nr | awk '$8 > 0.4' | cut -f1-3 > complex.r1uncertain.bed 
bedtools intersect -a <(grep -E "simple|CC1" ~/project/APG/data/complex.bed | bedtools merge -i - -d 1000) -b complex.r1uncertain.bed | awk '$3-$2 > 5000'> simple.r2.bed
cp complex.r1uncertain.bed init.bed
numm=`cat complex.r1uncertain.bed | wc -l`
while true
do
sleep 1
bedtools subtract -a init.bed -b simple.r2.bed -f 0.1 > temp.bed
num=`cat temp.bed | wc -l`
if [ $num == $numm ]; then
 break
fi
numm=${num}
cp temp.bed init.bed
done
bedtools coverage -a init.bed -b <(grep -E "simple|CC1" ~/project/APG/data/complex.bed | bedtools merge -i -) | bedtools coverage -a - -b <(grep -E "CC2" ~/project/APG/data/complex.bed | bedtools merge -i -) | cut -f1-3 > complex.r2certain.bed 
cat complex.r1certain.bed complex.r2certain.bed | cut -f1-3 | awk '$3 -$2 > 10000'> complex.r3.10k.bed
bedtools coverage -a complex.r3.10k.bed -b <(grep -E "simple|CC1" ~/project/APG/data/complex.bed | bedtools merge -i -| bedtools subtract -a - -b <(grep -Ev "simple|CC1" ~/project/APG/data/complex.bed)) | bedtools coverage -a - -b <(grep -E "CC2" ~/project/APG/data/complex.bed | bedtools merge -i -) | sort -k7,7nr | awk '$9 / ($6 - $5) < 1' | cut -f1-3 > final.CC2_1.10k.bed
