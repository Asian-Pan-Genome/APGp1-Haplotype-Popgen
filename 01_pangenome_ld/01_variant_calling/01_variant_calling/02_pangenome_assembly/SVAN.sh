#!/bin/bash
svan={path to SVAN}
export PYTHONPATH=SVAN/GAPI:$PYTHONPATH
ref={CHM13v2}.fasta

ins=CHM13-APGp1-HPRCp1-HGSVCp3_MC.autosome.ex_tc.final.dip.ins.SVs.vcf
del=CHM13-APGp1-HPRCp1-HGSVCp3_MC.autosome.ex_tc.final.dip.del.SVs.vcf
outDir=svan

set -euo pipefail

date

python $svan/scripts/ins2fasta.py ${ins} ${outDir}

trf ${outDir}/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs 1 > ${outDir}/ins_trf.out

python $svan/scripts/del2fasta.py ${del} ${outDir}

trf ${outDir}/deletions_seq.fa  2 7 7 80 10 10 500 -h -d -ngs 1 > ${outDir}/del_trf.out

if [ $? -ne 0 ]; then
    echo "Error: step1 failed."
    exit 1
fi

python $svan/SVAN-INS.py --replace-info -o ${outDir} ${ins} ./ins_trf.out \
 $svan/chm13/VNTR_chm13.bed $svan/chm13/EXONS_chm13.bed $svan/chm13/REPEATS_chm13.bed $svan/chm13/CONSENSUS.fa ${ref} ins.SAMPLEID

python $svan/SVAN-DEL.py --replace-info -o ${outDir} ${del} ./del_trf.out \
 $svan/chm13/VNTR_chm13.bed $svan/chm13/EXONS_chm13.bed $svan/chm13/REPEATS_chm13.bed $svan/chm13/CONSENSUS.fa ${ref} del.SAMPLEID

status=$?

echo

if [ $status -eq 0 ]; then
    echo "Steps all completed successfully."
    echo
    echo "finish!"
    echo
    echo $(date)
else
    echo "Step exited with error."
    echo
    echo "exit with error!!"
    echo
    echo $(date)
fi
