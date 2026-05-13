python3 merge.py -p pgrtkhap.info -m m.info > hap.info
python3 sv_ana.py -split -org EAS -info hap.info -samorg ~/project/APG/sample.info | cut -f1-5,7 > complex.EAS.info
python3 sv_ana.py -split -org AFR -info hap.info -samorg ~/project/APG/sample.info | cut -f1-5,7 > complex.AFR.info
awk 'OFS="\t" {print "NA","EAS",$4,$2,$3,$4,$5,"+","NA",$6}' complex.EAS.info > temp.EAS.info
awk 'OFS="\t" {print "NA","AFR",$4,$2,$3,$4,$5,"+","NA",$6}' complex.AFR.info > temp.AFR.info
python3 par_cal.py -info1 temp.EAS.info -info2 temp.AFR.info -Fst
