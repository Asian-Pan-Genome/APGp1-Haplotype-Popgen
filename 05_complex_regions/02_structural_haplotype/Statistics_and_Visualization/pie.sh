#!/bin/bash
python3  sv_ana.py -split -org EAS -info hap.info -samorg ~/project/APG/sample.info | cut -f4,5 | sort -k1,1 > complex.EAS.info
python3  sv_ana.py -split -org AFR -info hap.info -samorg ~/project/APG/sample.info | cut -f4,5 | sort -k1,1 > complex.AFR.info
python3  sv_ana.py -split -org SAS -info hap.info -samorg ~/project/APG/sample.info | cut -f4,5 | sort -k1,1 > complex.SAS.info
python3  sv_ana.py -split -org EUR -info hap.info -samorg ~/project/APG/sample.info | cut -f4,5 | sort -k1,1 > complex.EUR.info
python3  sv_ana.py -split -org AMR -info hap.info -samorg ~/project/APG/sample.info | cut -f4,5 | sort -k1,1 > complex.AMR.info
paste complex.EAS.info complex.AFR.info complex.SAS.info complex.EUR.info complex.AMR.info | cut -f1,2,4,6,8,10 | join - <(sort -k1,1 hapTohap.info) | sort -k7,7V | awk 'BEGIN{print "EAS\tAFR\tSAS\tEUR\tAMR\tHAP"; OFS="\t"}{print $2,$3,$4,$5,$6,$7}' > pop.stat
Rscript pie.R
Rscript struc.R
