#!/bin/bash
complex=$1
kmer=$2
win=$3
span=$4
bt=$5
rv=`awk 'BEGIN{OFMT = "%i"}{sum = sum + $2; t =t + 1}END{ava = sum / t; s = sqrt(ava/50000); if(s <2){print "2"} else if(s > 12){print "12"} else{print s}}' ${complex}.fa.fai`
name=`echo $complex | cut -f2 -d "."`
pgr-pbundle-decomp -k ${kmer} -w ${win} -r ${rv} --min-span ${span} --min-branch-size ${bt} --bundle-length-cutoff 1000 ${complex}.fa "pgrtkfilter_${name}.${kmer}.${win}.${span}.${bt}."
