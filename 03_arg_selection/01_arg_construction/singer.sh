#!/bin/bash -l



#source pyenv/bin/activate



python 03_arg_selection/01_arg_construction/run_singer_convert.py \
 -mut_map ratemaps/fake_mutmap_chr6.txt \
 -recomb_map ratemaps/pyrho.chm13/ratemap.avgmask.chr6.txt \
 -start 29e6 \
 -end 34e6 \
 -vcf_header 100_lr_data/2.9e7_3.4e7_CHM13-APGp1_MC.chr6.filtered_missing10.polarized.no_missing \
 -arg_save singer_arg/100_lr/2.9e7_3.4e7_CHM13-APGp1_MC.chr6.filtered_missing10.polarized.no_missing \
 -n 100 \
 -thin 100 \