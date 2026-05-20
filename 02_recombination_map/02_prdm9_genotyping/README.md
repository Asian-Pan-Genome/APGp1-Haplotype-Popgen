# Requirements
- RepeatMasker
- fxTools
- Biopython


# Extracting ZF sequences
For reference T2T-CHM13v2, we first compiled a ZFs library from [Alleva et al.](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2021.675286/full), which was ued by RepeatMasker to locate the range of ZFs:
```shell
RepeatMasker -pa 6 -nolow -lib ZF.fa CHM13v2.PRDM9.exon11.fa
```

For other assemblies, ZF sequences were extracted from the pangenome HAL file using the command:
```shell
halLiftover CHM13v2_APGp1-HPRCp1-HGSVCp3_MC_chr5.full.hal CHM13_Hap_chr5 <(echo -e "CHM13_Hap_chr5\t23634386\t23635478") ${sample_name} ${sample_name}.bed
start=`head -1 ${sample_name}.bed | cut -f 2`
end=`tail -1 ${sample_name}.bed | cut -f 3`
chrom=`rg '>' $fa_file | sed 's/>//g'`
fxTools getseq $fa_file -r <(echo -e "$chrom\t$start\t$end") > ${sample_name}.fa
```

# ZF and PRDM9 alleles genotyping
We assigned known ZFs from [Alleva et al.](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2021.675286/full) to each 84-bp unit to annotate the ZF and PRDM9 alleles.
```shell
python assign_PRDM9_alleles.clean.py \
  input.samples.ZFs.fa \
  PRDM9.alleles.lib.tsv \
  --allele-out samples.PRDM9_allele.list \
  --new-zf-out ZF.new.fa
```
