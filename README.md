# APGp1-Haplotype-Popgen
This repository hosts a comprehensive data and workflow resource for **assembly-based population genetic analysis** in the **APGp1** genome assembly.

---

## Overview

Main contents:

* [01](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/01_pangenome_ld) - variants and phasing comparison
* [02](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/02_recombination_map) - recombination landscape and PRDM9 genotyping
* [03](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/03_arg_selection) - ARG reconstruction and ARG-based selection scan
* [04](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/04_adaptive_sv) - positively selected SVs and functional relevance
* [05](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/05_complex_regions) - LD, recombination feature, and selection in genomic complex regions
* [06](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/06_genotyping_panel) - pangenome-based genotyping for short-read data

See each subfolder's `README.md` for details.

---

## Data Access

All **source data** and available download/request links are shown in the corresponding subfolders.
Navigation:

- Genomes & Variant callsets:
  + [GRCh38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/)
  + [CHM13v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)
  + [CN1v1.0](https://genome.zju.edu.cn/Downloads)
  + [HG002](https://github.com/marbl/HG002)
  + [APGp1 assemblies](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA030428) (NOTE: Access to the APGp1 assemblies is restricted to protect participant confidentiality. For details, see the [APGp1 repository](https://github.com/Asian-Pan-Genome/).)
  + [HPRCy1 assemblies](https://data.humanpangenome.org/assemblies)
  + [HGSVC3 assemblies](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/)
  + [APGp1 + HPRCy1 + HGSVC3 pangenome graph](https://genome.zju.edu.cn/APG/Resources#graphs)
  + [Great apes T2T assemblies](https://github.com/marbl/Primates)
  + [CHM13 based 1kGP 3202 reference panel](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/)
  + [GRCh38 based 1kGP 3202 reference panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
  + [GnomAD SVs](https://registry.opendata.aws/broad-gnomad) (v4.1.0)

- Recombination map: 
  + [1KGP T2T-CHM13v2.0 recombination maps](https://zenodo.org/records/14891074)

- Publicly available functional summary statistics:
  + [QTL sites of GTEx v10](https://www.gtexportal.org/home/downloads/adult-gtex/qtl)
  + [GWAS catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads)
  + [Regulatory elements annotation](https://downloads.wenglab.org/GRCh38-cCREs.Lifted-hs1.bed)


- NGS BAM files for genotyping: [sample list](https://github.com/Asian-Pan-Genome/APGp1-Haplotype-Popgen/tree/main/06_genotyping_panel)


## Contact

For questions or issues, please contact the corresponding author or open an issue on the GitHub repository.
