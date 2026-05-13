# prepare sample list
bcftools view CHM13-APGp1-HPRCp1-HGSVCp3_MC.chr19.ex_tc.final.SNPs.biallelic.addtag.filter.vcf.gz \
    | grep '^#CHROM' \
    | cut -f10- \
    | tr '\t' '\n' \
    > CHA_sample.txt

# index the VCF files
for i in {1..22}
do
bcftools index -f CHM13-APGp1-HPRCp1-HGSVCp3_MC.chr${i}.ex_tc.final.SNPs.biallelic.addtag.filter.vcf.gz
done





# vcf2smc + smcpp

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

CORES=24

# population name
pop="CHA"

# get comma-separated sample list for CHA
samples=$(awk '{print $1}' CHA_sample.txt | paste -sd, -)
echo "Samples for $pop: $samples"


# loop over chromosomes
for chr in {1..22}; do
    echo "  Running vcf2smc for chr $chr"
    VCF=CHM13-APGp1-HPRCp1-HGSVCp3_MC.chr${chr}.ex_tc.final.SNPs.biallelic.addtag.filter.vcf.gz
    singularity exec -B "$PWD":/mnt /share/home/zhanglab/user/weixinzhu/Lin/smcpp/smcpp_latest.sif \
        smc++ vcf2smc \
        --ignore-missing \
        --missing-cutoff 100000 \
        --cores $CORES \
        $VCF \
        ${pop}_vcf2smc_chr${chr} \
        chr${chr} \
        "$pop:$samples"
done


echo "Job of smcpp estimate started on $(hostname) at $(date)"

smc_files=$(for chr in {1..22}; do echo ${pop}_vcf2smc_chr${chr}; done | paste -sd ' ' -)

echo " Running smc++ estimate for $pop"
singularity exec -B "$PWD":/mnt /share/home/zhanglab/user/weixinzhu/Lin/smcpp/smcpp_latest.sif \
    smc++ estimate \
    --polarization-error 0.5 \
    --cores $CORES \
    --nonseg-cutoff 100000 \
    -o ./${pop}_pchip \
    --spline pchip \
    --timepoints 100 500000 \
    1.25e-8 \
    $smc_files

echo "Job finished at $(date)"


