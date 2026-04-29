set -euo pipefail

POP=CHA
BPEN=35
W=60

echo "Job started on $(hostname) at $(date)"
echo "Population: $POP, bpen: $BPEN, w: $W"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "CPU threads: $SLURM_CPUS_PER_TASK"

# Load conda
source /share/home/zhanglab/user/weixinzhu/tools/miniconda3/etc/profile.d/conda.sh
conda activate pyrho

# Loop over chromosomes
for chr in {1..22}; do
    echo "Processing chromosome $chr..."
    pyrho optimize \
        --tablefile ${POP}_lookuptable.hdf \
        --vcffile CHM13-APGp1-HPRCp1-HGSVCp3_MC.chr${chr}.ex_tc.final.SNPs.biallelic.addtag.filter.vcf.gz \
        --outfile ${POP}_recombmap_chr${chr}_bp${BPEN}w${W} \
        -bpen $BPEN \
        -w $W \
        --numthreads $SLURM_CPUS_PER_TASK \
        --ploidy 2 \
        --logfile ${POP}_chr${chr}.log
done

echo "Job finished on $(hostname) at $(date)"

