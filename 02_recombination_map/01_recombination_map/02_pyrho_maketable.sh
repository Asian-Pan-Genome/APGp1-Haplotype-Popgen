set -euo pipefail

echo "Job started on $(hostname) at $(date)"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "cpu: $SLURM_CPUS_PER_TASK"
source /share/home/zhanglab/user/weixinzhu/tools/miniconda3/etc/profile.d/conda.sh
conda activate pyrho

pop=CHA
n=198
N=300

pyrho make_table -n ${n} -N ${N} --mu 1.25e-8 --logfile ${pop}_make_table.log --outfile ${pop}_lookuptable.hdf --approx --smcpp_file ${pop}.csv --decimate_rel_tol 0.1 --numthreads 12

echo "Job finished on $(hostname) at $(date)"