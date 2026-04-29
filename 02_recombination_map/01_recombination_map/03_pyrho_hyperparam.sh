pop=CHA
n=198
N=300
pyrho hyperparam --samplesize ${n} --tablefile ${pop}_lookuptable.hdf \
--mu 1.25e-8 --num_sims 300 --logfile ${pop}_hyperparam.log --ploidy 2 \
--smcpp_file ${pop}.csv --outfile ${pop}_hyperparam_results.txt \
--numthreads 20