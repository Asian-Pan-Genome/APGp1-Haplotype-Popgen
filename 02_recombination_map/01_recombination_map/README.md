# 01_recombination_map

Build a fine-scale recombination map for the CHA cohort from the pangenome VCF.

## Scripts

- `01_smcpp.sh` — prepare per-chromosome SMC++ inputs from the VCFs and run SMC++ to estimate demographic history.
- `02_pyrho_maketable.sh` — build the pyrho two-locus likelihood lookup table using the SMC++ size history.
- `03_pyrho_hyperparam.sh` — run `pyrho hyperparam` to choose the `bpen`/`w` hyperparameters for the optimize step.
- `04_pyrho_optimize.sh` — run `pyrho optimize` on each chromosome to produce the final recombination map.

