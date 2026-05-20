#!/usr/bin/env python3
"""
Step 3: Prepare Data for PheWAS in R
============================================
Purpose: Integrate SV genotypes and QC-passed covariates for the R PheWAS package.
Operations:
  1. Read sv_genotype_matrix.tsv.gz.
  2. Read UKB.all.phenotypes.tsv and apply strict filtering.
  3. Merge genotypes with covariates (Age, Sex, PC1-10).
  4. Filter icd10_long_format.csv to retain only QC-passed samples.
"""

import argparse
import os
import pandas as pd
import time

def main():
    parser = argparse.ArgumentParser(description='Prepare PheWAS Covariates and Genotypes')
    parser.add_argument('--genotype', required=True, help='SV genotype matrix (sv_genotype_matrix.tsv.gz)')
    parser.add_argument('--phenotype', required=True, help='UKB.all.phenotypes.tsv')
    parser.add_argument('--icd10', required=True, help='icd10_long_format.csv')
    parser.add_argument('--outdir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load SV genotypes
    print(f"[{time.strftime('%H:%M:%S')}] Loading SV genotypes...")
    df_gt = pd.read_csv(args.genotype, sep='\t', compression='gzip')
    print(f"  Genotype shape: {df_gt.shape}")

    # 2. Load phenotype data and apply strict QC
    print(f"[{time.strftime('%H:%M:%S')}] Loading phenotype data and applying QC...")
    df_pheno = pd.read_csv(args.phenotype, sep='\t', usecols=lambda c: '-0' in c or c == 'eid')
    
    # Strict QC filtering
    df_pheno = df_pheno[
        (df_pheno['31-0.0'] == df_pheno['22001-0.0']) &
        (df_pheno['30079-0.0'] == 5) &
        (df_pheno['22006-0.0'] == 1) &
        (df_pheno['22019-0.0'].isnull()) &
        (df_pheno['22020-0.0'] == 1) &
        (df_pheno['22021-0.0'] != 10) &
        (df_pheno['22027-0.0'].isnull())
    ]
    
    # Rename covariates
    rename_dict = {'21022-0.0': 'Age', '31-0.0': 'Sex'}
    for i in range(1, 11):
        rename_dict[f'22009-0.{i}'] = f'PC{i}'
    df_pheno.rename(columns=rename_dict, inplace=True)
    
    covariates = ['eid', 'Age', 'Sex'] + [f'PC{i}' for i in range(1, 11)]
    df_cov = df_pheno[covariates].dropna()
    print(f"  Samples with complete QC and covariates: {df_cov.shape[0]}")

    # 3. Merge and export the file for R
    print(f"[{time.strftime('%H:%M:%S')}] Merging genotypes and covariates...")
    # Rename eid to id to match the default PheWAS convention
    df_cov.rename(columns={'eid': 'id'}, inplace=True)
    df_gt.rename(columns={'eid': 'id'}, inplace=True)
    
    df_merged = pd.merge(df_gt, df_cov, on='id', how='inner')
    out_merged = os.path.join(args.outdir, 'phewas_covariates_genotypes.tsv')
    df_merged.to_csv(out_merged, sep='\t', index=False)
    print(f"  → Saved: {out_merged} (N = {df_merged.shape[0]})")

    # 4. Filter ICD10 long-format data
    print(f"[{time.strftime('%H:%M:%S')}] Loading and filtering ICD10 data...")
    df_icd = pd.read_csv(args.icd10)
    # df_icd is expected to contain id, vocabulary_id and code by default
    valid_ids = set(df_merged['id'])
    df_icd_filtered = df_icd[df_icd['id'].isin(valid_ids)]
    
    out_icd = os.path.join(args.outdir, 'filtered_icd10.csv')
    df_icd_filtered.to_csv(out_icd, index=False)
    print(f"  → Saved: {out_icd} (rows: {df_icd_filtered.shape[0]})")
    
    print(f"[{time.strftime('%H:%M:%S')}] Done!")

if __name__ == '__main__':
    main()
