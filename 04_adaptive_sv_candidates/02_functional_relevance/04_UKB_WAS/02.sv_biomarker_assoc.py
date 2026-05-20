#!/usr/bin/env python3
"""
Step 2: SV x Biomarker Association Analysis
============================================
Purpose: Run association analyses between SVs and biomarkers (proteomics and metabolomics).
Model: RINT(Biomarker) ~ SV_Genotype + Age + Sex + PC1-10
Testing: Apply multiple-testing correction separately for proteomics and metabolomics (Bonferroni and FDR).
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import argparse
import sys
import time
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests
import multiprocessing
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

def rint(series):
    """Rank-based inverse normal transformation"""
    valid_idx = series.notna()
    if valid_idx.sum() < 10:
        return pd.Series(np.nan, index=series.index)
    ranks = rankdata(series[valid_idx], method="average")
    ranks = (ranks - 0.5) / len(ranks)
    transformed = pd.Series(index=series.index, dtype=float)
    transformed.loc[valid_idx] = norm.ppf(ranks)
    return transformed

# Global variables shared by multiprocessing workers
_df_global = None
_sv_col_global = None
_covariates_global = None

def init_worker(df, sv_col, covariates):
    global _df_global, _sv_col_global, _covariates_global
    _df_global = df
    _sv_col_global = sv_col
    _covariates_global = covariates

def run_association(bm):
    global _df_global, _sv_col_global, _covariates_global
    df = _df_global
    sv_col = _sv_col_global
    covs = _covariates_global

    if bm not in df.columns: return None
    
    # Extract required columns and remove missing values
    tmp_df = df[['eid', sv_col] + covs + [bm]].dropna()
    if len(tmp_df) < 100 or tmp_df[sv_col].nunique() < 2:
        return None

    tmp_df = tmp_df.copy()
    tmp_df['bm_rint'] = rint(tmp_df[bm])
    tmp_df = tmp_df.dropna(subset=['bm_rint'])
    
    if len(tmp_df) < 100: return None

    try:
        tmp_df = tmp_df.rename(columns={sv_col: 'SV_geno'})
        formula = "bm_rint ~ SV_geno + Age + Sex + " + " + ".join([f"PC{i}" for i in range(1, 11)])
        model = smf.ols(formula, data=tmp_df).fit()
        return {
            'Biomarker': bm,
            'Beta': model.params['SV_geno'],
            'SE': model.bse['SV_geno'],
            'P_value': model.pvalues['SV_geno'],
            'N_samples': len(tmp_df)
        }
    except Exception:
        return None

def correct_pvals(df_res, label):
    """Independent Bonferroni and FDR correction"""
    if len(df_res) == 0: return df_res
    mask = df_res['P_value'].notna()
    if mask.sum() == 0: return df_res
    
    pvals = df_res.loc[mask, 'P_value']
    _, bonf, _, _ = multipletests(pvals, alpha=0.05, method='bonferroni')
    _, fdr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    
    df_res.loc[mask, 'Bonferroni_P'] = bonf
    df_res.loc[mask, 'FDR_P'] = fdr
    df_res['Sig_Bonferroni'] = df_res['Bonferroni_P'] < 0.05
    df_res['Sig_FDR'] = df_res['FDR_P'] < 0.05
    
    print(f"  [{label}] Number of tests: {len(pvals)}, Bonferroni-significant: {df_res['Sig_Bonferroni'].sum()}, FDR-significant: {df_res['Sig_FDR'].sum()}")
    return df_res

def main():
    parser = argparse.ArgumentParser(description='SV x Biomarker Association')
    parser.add_argument('--genotype', required=True, help='SV genotype matrix (sv_genotype_matrix.tsv.gz)')
    parser.add_argument('--phenotype', required=True, help='UKB.all.phenotypes.tsv')
    parser.add_argument('--proteomics', required=True, help='UKB.all.proteomics.tsv')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--cpu', type=int, default=16, help='Number of CPUs')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] 1. Loading genotype matrix")
    df_gt = pd.read_csv(args.genotype, sep='\t', compression='gzip')
    sv_cols = [c for c in df_gt.columns if c != 'eid']
    print(f"  Samples: {df_gt.shape[0]}, SVs: {len(sv_cols)}")

    print(f"[{time.strftime('%H:%M:%S')}] 2. Loading phenotype data and applying strict QC")
    df_pheno = pd.read_csv(args.phenotype, sep='\t', usecols=lambda c: '-0' in c or c == 'eid')
    print(f"  Raw phenotype samples: {df_pheno.shape[0]}")
    
    # Strictly match the QC criteria used in the previous pipeline
    df_pheno = df_pheno[
        (df_pheno['31-0.0'] == df_pheno['22001-0.0']) &
        (df_pheno['30079-0.0'] == 5) &
        (df_pheno['22006-0.0'] == 1) &
        (df_pheno['22019-0.0'].isnull()) &
        (df_pheno['22020-0.0'] == 1) &
        (df_pheno['22021-0.0'] != 10) &
        (df_pheno['22027-0.0'].isnull())
    ]
    print(f"  QC-passed samples: {df_pheno.shape[0]}")

    print(f"[{time.strftime('%H:%M:%S')}] 3. Loading proteomics and metabolomics data")
    df_prot = pd.read_csv(args.proteomics, sep='\t')
    prot_all_cols = [c for c in df_prot.columns if c != 'eid']
    
    # Merge tables
    df = pd.merge(df_gt, df_prot, on='eid', how='inner')
    df = pd.merge(df, df_pheno, on='eid', how='inner')
    print(f"  Samples with complete merged data: {df.shape[0]}")

    # Rename covariates
    rename_dict = {'21022-0.0': 'Age', '31-0.0': 'Sex'}
    for i in range(1, 11):
        rename_dict[f'22009-0.{i}'] = f'PC{i}'
    df.rename(columns=rename_dict, inplace=True)
    
    covariates = ['Age', 'Sex'] + [f'PC{i}' for i in range(1, 11)]
    df = df.dropna(subset=covariates)
    
    prot_cols = [c for c in prot_all_cols if c in df.columns]
    metab_cols = [c for c in df.columns if (c.startswith('30') or c.startswith('23') or c.startswith('20')) and '-0.0' in c and c != '30079-0.0']
    
    print(f"  Identified proteomics traits: {len(prot_cols)} items; metabolomics traits: {len(metab_cols)} items")

    # 4. Run association analyses
    res_prot = []
    res_metab = []
    
    for i, sv in enumerate(sv_cols):
        print(f"\n[{time.strftime('%H:%M:%S')}] Analyzing SV: {sv} ({i+1}/{len(sv_cols)})")
        
        # Analyze proteomics traits in parallel
        if prot_cols:
            with multiprocessing.Pool(args.cpu, init_worker, (df, sv, covariates)) as pool:
                cur_prot = pool.map(run_association, prot_cols)
            cur_prot = [r for r in cur_prot if r]
            if cur_prot:
                pdf = pd.DataFrame(cur_prot)
                pdf.insert(0, 'SV_ID', sv)
                res_prot.append(pdf)
                
        # Analyze metabolomics traits in parallel
        if metab_cols:
            with multiprocessing.Pool(args.cpu, init_worker, (df, sv, covariates)) as pool:
                cur_metab = pool.map(run_association, metab_cols)
            cur_metab = [r for r in cur_metab if r]
            if cur_metab:
                mdf = pd.DataFrame(cur_metab)
                mdf.insert(0, 'SV_ID', sv)
                res_metab.append(mdf)
                
    # 5. Apply multiple-testing correction and save results
    print(f"\n[{time.strftime('%H:%M:%S')}] Applying stratified multiple-testing correction...")
    
    if res_prot:
        df_all_prot = pd.concat(res_prot, ignore_index=True)
        df_all_prot = correct_pvals(df_all_prot, "Proteomics")
        out_prot = os.path.join(args.outdir, 'sv_proteomics_results.tsv')
        df_all_prot.sort_values('P_value').to_csv(out_prot, sep='\t', index=False)
        print(f"  Saved -> {out_prot}")
        
    if res_metab:
        df_all_metab = pd.concat(res_metab, ignore_index=True)
        df_all_metab = correct_pvals(df_all_metab, "Metabolomics")
        out_metab = os.path.join(args.outdir, 'sv_metabolomics_results.tsv')
        df_all_metab.sort_values('P_value').to_csv(out_metab, sep='\t', index=False)
        print(f"  Saved -> {out_metab}")

    print(f"[{time.strftime('%H:%M:%S')}] Done!")

if __name__ == '__main__':
    main()
