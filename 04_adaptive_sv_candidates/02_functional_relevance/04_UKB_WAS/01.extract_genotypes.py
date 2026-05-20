#!/usr/bin/env python3
"""
Step 1: Extract SV genotype matrix from VCF
============================================
Purpose: Extract the genotype matrix for valid SVs from Step 0 and apply MAF filtering using QC-passed samples.
Logic:
  1. Read UKB.all.phenotypes.tsv and apply strict QC.
  2. Obtain the EID list for QC-passed samples (~330,000 individuals).
  3. Scan the VCF and retain only columns matching QC-passed samples.
  4. Calculate MAF in the QC-passed cohort and remove SVs with MAF < 5%.
  5. Export sv_genotype_matrix.tsv.gz containing only QC-passed samples.
"""

import argparse
import gzip
import os
import sys
import time
import numpy as np
import pandas as pd

def infer_svtype_and_len(ref, alt, info_str):
    svtype = "UNKNOWN"
    svlen = 0
    
    # Extract SVTYPE and SVLEN from INFO if available
    for field in info_str.split(';'):
        if field.startswith('SVTYPE='):
            svtype = field.split('=')[1].upper()
        elif field.startswith('SVLEN='):
            try:
                svlen = abs(int(field.split('=')[1]))
            except ValueError:
                pass

    first_alt = alt.split(',')[0]
    
    # Infer SVTYPE if missing
    if svtype == "UNKNOWN":
        if first_alt.startswith('<') and first_alt.endswith('>'):
            svtype = first_alt[1:-1].upper()
        elif first_alt != '.':
            if len(ref) > len(first_alt):
                svtype = "DEL"
            elif len(ref) < len(first_alt):
                svtype = "INS"
            else:
                svtype = "SNP_MNP"
                
    # Infer SVLEN if missing or 0
    if svlen == 0 and not first_alt.startswith('<') and first_alt != '.':
        svlen = abs(len(first_alt) - len(ref))
        
    return svtype, svlen

def parse_gt(gt_str):
    if gt_str in ('./.', '.|.', '.'): return 0
    sep = '|' if '|' in gt_str else '/'
    alleles = gt_str.split(sep)
    if len(alleles) != 2: return 0
    
    a1 = int(alleles[0]) if alleles[0] != '.' else 0
    a2 = int(alleles[1]) if alleles[1] != '.' else 0
    return a1 + a2

def get_qc_passed_eids(pheno_path):
    print(f"[{time.strftime('%H:%M:%S')}] Reading phenotype data and applying QC: {pheno_path}")
    df_pheno = pd.read_csv(pheno_path, sep='\t', usecols=lambda c: '-0' in c or c == 'eid')
    print(f"  Raw samples: {df_pheno.shape[0]}")
    
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
    return set(df_pheno['eid'].astype(int))

def main():
    parser = argparse.ArgumentParser(description='Extract SV genotype matrix for QC-passed samples')
    parser.add_argument('--vcf', required=True, help='UKB collapsed VCF (bgzipped)')
    parser.add_argument('--sv-list', required=True, help='valid_ukb_sv_ids.txt from Step 0')
    parser.add_argument('--phenotype', required=True, help='UKB.all.phenotypes.tsv for strict QC')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--min-maf', type=float, default=0.05, help='Minimum MAF in QC-passed cohort (default: 0.05)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Get the EID set for valid samples
    qc_eids = get_qc_passed_eids(args.phenotype)
    print('Total pass qc samples: {len(qc_eids)}')
    if not qc_eids:
        print("No QC-passed samples found. Exiting.")
        sys.exit(1)

    # 2. Read the valid SV list
    print(f"[{time.strftime('%H:%M:%S')}] Loading target SV list: {args.sv_list}")
    with open(args.sv_list) as f:
        target_svs = set(line.strip() for line in f if line.strip())
    print(f"  Number of SVs to extract: {len(target_svs)}")

    # 3. Scan the VCF
    print(f"[{time.strftime('%H:%M:%S')}] Scanning VCF and applying MAF filtering: {args.vcf}")
    sv_ids = []
    genotypes = []
    valid_indices = []
    valid_sample_names = []

    with gzip.open(args.vcf, 'rt') as f:
        for line in f:
            if line.startswith('##'): continue
            
            # Parse the header line
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                all_samples = cols[9:]
                # Find indices of QC-passed samples in the VCF (indices are relative to the sample fields)
                for i, sname in enumerate(all_samples):
                    if int(sname) in qc_eids:
                        valid_indices.append(i)
                        valid_sample_names.append(sname)
                print(f"  VCF contains {len(all_samples)} samples; matched QC-passed samples: {len(valid_indices)} ")
                continue

            # Parse data rows using a fast preliminary split
            fields_head = line.strip().split('\t', 8)
            if len(fields_head) < 9: continue
            
            sv_id = fields_head[2]
            if sv_id not in target_svs: continue

            # Target SV matched
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            svtype, svlen = infer_svtype_and_len(ref, alt, info)
            new_sv_id = f"{chrom}:{pos}:{svtype}:{svlen}"
            
            fmt_fields = fields[8].split(':')
            gt_idx = fmt_fields.index('GT') if 'GT' in fmt_fields else 0

            # Extract GT only for QC-passed samples
            gt_array = np.zeros(len(valid_indices), dtype=np.int8)
            for out_idx, in_idx in enumerate(valid_indices):
                gt_str = fields[9 + in_idx].split(':')[gt_idx]
                gt_array[out_idx] = parse_gt(gt_str)
            
            # Calculate MAF in this QC-passed sample set
            af = gt_array.sum() / (2.0 * len(valid_indices))
            maf = min(af, 1 - af)

            if maf >= args.min_maf:
                sv_ids.append(new_sv_id)
                genotypes.append(gt_array)
            #else:
                #print(f"  SV {new_sv_id} (original: {sv_id}) in QC-passed samples MAF={maf:.4f} < {args.min_maf}，filtered.")

    print(f"[{time.strftime('%H:%M:%S')}] VCF scan completed. Retained {len(sv_ids)} SVs passing filters.")
    if not sv_ids:
        print("No SVs passed the filters. Exiting.")
        sys.exit(0)

    # 4. Save the output file
    print(f"[{time.strftime('%H:%M:%S')}] Building and saving genotype matrix...")
    gt_matrix = np.array(genotypes).T
    df_gt = pd.DataFrame(gt_matrix, columns=sv_ids)
    df_gt.insert(0, 'eid', [int(s) for s in valid_sample_names])

    out_file = os.path.join(args.outdir, 'sv_genotype_matrix.tsv.gz')
    df_gt.to_csv(out_file, sep='\t', index=False, compression='gzip')
    print(f"  → Matrix saved to {out_file} (dimensions: {df_gt.shape[0]} rows x {df_gt.shape[1]-1} SV)")

if __name__ == '__main__':
    main()
