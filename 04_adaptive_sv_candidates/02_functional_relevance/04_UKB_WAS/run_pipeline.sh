#!/usr/bin/env bash
# ==============================================================================
# UKB SV PheWAS and Biomarker Association Master Pipeline
# ==============================================================================
#
# Dependencies:
#   Python 3: pandas, numpy, statsmodels, scipy
#   R: PheWAS, dplyr, tidyr, optparse
#
# Before running this pipeline, load or activate the required software environment.
#
# Example:
#   conda activate <ENV_NAME>
#   module load R/<VERSION>
#
# ==============================================================================

set -euo pipefail

# ==============================================================================
# 1. Directory and environment configuration
# ==============================================================================

# Project directories
BASE_DIR="<PROJECT_ROOT>"
SCRIPT_DIR="<SCRIPT_DIRECTORY>"
RESULT_DIR="<OUTPUT_DIRECTORY>"

# Input files
UKB_VCF="<PATH_TO_UKB_SV_VCF_GZ>"
TRUTH_VCF="<PATH_TO_TRUTH_SV_VCF_GZ>"
UKB_PHENOTYPES="<PATH_TO_UKB_ALL_PHENOTYPES_TSV>"
UKB_PROTEOMICS="<PATH_TO_UKB_ALL_PROTEOMICS_TSV>"
ICD10_DATA="<PATH_TO_ICD10_LONG_FORMAT_CSV>"

# Computational resources
CPU_CORES="<NUMBER_OF_CPU_CORES>"

# ==============================================================================
# 2. Create output directories
# ==============================================================================

mkdir -p "${RESULT_DIR}/00.filtered_svs"
mkdir -p "${RESULT_DIR}/01.genotypes"
mkdir -p "${RESULT_DIR}/02.biomarker_assoc"
mkdir -p "${RESULT_DIR}/03.phewas_data"
mkdir -p "${RESULT_DIR}/04.phewas"

echo "============================================================"
echo " Starting UKB SV PheWAS and Biomarker Pipeline"
echo " Date: $(date)"
echo "============================================================"

# ==============================================================================
# Step 0: Filter UKB SVs against the truth SV set
# Criterion: breakpoint distance within ±1 kb and matched SVTYPE
# ==============================================================================

echo -e "\n---> Step 0: Filtering UKB SVs against truth SVs"

python "${SCRIPT_DIR}/00.filter_ukb_svs.py" \
    --truth "${TRUTH_VCF}" \
    --ukb "${UKB_VCF}" \
    --outdir "${RESULT_DIR}/00.filtered_svs"

VALID_SVS="${RESULT_DIR}/00.filtered_svs/valid_ukb_sv_ids.txt"

# ==============================================================================
# Step 1: Extract SV genotypes
# Genotypes are encoded as 0/1/2 and filtered by MAF >= 0.05
# using QC-passed samples defined in the phenotype table
# ==============================================================================

echo -e "\n---> Step 1: Extracting genotypes for valid SVs"

python "${SCRIPT_DIR}/01.extract_genotypes.py" \
    --vcf "${UKB_VCF}" \
    --sv-list "${VALID_SVS}" \
    --phenotype "${UKB_PHENOTYPES}" \
    --outdir "${RESULT_DIR}/01.genotypes" \
    --min-maf 0.05

GENOTYPE_MAT="${RESULT_DIR}/01.genotypes/sv_genotype_matrix.tsv.gz"

# ==============================================================================
# Step 2, Step 3 and Step 4:
# Run biomarker association and PheWAS workflows in parallel
# ==============================================================================

echo -e "\n---> Starting Step 2 biomarker association and Step 3/4 PheWAS in parallel"

(
    echo "[Step 2] Running SV-by-biomarker association..."

    python "${SCRIPT_DIR}/02.sv_biomarker_assoc.py" \
        --genotype "${GENOTYPE_MAT}" \
        --phenotype "${UKB_PHENOTYPES}" \
        --proteomics "${UKB_PROTEOMICS}" \
        --outdir "${RESULT_DIR}/02.biomarker_assoc" \
        --cpu "${CPU_CORES}"

    echo "[Step 2] Finished."
) &

(
    echo "[Step 3] Preparing PheWAS input data..."

    python "${SCRIPT_DIR}/03.prepare_phewas_data.py" \
        --genotype "${GENOTYPE_MAT}" \
        --phenotype "${UKB_PHENOTYPES}" \
        --icd10 "${ICD10_DATA}" \
        --outdir "${RESULT_DIR}/03.phewas_data"

    PHEWAS_GENO_COV="${RESULT_DIR}/03.phewas_data/phewas_covariates_genotypes.tsv"
    PHEWAS_ICD10="${RESULT_DIR}/03.phewas_data/filtered_icd10.csv"

    echo "[Step 4] Running SV PheWAS in R..."

    Rscript "${SCRIPT_DIR}/04.sv_phewas.R" \
        --geno_cov "${PHEWAS_GENO_COV}" \
        --icd10 "${PHEWAS_ICD10}" \
        --outdir "${RESULT_DIR}/04.phewas"

    echo "[Step 4] Finished."
) &

# Wait for all background jobs to finish
wait

echo "============================================================"
echo " Pipeline finished successfully"
echo " Date: $(date)"
echo "============================================================"
