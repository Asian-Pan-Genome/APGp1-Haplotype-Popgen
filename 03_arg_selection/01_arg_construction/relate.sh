#!/bin/bash


PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"


Relate_path="/path/to/relate"
RELATE_BIN="${Relate_path}/bin/Relate"
RELATE_FILEFORMATS="${Relate_path}/bin/RelateFileFormats"
RELATE_CONVERT="${Relate_path}/relate_lib/bin/Convert"
RELATE_SCRIPTS="${Relate_path}/scripts"

DATA_DIR="${PROJECT_ROOT}/data/relate"
RESULT_DIR="${PROJECT_ROOT}/results/relate"
TS_DIR="${PROJECT_ROOT}/results/relate_ts"

mkdir -p "${RESULT_DIR}" "${TS_DIR}"


CHR=6
MUT=1.25e-8
N=40000
SEED=100

#for sr data
PREFIX="new_sr_full_chr${CHR}.polarized"
#PREFIX="new_lr_full_chr${CHR}.polarized"

HAPS="${DATA_DIR}/${PREFIX}.haps"
SAMPLE="${DATA_DIR}/${PREFIX}.sample"
INPUT_VCF="${DATA_DIR}/input.vcf"

MAP="${DATA_DIR}/map/chm13_chr${CHR}.map"
POPLABELS="${DATA_DIR}/${PREFIX}.poplabels"

OUT_PREFIX="${RESULT_DIR}/${PREFIX}"

#convert vcf to relate format
echo "ConvertFromVcf"
"${RELATE_FILEFORMATS}" \
    --mode ConvertFromVcf \
    --haps "${HAPS}" \
    --sample "${SAMPLE}" \
    -i "${INPUT_VCF}"


#run relate
echo "run relate"
cd "${RESULT_DIR}"
"${RELATE_BIN}" \
    --mode All \
    -m "${MUT}" \
    -N "${N}" \
    --haps "${HAPS}" \
    --sample "${SAMPLE}" \
    --map "${MAP}" \
    --seed "${SEED}" \
    -o "${PREFIX}"

#convert to ts
echo "ConvertToTreeSequence"

"${RELATE_CONVERT}" \
    --mode ConvertToTreeSequence \
    --anc "${OUT_PREFIX}.anc" \
    --mut "${OUT_PREFIX}.mut" \
    -o "${TS_DIR}/${PREFIX}"

# optional: estimate population size (and later do branch MCMC sampling)
echo "EstimatePopulationSize"
bash "${RELATE_SCRIPTS}/EstimatePopulationSize/EstimatePopulationSize.sh" \
    -i "${OUT_PREFIX}" \
    -m "${MUT}" \
    --poplabels "${POPLABELS}" \
    --seed "${SEED}" \
    -o "${PREFIX}" \
    --noplot

#branch MCMC sampling
echo "SampleBranchLengths"
bash "${RELATE_SCRIPTS}/SampleBranchLengths/ReEstimateBranchLengths.sh" \
    -i "${OUT_PREFIX}" \
    -o "${OUT_PREFIX}.re.branch" \
    -m "${MUT}" \
    --coal "${OUT_PREFIX}.coal" \
    --seed "${SEED}"

#relate ts after branch MCMC
echo "re-estimate ConvertToTreeSequence"
"${RELATE_CONVERT}" \
    --mode ConvertToTreeSequence \
    --anc "${OUT_PREFIX}.re.branch.anc" \
    --mut "${OUT_PREFIX}.re.branch.mut" \
    -o "${TS_DIR}/${PREFIX}.re.branch"

echo "done"