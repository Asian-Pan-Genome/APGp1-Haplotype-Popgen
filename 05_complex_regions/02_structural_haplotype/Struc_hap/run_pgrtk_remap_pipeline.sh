#!/usr/bin/env bash

set -euo pipefail

# Independent workflow: original PGR-TK bundles -> ranked iterative bundle library ->
# minimap2 remapping -> new PGR-compatible BED for downstream graph merging.
# No graph merging is performed here; optional TRF filtering/refinement precedes additions.

usage() {
    cat <<'EOF'
Usage:
  bash run_pgrtk_remap_pipeline.sh \
    --pgr-bed ORIGINAL_PGRTK.bed \
    --fasta POPULATION.fa \
    --output-dir OUTPUT_DIR [options]

Required:
  --pgr-bed FILE              Original four-column PGR-TK BED; supplies bundle
                              occurrences for consensus and original IDs/sizes.
  --fasta FILE                Population FASTA with the same path IDs as BED.
                              Repeated FASTA IDs are excluded entirely, along
                              with their BED/TRF/sample-map entries.
  --output-dir DIR            Separate output directory; inputs are not updated.

Options:
  --trf-bed FILE              Optional TRF annotation; empty/comment-only skips
                              all TRF operations. Omitting it at final-only
                              step 3 disables the final TRF filter.
  --trf-min-length BP         Minimum annotation length (default: 0; no filter).
  --trf-min-identity PCT      Minimum annotation identity percent (default: 80).
  --trf-full-coverage F       Initial single-merged-TRF-interval occurrence
                              coverage threshold (default: 0.90).
  --final-trf-full-coverage F Final-library single-merged-TRF-interval
                              occurrence coverage threshold (default: 0.95).
  --trf-bundle-removal-frequency F
                              Initial occurrence and final carrier-sequence
                              removal threshold (default: 0.50).
  --trf-split-min-bp BP       Long divergent repeat threshold (default: 1000).
  --trf-min-sequence-coverage F
                              Repeat presence coverage threshold (default: 0.50).
  --trf-min-population-presence F
                              Structural divergence threshold (default: 0.50).
  --trf-min-length-divergence F
                              Mean pairwise TR-length divergence threshold
                              (default: 0.05).
  --trf-min-nucleotide-divergence F
                              Nucleotide divergence threshold (default: 0.05).
  --trf-min-fragment-sequence-bp BP
                              Real bases required per child sequence (default: 1000).
  --trf-min-fragment-count INT Minimum retained child support (default: 2).
  --trf-min-fragment-frequency F
                              Minimum retained child frequency (default: 0.02).
  --start-step INT            Start at 1--3 and run through step 3 (default: 1).
                              Reuse earlier outputs under --output-dir.
  --sample-map FILE           Two-column path_id/sample_id map. Rank by distinct
                              samples; without a map, each path is one sample.
  --threads INT               MAFFT/minimap2 threads (default: 8).
  --mafft PATH                MAFFT executable (default: mafft; steps 2 and 3).
  --minimap2 PATH             minimap2 executable (default: minimap2).
  --needle PATH               EMBOSS needle executable used by the post-filter
                              overlap refinement (default: needle).
  --no-post-filter-needle-refinement
                              Disable the default 03-stage Needle confirmation,
                              split/reconcile, remap, and refilter loop.
  --mafft-extra STRING        Additional quoted MAFFT arguments in both stages.
  --mafft-02-params STRING    MAFFT strategy/parameters for 02_bundle_consensus
                              (default: --auto; replaces --auto).
  --mafft-02-quick            Window-align every source FASTA. Sequences at least
                              70% of the longest occurrence form the backbone;
                              shorter sequences join matching windows locally.
                              Skip truncated records <=10% of that reference
                              AND <1000 bp before alignment.
                              Off by default.
  --mafft-02-window-bp BP     Reference window width in 02 quick mode
                              (default: 10000).
  --mafft-03-params STRING    MAFFT strategy/parameters for the 03 common-indel
                              MSA (default: --auto; replaces --auto).
  --resume                    Reuse a completed iterative library only when
                              its input/parameter fingerprint is unchanged.
  --resume-ranks              Resume after the last completed rank checkpoint;
                              starts at step 2 and reruns any interrupted rank.
  --keep-intermediates        Keep disposable MAFFT/rank FASTAs and every self
                              snapshot; default cleans successful-stage temporaries.
  --min-column-occupancy F    Consensus column occupancy (default: 0.50).
  --min-base-frequency F      Consensus majority-base frequency (default: 0.50).
  --merge-identity F          Shared overlap identity >F: one representative;
                              otherwise retain same-class alleles (default: 0.90).
                              Also used for same-ID and common-deletion children.
  --addition-min-bp BP        Minimum uncovered fragment to add to the library
                              (default: 1; positive integer).
  --rank1-self-min-overlap-fraction F
                              Split self-overlaps at rank 1 and after additions
                              and reconcile again after the final split, exceeding either
                              bundle's length fraction (default: 0.20).
  --rank1-self-min-overlap-bp BP
                              Also split shared spans exceeding BP at all
                              self-alignment stages (default: 10000).
                              Fixed self terminal guard: suppress an outer cut
                              only if fragment <200 bp AND <5% of current
                              parent length; internal pieces have no minimum.
  --no-rank1-self-split       Disable iterative rank-1 bundle self-alignment
                              splitting; independent of --no-mapping-split.
  --no-addition-self-split    Disable iterative library self-splitting after new
                              lower-rank additions; other split stages unchanged.
  --no-self-overlap-reconciliation
                              Legacy self cuts; do not reconcile shared regions,
                              including after the final split.
  --min-bundle-coverage F     Chained alignment bundle coverage (default: 0.60).
  --min-identity F            Minimum PAF block identity (default: 0.0).
  --min-mapq INT              Minimum PAF block MAPQ (default: 0).
  --max-bundle-gap BP         Maximum collinear bundle gap/overlap (default: 1000).
  --max-genome-gap BP         Maximum collinear genomic gap/overlap (default: 1000).
  --max-overlap-fraction F   Maximum overlap on either axis / bundle consensus
                              length (default: 0.10).
  --fixed-bundle-min-coverage F
                              Single-PAF coverage to fix a block and exclude
                              it from chaining (default: 0.80).
  --locus-dedup-min-bundle-coverage F
                              Minimum completed-chain union bundle coverage
                              admitted to final remapping (default: 0.70).
  --max-fixed-overlap-fraction F
                              Maximum candidate-span overlap / fixed bundle
                              consensus length (default: 0.20).
  --max-fixed-overlap-bp BP   Maximum overlap with a different fixed class
                              (default: 10000).
  --mapping-split-min-truncated-count INT
                              Gap-supported short-occurrence count (default: 5).
  --mapping-split-truncated-fraction F
                              Strict truncation fraction (default: 0.60).
  --mapping-split-max-unmapped-bp BP
                              Long unaligned bundle report threshold (default: 10000).
  --mapping-split-cut-merge-bp BP
                              Merge nearby candidate cuts (default: 100).
  --mapping-split-min-cut-support INT
                              Short-occurrence support per cut (default: 1).
  --mapping-split-min-terminal-bp BP
                              Genome-remapping-only terminal fragment guard;
                              does not apply to bundle self-alignment.
                              (default: 1000; 0 disables; internal pieces retained).
  --mapping-split-gap-min-mapping-bp BP
                              Selected PAF span in a preliminary-BED gap must
                              exceed BP (default: 10000).
  --mapping-split-gap-max-overlap-bp BP
                              Selected pairs must overlap by less than BP
                              inside the gap (default: 1000; 0: no overlap).
                              Must not exceed --mapping-split-gap-min-mapping-bp.
  --no-mapping-split          Disable mapping-driven splitting in ranked rounds
                              and final step 3; addition overlap clipping remains.
  --common-deletion-min-bp BP Final-MSA internal deletion/insertion length must
                              be strictly greater than BP (default: 1000).
  --common-deletion-min-count INT
                              Minimum distinct samples (with --sample-map), else
                              paths, supporting a common indel (default: 5).
  --common-deletion-min-frequency F
                              Minimum mapped population-unit frequency supporting
                              a common indel (default: 0.01).
  --common-deletion-event-overlap F
                              Reciprocal overlap for clustering deletion events
                              (default: 0.80).
  --no-common-deletion-split  Skip final occurrence MSA/indel splitting.
  -h, --help                  Show this help.

Steps:
  1  Optionally TRF-filter occurrences, group signed bundle orders and rank
     haplotypes by distinct sample count.
  2  Initialize from rank 1 and iteratively split bundle self-overlaps; then
     map to each next rank, split supported library bundles, and add
     optional TRF-refined consensuses of uncovered source-bundle fragments,
     then self-split the updated library before advancing to the next rank.
  3  Remap to all population paths, apply mapping-driven splitting, build a
     preliminary locus-deduplicated BED, detect common internal indels by
     MAFFT, self-align the final split library, reconcile Needle-confirmed
     overlaps, apply the final per-sequence TRF filter, remap any corrected
     library, and publish pgrtk_remapped.bed.

This is an initial PGR-TK refinement step, before downstream graph merging.
No previous merged BED is required and no novelty comparison is performed.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pgr_bed=""
population_fasta=""
trf_bed=""
trf_min_length=0
trf_min_identity=80
trf_full_coverage=0.90
final_trf_full_coverage=0.95
trf_bundle_removal_frequency=0.50
trf_split_min_bp=1000
trf_min_sequence_coverage=0.50
trf_min_population_presence=0.50
trf_min_length_divergence=0.05
trf_min_nucleotide_divergence=0.05
trf_min_fragment_sequence_bp=1000
trf_min_fragment_count=2
trf_min_fragment_frequency=0.02
sample_map=""
output_dir=""
start_step=1
threads=8
mafft_cmd="mafft"
minimap2_cmd="minimap2"
needle_cmd="needle"
post_filter_needle_refinement=1
mafft_extra=""
mafft_02_params="--auto"
mafft_03_params="--auto"
mafft_02_quick=0
mafft_02_window_bp=10000
resume=0
resume_ranks=0
min_column_occupancy=0.50
min_base_frequency=0.50
merge_identity=0.90
addition_min_bp=1
rank1_self_min_overlap_fraction=0.20
rank1_self_min_overlap_bp=10000
rank1_self_split_args=()
addition_self_split_args=()
self_overlap_reconciliation_args=()
intermediate_args=()
min_bundle_coverage=0.60
min_identity=0.0
min_mapq=0
max_bundle_gap=1000
max_genome_gap=1000
max_overlap_fraction=0.10
fixed_bundle_min_coverage=0.80
locus_dedup_min_bundle_coverage=0.70
max_fixed_overlap_fraction=0.20
max_fixed_overlap_bp=10000
split_gap_min_mapping_bp=10000
split_gap_max_overlap_bp=1000
split_min_terminal_bp=1000
mapping_split_args=()
common_deletion_split=1
common_deletion_min_bp=1000
common_deletion_min_count=5
common_deletion_min_frequency=0.01
common_deletion_event_overlap=0.80

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pgr-bed) pgr_bed="${2:?missing --pgr-bed value}"; shift 2 ;;
        --fasta) population_fasta="${2:?missing --fasta value}"; shift 2 ;;
        --trf-bed) trf_bed="${2:?missing --trf-bed value}"; shift 2 ;;
        --trf-min-length|--trf-min-identity|--trf-full-coverage|--trf-bundle-removal-frequency|--trf-split-min-bp|--trf-min-sequence-coverage|--trf-min-population-presence|--trf-min-length-divergence|--trf-min-nucleotide-divergence|--trf-min-fragment-sequence-bp|--trf-min-fragment-count|--trf-min-fragment-frequency)
            trf_option="${1#--}"
            trf_option="${trf_option//-/_}"
            printf -v "$trf_option" '%s' "${2:?missing TRF option value}"
            shift 2 ;;
        --final-trf-full-coverage)
            final_trf_full_coverage="${2:?missing final TRF full-coverage value}"; shift 2 ;;
        --sample-map) sample_map="${2:?missing --sample-map value}"; shift 2 ;;
        --output-dir) output_dir="${2:?missing --output-dir value}"; shift 2 ;;
        --start-step) start_step="${2:?missing --start-step value}"; shift 2 ;;
        --threads) threads="${2:?missing --threads value}"; shift 2 ;;
        --mafft) mafft_cmd="${2:?missing --mafft value}"; shift 2 ;;
        --minimap2) minimap2_cmd="${2:?missing --minimap2 value}"; shift 2 ;;
        --mafft-extra) mafft_extra="${2:?missing --mafft-extra value}"; shift 2 ;;
        --mafft-02-params) mafft_02_params="${2:?missing --mafft-02-params value}"; shift 2 ;;
        --mafft-03-params) mafft_03_params="${2:?missing --mafft-03-params value}"; shift 2 ;;
        --mafft-02-quick) mafft_02_quick=1; shift ;;
        --mafft-02-window-bp) mafft_02_window_bp="${2:?missing --mafft-02-window-bp value}"; shift 2 ;;
        --needle) needle_cmd="${2:?missing --needle value}"; shift 2 ;;
        --no-post-filter-needle-refinement) post_filter_needle_refinement=0; shift ;;
        --min-column-occupancy) min_column_occupancy="${2:?missing value}"; shift 2 ;;
        --min-base-frequency) min_base_frequency="${2:?missing value}"; shift 2 ;;
        --merge-identity) merge_identity="${2:?missing value}"; shift 2 ;;
        --addition-min-bp) addition_min_bp="${2:?missing value}"; shift 2 ;;
        --rank1-self-min-overlap-fraction)
            rank1_self_min_overlap_fraction="${2:?missing rank1 self-split value}"; shift 2 ;;
        --rank1-self-min-overlap-bp)
            rank1_self_min_overlap_bp="${2:?missing rank1 self-split value}"; shift 2 ;;
        --no-rank1-self-split) rank1_self_split_args+=("$1"); shift ;;
        --no-addition-self-split) addition_self_split_args+=("$1"); shift ;;
        --no-self-overlap-reconciliation) self_overlap_reconciliation_args+=("$1"); shift ;;
        --keep-intermediates) intermediate_args+=("$1"); shift ;;
        --min-bundle-coverage) min_bundle_coverage="${2:?missing value}"; shift 2 ;;
        --min-identity) min_identity="${2:?missing value}"; shift 2 ;;
        --min-mapq) min_mapq="${2:?missing value}"; shift 2 ;;
        --max-bundle-gap) max_bundle_gap="${2:?missing value}"; shift 2 ;;
        --max-genome-gap) max_genome_gap="${2:?missing value}"; shift 2 ;;
        --max-overlap-fraction) max_overlap_fraction="${2:?missing value}"; shift 2 ;;
        --fixed-bundle-min-coverage) fixed_bundle_min_coverage="${2:?missing value}"; shift 2 ;;
        --locus-dedup-min-bundle-coverage) locus_dedup_min_bundle_coverage="${2:?missing value}"; shift 2 ;;
        --max-fixed-overlap-fraction) max_fixed_overlap_fraction="${2:?missing value}"; shift 2 ;;
        --max-fixed-overlap-bp) max_fixed_overlap_bp="${2:?missing value}"; shift 2 ;;
        --mapping-split-gap-min-mapping-bp)
            split_gap_min_mapping_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_gap_min_mapping_bp"); shift 2 ;;
        --mapping-split-gap-max-overlap-bp)
            split_gap_max_overlap_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_gap_max_overlap_bp"); shift 2 ;;
        --mapping-split-min-terminal-bp)
            split_min_terminal_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_min_terminal_bp"); shift 2 ;;
        --mapping-split-min-truncated-count|--mapping-split-truncated-fraction|--mapping-split-max-unmapped-bp|--mapping-split-cut-merge-bp|--mapping-split-min-cut-support)
            mapping_split_args+=("$1" "${2:?missing mapping-split value}"); shift 2 ;;
        --no-mapping-split) mapping_split_args+=("$1"); shift ;;
        --common-deletion-min-bp)
            common_deletion_min_bp="${2:?missing common-deletion value}"; shift 2 ;;
        --common-deletion-min-count)
            common_deletion_min_count="${2:?missing common-deletion value}"; shift 2 ;;
        --common-deletion-min-frequency)
            common_deletion_min_frequency="${2:?missing common-deletion value}"; shift 2 ;;
        --common-deletion-event-overlap)
            common_deletion_event_overlap="${2:?missing common-deletion value}"; shift 2 ;;
        --no-common-deletion-split) common_deletion_split=0; shift ;;
        --resume) resume=1; shift ;;
        --resume-ranks) resume_ranks=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "error: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$pgr_bed" || -z "$population_fasta" || -z "$output_dir" ]]; then
    echo "error: --pgr-bed, --fasta, and --output-dir are required" >&2
    usage >&2
    exit 2
fi
if [[ ! "$start_step" =~ ^[1-3]$ ]]; then
    echo "error: --start-step must be 1, 2, or 3" >&2
    exit 2
fi
if [[ "$resume_ranks" -eq 1 ]]; then
    if [[ "$start_step" -eq 3 ]]; then
        echo "error: --resume-ranks requires the step-2 library build" >&2
        exit 2
    fi
    start_step=2
fi
if [[ ! "$threads" =~ ^[0-9]+$ ]] || [[ "$threads" -le 0 ]]; then
    echo "error: --threads must be a positive integer" >&2
    exit 2
fi
if [[ ! "$mafft_02_window_bp" =~ ^[0-9]+$ ]] || [[ "$mafft_02_window_bp" -le 0 ]]; then
    echo "error: --mafft-02-window-bp must be a positive integer" >&2
    exit 2
fi

require_nonempty_file() {
    if [[ ! -f "$1" || ! -s "$1" ]]; then
        echo "error: $2 is missing or empty: $1" >&2
        exit 2
    fi
}
require_executable() {
    if [[ "$1" == */* ]]; then
        if [[ ! -f "$1" || ! -x "$1" ]]; then
            echo "error: $2 is not executable: $1" >&2
            exit 2
        fi
    elif ! command -v "$1" >/dev/null 2>&1; then
        echo "error: $2 was not found: $1" >&2
        exit 2
    fi
}

alignment_dir="$output_dir/01_bundle_alignments"
consensus_dir="$output_dir/02_bundle_consensus"
refine_dir="$output_dir/03_refined_bundle_mapping"
final_bed="$output_dir/pgrtk_remapped.bed"
library_script="$script_dir/scripts/build_ranked_bundle_library.py"
refine_script="$script_dir/scripts/refine_v2.sh"

require_nonempty_file "$pgr_bed" "original PGR BED"
require_nonempty_file "$population_fasta" "population FASTA"
if [[ -n "$trf_bed" && ! -f "$trf_bed" ]]; then
    echo "error: --trf-bed file is missing: $trf_bed" >&2
    exit 2
fi
if [[ -n "$sample_map" ]]; then require_nonempty_file "$sample_map" "sample map"; fi
require_executable python3 Python3
require_executable "$minimap2_cmd" minimap2
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    require_executable "$needle_cmd" "EMBOSS needle"
fi
if [[ "$start_step" -le 2 || "$common_deletion_split" -eq 1 ]]; then
    require_executable "$mafft_cmd" MAFFT
fi
if [[ "$start_step" -eq 2 ]]; then
    require_nonempty_file "$alignment_dir/ranking_summary.json" "step-2 haplotype ranking"
elif [[ "$start_step" -eq 3 ]]; then
    require_nonempty_file "$consensus_dir/bundle.fa" "step-3 ranked bundle library"
    require_nonempty_file "$consensus_dir/library_summary.json" "step-3 completed library audit; rerun steps 1/2 for an old library"
    require_nonempty_file "$consensus_dir/bundle_metadata.bed" "step-3 ranked library metadata"
fi
require_nonempty_file "$library_script" "ranked library helper/audit checker"
require_nonempty_file "$refine_script" "refinement workflow"
require_nonempty_file "$script_dir/scripts/refine_bundle_mappings.py" "refinement helper"

# Validate all numerical parameters before expensive alignment or output writes.
python3 - "$min_column_occupancy" "$min_base_frequency" "$merge_identity" \
    "$min_bundle_coverage" "$min_identity" "$max_overlap_fraction" \
    "$fixed_bundle_min_coverage" "$locus_dedup_min_bundle_coverage" \
    "$max_fixed_overlap_fraction" \
    "$min_mapq" "$max_bundle_gap" "$max_genome_gap" "$max_fixed_overlap_bp" \
    "$split_gap_min_mapping_bp" "$split_gap_max_overlap_bp" "$split_min_terminal_bp" \
    "$rank1_self_min_overlap_fraction" "$rank1_self_min_overlap_bp" \
    "$final_bed" "$pgr_bed" "$population_fasta" "$trf_bed" <<'PY'
import sys
from pathlib import Path

fractions = ("min-column-occupancy", "min-base-frequency", "merge-identity",
             "min-bundle-coverage", "min-identity", "max-overlap-fraction",
             "fixed-bundle-min-coverage", "locus-dedup-min-bundle-coverage",
             "max-fixed-overlap-fraction")
integers = ("min-mapq", "max-bundle-gap", "max-genome-gap", "max-fixed-overlap-bp",
            "mapping-split-gap-min-mapping-bp", "mapping-split-gap-max-overlap-bp",
            "mapping-split-min-terminal-bp")
try:
    for name, value in zip(fractions, sys.argv[1:10]):
        if not 0 <= float(value) <= 1:
            raise ValueError(f"--{name} must be in [0, 1]")
    for name, value in zip(integers, sys.argv[10:17]):
        if int(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    if int(sys.argv[15]) > int(sys.argv[14]):
        raise ValueError("--mapping-split-gap-max-overlap-bp must not exceed "
                         "--mapping-split-gap-min-mapping-bp")
    if not 0 <= float(sys.argv[17]) <= 1:
        raise ValueError("--rank1-self-min-overlap-fraction must be in [0, 1]")
    if int(sys.argv[18]) < 0:
        raise ValueError("--rank1-self-min-overlap-bp must be non-negative")
    if Path(sys.argv[19]).resolve() in {Path(p).resolve() for p in sys.argv[20:] if p}:
        raise ValueError("final output BED must not overwrite an input file")
except ValueError as error:
    sys.exit(f"error: {error}")
PY

python3 - "$trf_min_length" "$trf_min_identity" "$trf_full_coverage" \
    "$trf_bundle_removal_frequency" "$trf_split_min_bp" \
    "$trf_min_sequence_coverage" "$trf_min_population_presence" \
    "$trf_min_length_divergence" "$trf_min_nucleotide_divergence" \
    "$trf_min_fragment_sequence_bp" \
    "$trf_min_fragment_count" "$trf_min_fragment_frequency" \
    "$final_trf_full_coverage" <<'PY'
import sys
names = ("trf-min-length", "trf-min-identity", "trf-full-coverage",
         "trf-bundle-removal-frequency", "trf-split-min-bp", "trf-min-sequence-coverage",
         "trf-min-population-presence", "trf-min-length-divergence",
         "trf-min-nucleotide-divergence",
         "trf-min-fragment-sequence-bp", "trf-min-fragment-count", "trf-min-fragment-frequency",
         "final-trf-full-coverage")
try:
    if int(sys.argv[1]) < 0:
        raise ValueError("--trf-min-length must be non-negative")
    for index in (4, 9, 10):
        if int(sys.argv[index + 1]) < 1:
            raise ValueError(f"--{names[index]} must be positive")
    if not 0 <= float(sys.argv[2]) <= 100:
        raise ValueError("--trf-min-identity must be in [0, 100] (percent)")
    for index in (2, 3, 5, 6, 7, 8, 11, 12):
        if not 0 <= float(sys.argv[index + 1]) <= 1:
            raise ValueError(f"--{names[index]} must be in [0, 1]")
    for index in (2, 3, 12):
        if float(sys.argv[index + 1]) <= 0:
            raise ValueError(f"--{names[index]} must be in (0, 1]")
except ValueError as error:
    sys.exit(f"error: {error}")
PY

python3 - "$common_deletion_min_bp" "$common_deletion_min_count" \
    "$common_deletion_min_frequency" "$common_deletion_event_overlap" <<'PY'
import sys
try:
    if int(sys.argv[1]) < 0:
        raise ValueError("--common-deletion-min-bp must be non-negative")
    if int(sys.argv[2]) < 1:
        raise ValueError("--common-deletion-min-count must be positive")
    for name, value in zip(("common-deletion-min-frequency",
                            "common-deletion-event-overlap"), sys.argv[3:5]):
        if not 0 < float(value) <= 1:
            raise ValueError(f"--{name} must be in (0, 1]")
except ValueError as error:
    sys.exit(f"error: {error}")
PY

if [[ ! "$addition_min_bp" =~ ^[0-9]+$ ]] || [[ "$addition_min_bp" -le 0 ]]; then
    echo "error: --addition-min-bp must be a positive integer" >&2
    exit 2
fi

# A repeated FASTA ID is ambiguous as a population path. Exclude every copy
# and its linked annotations before ranking, remapping, or resume checks.
duplicate_filter_dir="$output_dir/00_input_qc"
duplicate_filter_args=(
    --fasta "$population_fasta" --bed "$pgr_bed"
    --output-dir "$duplicate_filter_dir"
)
if [[ -n "$trf_bed" ]]; then duplicate_filter_args+=(--trf-bed "$trf_bed"); fi
if [[ -n "$sample_map" ]]; then duplicate_filter_args+=(--sample-map "$sample_map"); fi
duplicate_id_count=$(python3 "$script_dir/scripts/filter_duplicate_fasta_ids.py" "${duplicate_filter_args[@]}")
if [[ "$duplicate_id_count" -gt 0 ]]; then
    population_fasta="$duplicate_filter_dir/population.unique_ids.fa"
    pgr_bed="$duplicate_filter_dir/pgrtk.unique_ids.bed"
    if [[ -n "$trf_bed" ]]; then trf_bed="$duplicate_filter_dir/trf.unique_ids.bed"; fi
    if [[ -n "$sample_map" ]]; then sample_map="$duplicate_filter_dir/sample_map.unique_ids.tsv"; fi
    echo "Excluded $duplicate_id_count duplicated FASTA IDs; see $duplicate_filter_dir/duplicate_fasta_ids.tsv"
fi

library_args=(
    --bed "$pgr_bed" --fasta "$population_fasta"
    --output-dir "$consensus_dir" --ranking-dir "$alignment_dir"
    --threads "$threads" --mafft "$mafft_cmd" --minimap2 "$minimap2_cmd"
    "--mafft-params=$mafft_02_params"
    --mafft-window-bp "$mafft_02_window_bp"
    --min-column-occupancy "$min_column_occupancy" --min-base-frequency "$min_base_frequency"
    --merge-identity "$merge_identity" --addition-min-bp "$addition_min_bp"
    --rank1-self-min-overlap-fraction "$rank1_self_min_overlap_fraction"
    --rank1-self-min-overlap-bp "$rank1_self_min_overlap_bp"
    --min-bundle-coverage "$min_bundle_coverage" --min-identity "$min_identity" --min-mapq "$min_mapq"
    --max-bundle-gap "$max_bundle_gap" --max-genome-gap "$max_genome_gap"
    --max-overlap-fraction "$max_overlap_fraction" --fixed-bundle-min-coverage "$fixed_bundle_min_coverage"
    --max-fixed-overlap-fraction "$max_fixed_overlap_fraction" --max-fixed-overlap-bp "$max_fixed_overlap_bp"
)
trf_args=(
    --trf-min-length "$trf_min_length" --trf-min-identity "$trf_min_identity"
    --trf-full-coverage "$trf_full_coverage"
    --trf-bundle-removal-frequency "$trf_bundle_removal_frequency"
    --trf-split-min-bp "$trf_split_min_bp"
    --trf-min-sequence-coverage "$trf_min_sequence_coverage"
    --trf-min-population-presence "$trf_min_population_presence"
    --trf-min-length-divergence "$trf_min_length_divergence"
    --trf-min-nucleotide-divergence "$trf_min_nucleotide_divergence"
    --trf-min-fragment-sequence-bp "$trf_min_fragment_sequence_bp"
    --trf-min-fragment-count "$trf_min_fragment_count"
    --trf-min-fragment-frequency "$trf_min_fragment_frequency"
)
if [[ -n "$trf_bed" ]]; then trf_args+=(--trf-bed "$trf_bed"); fi
library_args+=("${trf_args[@]}")
if [[ -n "$sample_map" ]]; then library_args+=(--sample-map "$sample_map"); fi
if [[ -n "$mafft_extra" ]]; then library_args+=("--mafft-extra=$mafft_extra"); fi
if [[ "$mafft_02_quick" -eq 1 ]]; then library_args+=(--mafft-quick); fi
if [[ "$resume" -eq 1 ]]; then library_args+=(--resume); fi
if [[ "$resume_ranks" -eq 1 ]]; then library_args+=(--resume-ranks); fi
library_args+=(${mapping_split_args[@]+"${mapping_split_args[@]}"})
library_args+=(${rank1_self_split_args[@]+"${rank1_self_split_args[@]}"})
library_args+=(${addition_self_split_args[@]+"${addition_self_split_args[@]}"})
library_args+=(${self_overlap_reconciliation_args[@]+"${self_overlap_reconciliation_args[@]}"})
library_args+=(${intermediate_args[@]+"${intermediate_args[@]}"})

echo "Starting at step $start_step; original bundle metadata: $pgr_bed"
echo "Initial PGR-TK refinement; no merged-BED coverage comparison"
if [[ "$start_step" -eq 1 ]]; then
    echo "[1/3] Grouping signed bundle orders and ranking haplotypes by sample count"
    python3 "$library_script" "${library_args[@]}" --prepare-only
else
    echo "[1/3] Skipped; reusing the haplotype ranking"
fi

if [[ "$start_step" -le 2 ]]; then
    echo "[2/3] Building the bundle library by iterative ranked remapping and additions"
    python3 "$library_script" "${library_args[@]}"
else
    echo "[2/3] Skipped; reusing the completed ranked bundle library"
fi

python3 - "$script_dir/scripts" "$consensus_dir/bundle.fa" "$consensus_dir/library_summary.json" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from build_ranked_bundle_library import validate_ranked_library
try:
    validate_ranked_library(Path(sys.argv[2]), Path(sys.argv[3]))
except (OSError, ValueError) as error:
    sys.exit(f"error: {error}")
PY

echo "[3/3] Refining raw PGR-TK bundles by remapping and mapping-driven splitting"
common_deletion_args=()
if [[ "$common_deletion_split" -eq 1 ]]; then
    common_deletion_args=(
        --common-deletion-split
        --common-deletion-min-bp "$common_deletion_min_bp"
        --common-deletion-min-count "$common_deletion_min_count"
        --common-deletion-min-frequency "$common_deletion_min_frequency"
        --common-deletion-event-overlap "$common_deletion_event_overlap"
        --common-deletion-merge-identity "$merge_identity"
        --mafft "$mafft_cmd"
        --mafft-params "$mafft_03_params"
    )
    if [[ -n "$sample_map" ]]; then
        common_deletion_args+=(--sample-map "$sample_map")
    fi
    if [[ -n "$mafft_extra" ]]; then
        common_deletion_args+=("--mafft-extra=$mafft_extra")
    fi
fi
post_filter_needle_args=()
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    post_filter_needle_args=(--post-filter-needle-refinement --needle "$needle_cmd")
fi
final_trf_args=()
if [[ -n "$trf_bed" ]] && awk '
    /^[[:space:]]*($|#)/ {next}
    {found = 1; exit}
    END {exit(found ? 0 : 1)}
' "$trf_bed"; then
    final_trf_args=(
        --final-trf-bed "$trf_bed"
        --final-trf-min-length "$trf_min_length"
        --final-trf-min-identity "$trf_min_identity"
        --final-trf-full-coverage "$final_trf_full_coverage"
        --final-trf-removal-frequency "$trf_bundle_removal_frequency"
    )
fi
bash "$refine_script" \
    --no-final-bundle-filter \
    --bundle-fasta "$consensus_dir/bundle.fa" \
    --genome-fasta "$population_fasta" \
    --bundle-metadata-bed "$consensus_dir/bundle_metadata.bed" \
    --output-dir "$refine_dir" \
    --minimap2 "$minimap2_cmd" --threads "$threads" \
    --min-bundle-coverage "$min_bundle_coverage" \
    --min-identity "$min_identity" --min-mapq "$min_mapq" \
    --max-bundle-gap "$max_bundle_gap" --max-genome-gap "$max_genome_gap" \
    --max-overlap-fraction "$max_overlap_fraction" \
    --fixed-bundle-min-coverage "$fixed_bundle_min_coverage" \
    --locus-dedup-min-bundle-coverage "$locus_dedup_min_bundle_coverage" \
    --max-fixed-overlap-fraction "$max_fixed_overlap_fraction" \
    --max-fixed-overlap-bp "$max_fixed_overlap_bp" \
    --final-self-min-overlap-fraction "$rank1_self_min_overlap_fraction" \
    --final-self-min-overlap-bp "$rank1_self_min_overlap_bp" \
    ${post_filter_needle_args[@]+"${post_filter_needle_args[@]}"} \
    ${final_trf_args[@]+"${final_trf_args[@]}"} \
    ${self_overlap_reconciliation_args[@]+"${self_overlap_reconciliation_args[@]}"} \
    ${common_deletion_args[@]+"${common_deletion_args[@]}"} \
    ${intermediate_args[@]+"${intermediate_args[@]}"} \
    ${mapping_split_args[@]+"${mapping_split_args[@]}"}

# Publish only after successful post-processing; original PGR-TK inputs are untouched.
cp "$refine_dir/refined_bundle_mappings.pgr.bed" "$final_bed"
echo "Pipeline complete"
echo "Per-round PAF/BED/truncation index: $consensus_dir/round_artifacts.tsv"
echo "New PGR-compatible BED (all retained mappings): $final_bed"
echo "Final split bundle consensus: $refine_dir/final_bundle.fa"
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    echo "Post-filter Needle audit: $refine_dir/post_filter_needle"
fi
if [[ "$common_deletion_split" -eq 1 ]]; then
    echo "Common-indel audit: $refine_dir/common_deletion_split/common_deletion_summary.json"
fi
