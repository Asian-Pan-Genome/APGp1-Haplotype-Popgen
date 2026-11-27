#!/usr/bin/env bash

set -euo pipefail

# End-to-end workflow:
# raw PGR-TK BED -> population graph cleanup/merge -> optional recursive TRF
# splitting of selected bundles -> sample-ranked iterative consensus library
# FASTA -> minimap2 remapping and previous-coverage comparison.

usage() {
    cat <<'EOF'
Usage:
  bash run_bundle_consensus_pipeline.sh \
    --pgr-bed INPUT.bed \
    --fasta POPULATION.fa \
    --trf-bed TRF.bed \
    --output-dir OUTPUT_DIR [options]

Required:
  --fasta FILE                FASTA whose record IDs match BED path IDs.
  --output-dir DIR            Pipeline output directory.

Required for step 1:
  --pgr-bed FILE              Raw four-column pgr-pbundle-decomp BED.
  --trf-bed FILE              TRF BED using the same paths and coordinates.
                              Also needed for requested TRF splits in step 2.
                              An empty/comment-only BED disables TRF steps.
                              Re-supply for ranking/library restarts at steps 3/4;
                              omitting it at final-only step 5 disables the
                              final TRF filter.

Optional:
  --start-step INT            Start at step 1--5 and run through step 5
                              (default: 1). Reuse earlier outputs in --output-dir.
  --sample-map FILE           Two-column path_id/sample_id map. Rank haplotypes
                              by distinct samples (default: one sample per path).
  --graph-coverage-denominator MODE
                              Step-1 short-bundle coverage denominator: auto,
                              fasta-length, or input-bed-union (default: auto).
                              Auto uses the supplied population FASTA lengths.
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
  --threads INT               Threads for MAFFT and minimap2 (default: 8).
  --mafft PATH                MAFFT executable (default: mafft).
  --minimap2 PATH             minimap2 executable for consensus remapping
                              (default: minimap2).
  --mafft-extra STRING        Extra MAFFT arguments as one quoted string.
  --split-bundle-id ID        Recursively TRF-split this numeric bundle ID.
                              May be repeated; the ID must occur at most once
                              per path. No bundle is MSA-split by default.
  --min-column-occupancy F    Consensus column occupancy (default: 0.50).
  --min-base-frequency F      Majority-base frequency or emit N (default: 0.50).
  --merge-identity F          Shared overlap identity >F: one representative;
                              otherwise retain same-class alleles (default: 0.90).
                              Also used for same-bundle and common-deletion children.
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
  --novel-coverage-fraction F Report a remapped bundle when previously uncovered
                              bp reach this bundle-length fraction (default: 0.25).
  --min-bundle-coverage F     Chained alignment bundle coverage (default: 0.60).
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
  --no-mapping-split          Disable mapping-driven splits in ranked rounds and
                              step 5; addition overlap clipping remains active.
  --common-deletion-min-bp BP Final-MSA internal deletion/insertion length must
                              be strictly greater than BP (default: 1000).
  --common-deletion-min-count INT
                              Minimum distinct samples (with --sample-map), else
                              paths, supporting one indel (default: 5).
  --common-deletion-min-frequency F
                              Minimum mapped population-unit frequency supporting
                              one indel (default: 0.01).
  --common-deletion-event-overlap F
                              Reciprocal overlap for clustering deletion events
                              (default: 0.80).
  --no-common-deletion-split  Skip final occurrence MSA/indel splitting.
  --resume                    Reuse a completed iterative library only when its
                              input/parameter fingerprint is unchanged. Also
                              reuse unchanged requested TRF-split alignments.
  --resume-ranks              Resume after the last completed rank checkpoint;
                              starts at step 4 and reruns any interrupted rank.
  --keep-intermediates        Keep disposable MAFFT/rank FASTAs and every self
                              snapshot; default cleans successful-stage temporaries.
  --no-merge-common-short-paths
                              Preserve common short alternative paths.
  -h, --help                  Show this message.

For the current COMPLEX_17 analysis, add: --split-bundle-id 89

Steps:
  1  Graph cleanup/merge (retains unmerged non-rare short bundles)
  2  Optional recursive TRF splitting (restart from the graph-merged BED)
  3  Optionally TRF-filter, group signed orders and rank haplotypes by sample count
  4  TRF-refine candidates; self-split rank 1 and each newly extended library
  5  minimap2 refinement, final occurrence MSA/common-indel splitting and
     final bundle self-alignment/remapping, final per-sequence TRF filtering,
     then coverage comparison and final BED class filtering
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_cmd="python3"
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
trf_enabled=0
trf_empty=0
output_dir=""
start_step=1
sample_map=""
graph_coverage_denominator="auto"
threads=8
mafft_cmd="mafft"
minimap2_cmd="minimap2"
mafft_extra=""
resume=0
resume_ranks=0
merge_common_short_paths=1
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
novel_coverage_fraction=0.25
min_bundle_coverage=0.60
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
split_bundle_ids=""
mapping_split_args=()
common_deletion_split=1
common_deletion_min_bp=1000
common_deletion_min_count=5
common_deletion_min_frequency=0.01
common_deletion_event_overlap=0.80

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pgr-bed)
            pgr_bed="${2:?missing value for --pgr-bed}"
            shift 2
            ;;
        --fasta)
            population_fasta="${2:?missing value for --fasta}"
            shift 2
            ;;
        --trf-bed)
            trf_bed="${2:?missing value for --trf-bed}"
            shift 2
            ;;
        --trf-min-length|--trf-min-identity|--trf-full-coverage|--trf-bundle-removal-frequency|--trf-split-min-bp|--trf-min-sequence-coverage|--trf-min-population-presence|--trf-min-length-divergence|--trf-min-nucleotide-divergence|--trf-min-fragment-sequence-bp|--trf-min-fragment-count|--trf-min-fragment-frequency)
            trf_option="${1#--}"
            trf_option="${trf_option//-/_}"
            printf -v "$trf_option" '%s' "${2:?missing TRF option value}"
            shift 2
            ;;
        --final-trf-full-coverage)
            final_trf_full_coverage="${2:?missing final TRF full-coverage value}"
            shift 2
            ;;
        --output-dir)
            output_dir="${2:?missing value for --output-dir}"
            shift 2
            ;;
        --start-step)
            start_step="${2:?missing value for --start-step}"
            shift 2
            ;;
        --sample-map)
            sample_map="${2:?missing value for --sample-map}"
            shift 2
            ;;
        --graph-coverage-denominator)
            graph_coverage_denominator="${2:?missing value for --graph-coverage-denominator}"
            shift 2
            ;;
        --threads)
            threads="${2:?missing value for --threads}"
            shift 2
            ;;
        --mafft)
            mafft_cmd="${2:?missing value for --mafft}"
            shift 2
            ;;
        --minimap2)
            minimap2_cmd="${2:?missing value for --minimap2}"
            shift 2
            ;;
        --mafft-extra)
            mafft_extra="${2:?missing value for --mafft-extra}"
            shift 2
            ;;
        --split-bundle-id)
            if [[ -z "$split_bundle_ids" ]]; then
                split_bundle_ids="${2:?missing value for --split-bundle-id}"
            else
                split_bundle_ids="$split_bundle_ids ${2:?missing value for --split-bundle-id}"
            fi
            shift 2
            ;;
        --min-column-occupancy)
            min_column_occupancy="${2:?missing value for --min-column-occupancy}"
            shift 2
            ;;
        --min-base-frequency)
            min_base_frequency="${2:?missing value for --min-base-frequency}"
            shift 2
            ;;
        --merge-identity)
            merge_identity="${2:?missing value for --merge-identity}"
            shift 2
            ;;
        --addition-min-bp)
            addition_min_bp="${2:?missing value for --addition-min-bp}"
            shift 2
            ;;
        --rank1-self-min-overlap-fraction)
            rank1_self_min_overlap_fraction="${2:?missing rank1 self-split value}"
            shift 2
            ;;
        --rank1-self-min-overlap-bp)
            rank1_self_min_overlap_bp="${2:?missing rank1 self-split value}"
            shift 2
            ;;
        --no-rank1-self-split)
            rank1_self_split_args+=("$1")
            shift
            ;;
        --no-addition-self-split)
            addition_self_split_args+=("$1")
            shift
            ;;
        --no-self-overlap-reconciliation)
            self_overlap_reconciliation_args+=("$1")
            shift
            ;;
        --keep-intermediates)
            intermediate_args+=("$1")
            shift
            ;;
        --novel-coverage-fraction)
            novel_coverage_fraction="${2:?missing value for --novel-coverage-fraction}"
            shift 2
            ;;
        --min-bundle-coverage)
            min_bundle_coverage="${2:?missing value for --min-bundle-coverage}"
            shift 2
            ;;
        --max-bundle-gap)
            max_bundle_gap="${2:?missing value for --max-bundle-gap}"
            shift 2
            ;;
        --max-genome-gap)
            max_genome_gap="${2:?missing value for --max-genome-gap}"
            shift 2
            ;;
        --max-overlap-fraction)
            max_overlap_fraction="${2:?missing value for --max-overlap-fraction}"
            shift 2
            ;;
        --fixed-bundle-min-coverage)
            fixed_bundle_min_coverage="${2:?missing value for --fixed-bundle-min-coverage}"
            shift 2
            ;;
        --locus-dedup-min-bundle-coverage)
            locus_dedup_min_bundle_coverage="${2:?missing value for --locus-dedup-min-bundle-coverage}"
            shift 2
            ;;
        --max-fixed-overlap-fraction)
            max_fixed_overlap_fraction="${2:?missing value for --max-fixed-overlap-fraction}"
            shift 2
            ;;
        --max-fixed-overlap-bp)
            max_fixed_overlap_bp="${2:?missing value for --max-fixed-overlap-bp}"
            shift 2
            ;;
        --mapping-split-gap-min-mapping-bp)
            split_gap_min_mapping_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_gap_min_mapping_bp")
            shift 2
            ;;
        --mapping-split-gap-max-overlap-bp)
            split_gap_max_overlap_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_gap_max_overlap_bp")
            shift 2
            ;;
        --mapping-split-min-terminal-bp)
            split_min_terminal_bp="${2:?missing mapping-split value}"
            mapping_split_args+=("$1" "$split_min_terminal_bp")
            shift 2
            ;;
        --mapping-split-min-truncated-count|--mapping-split-truncated-fraction|--mapping-split-max-unmapped-bp|--mapping-split-cut-merge-bp|--mapping-split-min-cut-support)
            mapping_split_args+=("$1" "${2:?missing mapping-split value}")
            shift 2
            ;;
        --no-mapping-split)
            mapping_split_args+=("$1")
            shift
            ;;
        --common-deletion-min-bp)
            common_deletion_min_bp="${2:?missing common-deletion value}"
            shift 2
            ;;
        --common-deletion-min-count)
            common_deletion_min_count="${2:?missing common-deletion value}"
            shift 2
            ;;
        --common-deletion-min-frequency)
            common_deletion_min_frequency="${2:?missing common-deletion value}"
            shift 2
            ;;
        --common-deletion-event-overlap)
            common_deletion_event_overlap="${2:?missing common-deletion value}"
            shift 2
            ;;
        --no-common-deletion-split)
            common_deletion_split=0
            shift
            ;;
        --resume)
            resume=1
            shift
            ;;
        --resume-ranks)
            resume_ranks=1
            shift
            ;;
        --no-merge-common-short-paths)
            merge_common_short_paths=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "error: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! "$start_step" =~ ^[1-5]$ ]]; then
    echo "error: --start-step must be an integer from 1 to 5" >&2
    exit 2
fi
if [[ "$resume_ranks" -eq 1 ]]; then
    if [[ "$start_step" -ne 1 && "$start_step" -ne 4 ]]; then
        echo "error: --resume-ranks requires the step-4 library build" >&2
        exit 2
    fi
    start_step=4
fi
case "$graph_coverage_denominator" in
    auto|fasta-length|input-bed-union) ;;
    *)
        echo "error: --graph-coverage-denominator must be auto, fasta-length, or input-bed-union" >&2
        exit 2
        ;;
esac
if [[ -z "$population_fasta" || -z "$output_dir" ]]; then
    echo "error: --fasta and --output-dir are required" >&2
    usage >&2
    exit 2
fi

if [[ ! "$threads" =~ ^[0-9]+$ ]] || [[ "$threads" -le 0 ]]; then
    echo "error: --threads must be a positive integer" >&2
    exit 2
fi

require_file() {
    local input_path="$1"
    local description="$2"
    if [[ -z "$input_path" || ! -f "$input_path" || ! -s "$input_path" ]]; then
        echo "error: $description is missing or empty: ${input_path:-<not supplied>}" >&2
        exit 2
    fi
}

merge_script="$script_dir/scripts/merge_pgr_bundles.py"
align_script="$script_dir/scripts/align_bundles_with_mafft.py"
split_script="$script_dir/scripts/split_bundle_by_trf_divergence.py"
library_script="$script_dir/scripts/build_ranked_bundle_library.py"
refine_script="$script_dir/scripts/refine_v2.sh"
graph_dir="$output_dir/01_graph_bundles"
split_dir="$output_dir/02_recursive_trf_split"
alignment_dir="$output_dir/03_bundle_alignments"
consensus_dir="$output_dir/04_bundle_consensus"
refine_dir="$output_dir/05_refined_bundle_mapping"
merge_prefix="$graph_dir/pgrtk_graph_merge_common"
final_bed="$split_dir/final_bundles.bed"

require_file "$population_fasta" "--fasta"
if [[ "$start_step" -eq 1 ]]; then
    require_file "$pgr_bed" "--pgr-bed required for step 1"
else
    require_file "$merge_prefix.bed" "previous graph-merged BED required from step $start_step"
fi
if [[ -n "$sample_map" ]]; then require_file "$sample_map" "--sample-map"; fi
if [[ "$start_step" -eq 1 || ( "$start_step" -eq 2 && -n "$split_bundle_ids" ) ]]; then
    if [[ -z "$trf_bed" ]]; then
        echo "error: --trf-bed required for step $start_step is missing: ${trf_bed:-<not supplied>}" >&2
        exit 2
    fi
fi
if [[ -n "$trf_bed" ]]; then
    if [[ ! -f "$trf_bed" ]]; then
        echo "error: --trf-bed file is missing: $trf_bed" >&2
        exit 2
    fi
    # Match the TRF readers: blank lines and comment lines are not annotations.
    if awk '
        /^[[:space:]]*($|#)/ {next}
        {found = 1; exit}
        END {exit(found ? 0 : 1)}
    ' "$trf_bed"; then
        trf_enabled=1
    else
        trf_data_status=$?
        if [[ "$trf_data_status" -eq 1 ]]; then
            trf_empty=1
        else
            echo "error: failed to read TRF BED: $trf_bed" >&2
            exit 2
        fi
    fi
fi
if [[ "$start_step" -eq 3 || "$start_step" -eq 4 ]]; then
    require_file "$final_bed" "final bundle BED required from step $start_step"
fi
if [[ "$start_step" -eq 4 ]]; then
    require_file "$alignment_dir/ranking_summary.json" "haplotype ranking required for step 4"
fi
if [[ "$start_step" -eq 5 ]]; then
    require_file "$consensus_dir/bundle.fa" "ranked bundle library required for step 5"
    require_file "$consensus_dir/library_summary.json" "completed library audit required for step 5; rerun steps 3/4 for an old library"
    require_file "$consensus_dir/bundle_metadata.bed" "ranked library metadata required for step 5"
fi

if ! command -v "$python_cmd" >/dev/null 2>&1; then
    echo "error: python3 was not found" >&2
    exit 2
fi
"$python_cmd" - "$max_bundle_gap" "$max_genome_gap" "$max_overlap_fraction" \
    "$fixed_bundle_min_coverage" "$locus_dedup_min_bundle_coverage" \
    "$max_fixed_overlap_fraction" \
    "$max_fixed_overlap_bp" "$min_bundle_coverage" \
    "$split_gap_min_mapping_bp" "$split_gap_max_overlap_bp" \
    "$split_min_terminal_bp" "$rank1_self_min_overlap_fraction" \
    "$rank1_self_min_overlap_bp" <<'PY'
import sys
try:
    for name, value in zip(("max-bundle-gap", "max-genome-gap"), sys.argv[1:3]):
        if int(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    for name, value in zip(("max-overlap-fraction", "fixed-bundle-min-coverage",
                            "locus-dedup-min-bundle-coverage",
                            "max-fixed-overlap-fraction"), sys.argv[3:7]):
        if not 0 <= float(value) <= 1:
            raise ValueError(f"--{name} must be in [0, 1]")
    if int(sys.argv[7]) < 0:
        raise ValueError("--max-fixed-overlap-bp must be non-negative")
    if not 0 <= float(sys.argv[8]) <= 1:
        raise ValueError("--min-bundle-coverage must be in [0, 1]")
    for name, value in zip(("mapping-split-gap-min-mapping-bp",
                            "mapping-split-gap-max-overlap-bp"), sys.argv[9:11]):
        if int(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    if int(sys.argv[10]) > int(sys.argv[9]):
        raise ValueError("--mapping-split-gap-max-overlap-bp must not exceed "
                         "--mapping-split-gap-min-mapping-bp")
    if int(sys.argv[11]) < 0:
        raise ValueError("--mapping-split-min-terminal-bp must be non-negative")
    if not 0 <= float(sys.argv[12]) <= 1:
        raise ValueError("--rank1-self-min-overlap-fraction must be in [0, 1]")
    if int(sys.argv[13]) < 0:
        raise ValueError("--rank1-self-min-overlap-bp must be non-negative")
except ValueError as error:
    sys.exit(f"error: {error}")
PY
"$python_cmd" - "$trf_min_length" "$trf_min_identity" "$trf_full_coverage" \
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
"$python_cmd" - "$common_deletion_min_bp" "$common_deletion_min_count" \
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
if [[ "$start_step" -le 4 || "$common_deletion_split" -eq 1 ]]; then
    if [[ "$mafft_cmd" == */* ]]; then
        if [[ ! -x "$mafft_cmd" ]]; then
            echo "error: MAFFT executable not found or not executable: $mafft_cmd" >&2
            exit 2
        fi
    elif ! command -v "$mafft_cmd" >/dev/null 2>&1; then
        echo "error: MAFFT was not found; install it or supply --mafft /path/to/mafft" >&2
        exit 2
    fi
fi
if [[ "$minimap2_cmd" == */* ]]; then
    if [[ ! -x "$minimap2_cmd" ]]; then
        echo "error: minimap2 executable not found or not executable: $minimap2_cmd" >&2
        exit 2
    fi
elif ! command -v "$minimap2_cmd" >/dev/null 2>&1; then
    echo "error: minimap2 was not found; supply --minimap2 /path/to/minimap2" >&2
    exit 2
fi
for bundle_id in $split_bundle_ids; do
    if [[ ! "$bundle_id" =~ ^[0-9]+$ ]]; then
        echo "error: --split-bundle-id must be numeric: $bundle_id" >&2
        exit 2
    fi
done

helpers=("$refine_script" "$library_script")
if [[ "$start_step" -eq 1 ]]; then
    helpers+=("$merge_script")
fi
if [[ "$start_step" -le 2 && "$trf_enabled" -eq 1 && -n "$split_bundle_ids" ]]; then
    helpers+=("$align_script")
fi
if [[ "$start_step" -le 2 && "$trf_enabled" -eq 1 && -n "$split_bundle_ids" ]]; then
    helpers+=("$split_script")
fi
for helper in "${helpers[@]}"; do
    if [[ ! -f "$helper" ]]; then
        echo "error: helper script not found: $helper" >&2
        exit 2
    fi
done

mkdir -p "$graph_dir" "$split_dir" "$alignment_dir" "$consensus_dir" "$refine_dir"

echo "Starting at step $start_step; earlier outputs are reused from: $output_dir"
if [[ "$trf_empty" -eq 1 ]]; then
    echo "TRF BED contains no annotation records; skipping all TRF filtering and splitting: $trf_bed"
fi
if [[ "$start_step" -le 1 ]]; then
    merge_args=(
        --fasta "$population_fasta"
        --coverage-denominator "$graph_coverage_denominator"
        --min-frequency 0.01
        --min-count 2
        --min-dup-frequency 0.01
        --min-dup-count 2
        --long-dup-break-bp 10000
        --long-dup-break-frequency 0.50
        --short-bp 10000
        --anchor-min-frequency 0.95
        --anchor-min-single-copy 0.95
        --min-path-frequency 0.02
        --min-path-count 2
        --max-trim-bp 10000
        --max-trim-bundles 5
        --subgraph-merge-min-shared-bp 10000
        --coordinate-overlap-seed 17
    )
    if [[ "$trf_enabled" -eq 1 ]]; then
        merge_args+=(
            --trf-bed "$trf_bed"
            --trf-min-length "$trf_min_length"
            --trf-min-identity "$trf_min_identity"
            --trf-full-coverage "$trf_full_coverage"
            --trf-bundle-removal-frequency "$trf_bundle_removal_frequency"
        )
    fi
    if [[ -n "$sample_map" ]]; then
        merge_args+=(--sample-map "$sample_map")
    fi
    if [[ "$merge_common_short_paths" -eq 1 ]]; then
        merge_args+=(--merge-common-short-paths)
    fi

    echo "[1/5] Cleaning and merging raw PGR-TK bundles"
    "$python_cmd" "$merge_script" "$pgr_bed" \
        --output-prefix "$merge_prefix" \
        "${merge_args[@]}"
else
    echo "[1/5] Skipped; reusing the graph-merged BED"
fi

if [[ "$start_step" -le 2 ]]; then
    # Step 2 always starts from the unsplit graph BED, not a previously split BED.
    cp "$merge_prefix.bed" "$final_bed"

    if [[ "$trf_empty" -eq 1 ]]; then
        echo "[2/5] Skipped; empty TRF BED, retaining the graph-cleaned BED"
    elif [[ -z "$split_bundle_ids" ]]; then
        echo "[2/5] Skipped; no --split-bundle-id supplied, retaining the graph-cleaned BED"
    else
        echo "[2/5] Applying requested MSA-guided recursive TRF splits"
        for bundle_id in $split_bundle_ids; do
            if ! awk -F '\t' -v target="$bundle_id" '
                $0 !~ /^#/ {split($4, value, ":"); if (value[1] == target) found = 1}
                END {exit(found ? 0 : 1)}
            ' "$final_bed"; then
                echo "warning: bundle $bundle_id is absent from $final_bed; skipping" >&2
                continue
            fi
            bundle_work_dir="$split_dir/bundle_$bundle_id"
            prealign_dir="$bundle_work_dir/prealignment"
            decision_dir="$bundle_work_dir/decision"
            align_args=(
                "$population_fasta"
                "$final_bed"
                --output-dir "$prealign_dir"
                --bundle-id "$bundle_id"
                --threads "$threads"
                --mafft "$mafft_cmd"
            )
            if [[ -n "$mafft_extra" ]]; then
                align_args+=(--mafft-extra "$mafft_extra")
            fi
            if [[ "$resume" -eq 1 ]]; then
                align_args+=(--resume)
            fi
            "$python_cmd" "$align_script" "${align_args[@]}"

            bundle_alignment="$prealign_dir/alignments/bundle_$bundle_id.aln.fa"
            "$python_cmd" "$split_script" \
                "$bundle_alignment" \
                "$trf_bed" \
                "$final_bed" \
                --output-dir "$decision_dir" \
                --trf-min-length "$trf_min_length" \
                --trf-min-identity "$trf_min_identity" \
                --split-repeat-min-bp "$trf_split_min_bp" \
                --min-sequence-coverage "$trf_min_sequence_coverage" \
                --min-population-presence "$trf_min_population_presence" \
                --min-length-divergence "$trf_min_length_divergence" \
                --min-nucleotide-divergence "$trf_min_nucleotide_divergence" \
                --min-fragment-sequence-bp "$trf_min_fragment_sequence_bp" \
                --min-fragment-count "$trf_min_fragment_count" \
                --min-fragment-frequency "$trf_min_fragment_frequency"

            split_bed="$decision_dir/bundle_$bundle_id.split.bed"
            if [[ -s "$split_bed" ]]; then
                next_bed="$split_dir/final_bundles.next.bed"
                cp "$split_bed" "$next_bed"
                mv "$next_bed" "$final_bed"
                echo "Applied recursive split for bundle $bundle_id"
            else
                echo "Bundle $bundle_id had no eligible long variable repeat"
            fi
        done
    fi
else
    echo "[2/5] Skipped; TRF splitting is not rerun"
fi

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

if [[ "$start_step" -le 3 ]]; then
    echo "[3/5] Grouping signed bundle orders and ranking haplotypes by sample count"
    ranking_args=(--bed "$final_bed" --fasta "$population_fasta"
                  --output-dir "$consensus_dir" --ranking-dir "$alignment_dir" --prepare-only)
    if [[ -n "$sample_map" ]]; then ranking_args+=(--sample-map "$sample_map"); fi
    ranking_args+=("${trf_args[@]}")
    ranking_args+=(--rank1-self-min-overlap-fraction "$rank1_self_min_overlap_fraction"
                   --rank1-self-min-overlap-bp "$rank1_self_min_overlap_bp")
    ranking_args+=(${rank1_self_split_args[@]+"${rank1_self_split_args[@]}"})
    ranking_args+=(${addition_self_split_args[@]+"${addition_self_split_args[@]}"})
    ranking_args+=(${self_overlap_reconciliation_args[@]+"${self_overlap_reconciliation_args[@]}"})
    ranking_args+=(${intermediate_args[@]+"${intermediate_args[@]}"})
    "$python_cmd" "$library_script" "${ranking_args[@]}"
else
    echo "[3/5] Skipped; reusing the haplotype ranking"
fi

if [[ "$start_step" -le 4 ]]; then
    echo "[4/5] Building the bundle library by iterative ranked remapping and additions"
    library_args=(
        --bed "$final_bed" --fasta "$population_fasta"
        --output-dir "$consensus_dir" --ranking-dir "$alignment_dir"
        --threads "$threads" --mafft "$mafft_cmd" --minimap2 "$minimap2_cmd"
        --min-column-occupancy "$min_column_occupancy"
        --min-base-frequency "$min_base_frequency"
        --merge-identity "$merge_identity"
        --addition-min-bp "$addition_min_bp"
        --rank1-self-min-overlap-fraction "$rank1_self_min_overlap_fraction"
        --rank1-self-min-overlap-bp "$rank1_self_min_overlap_bp"
        --min-bundle-coverage "$min_bundle_coverage"
        --max-bundle-gap "$max_bundle_gap" --max-genome-gap "$max_genome_gap"
        --max-overlap-fraction "$max_overlap_fraction"
        --fixed-bundle-min-coverage "$fixed_bundle_min_coverage"
        --max-fixed-overlap-fraction "$max_fixed_overlap_fraction"
        --max-fixed-overlap-bp "$max_fixed_overlap_bp"
    )
    if [[ -n "$sample_map" ]]; then library_args+=(--sample-map "$sample_map"); fi
    if [[ -n "$mafft_extra" ]]; then library_args+=("--mafft-extra=$mafft_extra"); fi
    if [[ "$resume" -eq 1 ]]; then library_args+=(--resume); fi
    if [[ "$resume_ranks" -eq 1 ]]; then library_args+=(--resume-ranks); fi
    library_args+=(${mapping_split_args[@]+"${mapping_split_args[@]}"})
    library_args+=(${rank1_self_split_args[@]+"${rank1_self_split_args[@]}"})
    library_args+=(${addition_self_split_args[@]+"${addition_self_split_args[@]}"})
    library_args+=(${self_overlap_reconciliation_args[@]+"${self_overlap_reconciliation_args[@]}"})
    library_args+=(${intermediate_args[@]+"${intermediate_args[@]}"})
    library_args+=("${trf_args[@]}")
    "$python_cmd" "$library_script" "${library_args[@]}"
else
    echo "[4/5] Skipped; reusing the completed ranked bundle library"
fi

"$python_cmd" - "$script_dir/scripts" "$consensus_dir/bundle.fa" \
    "$consensus_dir/library_summary.json" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from build_ranked_bundle_library import validate_ranked_library
try:
    validate_ranked_library(Path(sys.argv[2]), Path(sys.argv[3]))
except (OSError, ValueError) as error:
    sys.exit(f"error: {error}")
PY

echo "[5/5] Remapping the completed ranked bundle library and checking novel coverage"
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
    )
    if [[ -n "$sample_map" ]]; then
        common_deletion_args+=(--sample-map "$sample_map")
    fi
    if [[ -n "$mafft_extra" ]]; then
        common_deletion_args+=("--mafft-extra=$mafft_extra")
    fi
fi
final_trf_args=()
if [[ "$trf_enabled" -eq 1 ]]; then
    final_trf_args=(
        --final-trf-bed "$trf_bed"
        --final-trf-min-length "$trf_min_length"
        --final-trf-min-identity "$trf_min_identity"
        --final-trf-full-coverage "$final_trf_full_coverage"
        --final-trf-removal-frequency "$trf_bundle_removal_frequency"
    )
fi
bash "$refine_script" \
    --bundle-fasta "$consensus_dir/bundle.fa" \
    --filter-final-bundles \
    --final-max-bundle-bp 500 \
    --genome-fasta "$population_fasta" \
    --merged-bed "$merge_prefix.bed" \
    --bundle-metadata-bed "$consensus_dir/bundle_metadata.bed" \
    --output-dir "$refine_dir" \
    --minimap2 "$minimap2_cmd" \
    --threads "$threads" \
    --novel-fraction "$novel_coverage_fraction" \
    --min-bundle-coverage "$min_bundle_coverage" \
    --max-bundle-gap "$max_bundle_gap" --max-genome-gap "$max_genome_gap" \
    --max-overlap-fraction "$max_overlap_fraction" \
    --fixed-bundle-min-coverage "$fixed_bundle_min_coverage" \
    --locus-dedup-min-bundle-coverage "$locus_dedup_min_bundle_coverage" \
    --max-fixed-overlap-fraction "$max_fixed_overlap_fraction" \
    --max-fixed-overlap-bp "$max_fixed_overlap_bp" \
    --final-self-min-overlap-fraction "$rank1_self_min_overlap_fraction" \
    --final-self-min-overlap-bp "$rank1_self_min_overlap_bp" \
    ${self_overlap_reconciliation_args[@]+"${self_overlap_reconciliation_args[@]}"} \
    ${common_deletion_args[@]+"${common_deletion_args[@]}"} \
    ${final_trf_args[@]+"${final_trf_args[@]}"} \
    ${intermediate_args[@]+"${intermediate_args[@]}"} \
    ${mapping_split_args[@]+"${mapping_split_args[@]}"}

echo "Pipeline complete"
if [[ "$start_step" -le 4 ]]; then
    echo "Final bundle BED: $final_bed"
    echo "Combined consensus FASTA: $consensus_dir/bundle_consensus.fa"
    echo "Ranked library summary: $consensus_dir/library_summary.json"
fi
echo "Completed ranked bundle library: $consensus_dir/bundle.fa"
echo "Per-round PAF/BED/truncation index: $consensus_dir/round_artifacts.tsv"
echo "Refined bundle mappings: $refine_dir/refined_bundle_mappings.bed"
echo "PGR-compatible refined mappings: $refine_dir/refined_bundle_mappings.pgr.bed"
echo "Final split bundle consensus: $refine_dir/final_bundle.fa"
if [[ "$common_deletion_split" -eq 1 ]]; then
    echo "Common-indel audit: $refine_dir/common_deletion_split/common_deletion_summary.json"
fi
echo "Final BED filtering audit: $refine_dir/final_bundle_filter.tsv"
echo "Signed haplotypes before/after: $refine_dir/final_haplotypes.before.tsv / final_haplotypes.after.tsv"
echo "Novel-coverage report: $refine_dir/novel_bundle_coverage.tsv"
