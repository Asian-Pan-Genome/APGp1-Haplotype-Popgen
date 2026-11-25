#!/usr/bin/env bash

set -euo pipefail

# Remap rank-priority bundle consensuses to population genomes. Optionally compare
# coverage with a graph-merged PGR bundle BED when --merged-bed is supplied.

usage() {
    cat <<'EOF'
Usage:
  bash scripts/refine_v2.sh \
    --bundle-fasta bundle.fa \
    --genome-fasta COMPLEX_17.rename.fa \
    --bundle-metadata-bed original_pgrtk.bed \
    --output-dir refine_results [--merged-bed PREVIOUS.bed] [options]

Options:
  --bundle-fasta FILE         Bundle consensus FASTA (default: bundle.fa).
  --genome-fasta FILE         Population genome/haplotype FASTA
                              (default: COMPLEX_17.rename.fa).
  --merged-bed FILE           Optional previous merged BED for coverage comparison.
                              Omit to refine without novelty assessment/reports.
  --bundle-metadata-bed FILE  PGR BED supplying bundle IDs/sizes independently
                              of coverage (default: --merged-bed); required
                              when --merged-bed is omitted.
  --output-dir DIR            Output directory (default: refine_v2_results).
  --minimap2 PATH             minimap2 executable (default: minimap2).
  --threads INT               minimap2 threads (default: 8).
  --min-bundle-coverage F     Minimum aligned bundle fraction (default: 0.60).
  --min-identity F            Minimum PAF alignment identity (default: 0.0).
  --min-mapq INT              Minimum PAF MAPQ (default: 0).
  --max-bundle-gap BP         Maximum bundle-coordinate gap/overlap for joining
                              collinear PAF blocks (default: 1000).
  --max-genome-gap BP         Maximum population-genome gap/overlap for joining
                              blocks (default: 1000).
  --max-overlap-fraction F   Maximum overlap on either axis as a fraction of
                              the bundle consensus length (default: 0.10).
  --fixed-bundle-min-coverage F
                              Single-PAF coverage to fix a block and exclude
                              it from chaining (default: 0.80).
  --locus-dedup-min-bundle-coverage F
                              Minimum completed-chain union bundle coverage
                              admitted to final remap/locus dedup (default: 0.70).
  --max-fixed-overlap-fraction F
                              Maximum candidate-span overlap / fixed bundle
                              consensus length (default: 0.20).
  --max-fixed-overlap-bp BP   Maximum candidate-span overlap with a different
                              fixed class (default: 10000).
  --novel-fraction F          Report novel coverage when uncovered bp is at
                              least this fraction of bundle length (default: 0.25).
  --mapping-split-min-truncated-count INT
                              Gap-supported short-occurrence count (default: 5).
  --mapping-split-truncated-fraction F
                              Strict truncation length fraction (default: 0.60).
  --mapping-split-max-unmapped-bp BP
                              Report longer globally unaligned bundle intervals
                              (strictly greater; default: 10000).
  --mapping-split-cut-merge-bp BP
                              Cluster nearby candidate cuts (default: 100).
  --mapping-split-min-cut-support INT
                              Short-occurrence support per cut (default: 1;
                              still requires the bundle truncation count gate).
  --mapping-split-min-terminal-bp BP
                              Suppress cuts producing shorter terminal fragments
                              (default: 1000; 0 disables; internal pieces retained).
  --mapping-split-gap-min-mapping-bp BP
                              Require each selected PAF's mapped span inside
                              a preliminary-BED gap to exceed BP (default: 10000).
  --mapping-split-gap-max-overlap-bp BP
                              Require every selected pair's genomic overlap
                              to be less than BP (default: 1000; 0: no overlap).
                              Must not exceed --mapping-split-gap-min-mapping-bp.
  --no-mapping-split          Disable mapping-driven splitting/remapping.
  --common-deletion-split     After preliminary class-local locus deduplication,
                              align each
                              final bundle's retained occurrences with MAFFT and
                              split recurrent internal deletions or insertions. Off by default;
                              both top-level pipelines enable it.
  --common-deletion-min-bp BP Internal deletion/insertion length must be strictly
                              greater than BP (default: 1000).
  --common-deletion-min-count INT
                              Minimum distinct population paths supporting one
                              deletion/insertion (default: 5).
  --common-deletion-min-frequency F
                              Minimum fraction of mapped population paths
                              supporting one deletion/insertion (default: 0.01).
  --common-deletion-event-overlap F
                              Minimum overlap fraction of both event lengths used
                              to cluster shifted breakpoints (default: 0.80).
  --mafft PATH                MAFFT executable for common-indel MSA
                              (default: mafft).
  --mafft-extra STRING        Additional quoted MAFFT arguments.
  --mafft-params STRING       MAFFT parameters replacing --auto for the
                              common-indel MSA (default: --auto).
  --sample-map FILE           Optional path_id/sample_id map; common-indel
                              support is counted by sample instead of path.
  --common-deletion-merge-identity F
                              Child alleles with aligned identity >F collapse to
                              one sequence while retaining one class (default: 0.90).
  --final-self-min-overlap-fraction F
                              After the last successful split, iteratively
                              reconcile bundle self-overlaps covering more than
                              this fraction of either bundle (default: 0.20).
  --final-self-min-overlap-bp BP
                              Also reconcile final self-overlaps longer than BP
                              (default: 10000). Self minimap2 uses -g 500.
  --no-self-overlap-reconciliation
                              After the final split, use legacy self cuts instead
                              of overlap-family reconciliation.
  --post-filter-needle-refinement
                              After final mapping/coverage filtering, shortlist
                              genomically overlapping bundle pairs, confirm the
                              sequence overlap with EMBOSS Needle, reconcile one
                              pair, remap, and repeat to a fixed point.
  --needle PATH               EMBOSS needle executable (default: needle).
  --needle-gapopen F          Needle gap-open penalty (default: 10.0).
  --needle-gapextend F        Needle gap-extension penalty (default: 0.5).
  --post-filter-needle-max-rounds INT
                              Safety limit for remap/filter/Needle rounds
                              (default: 100).
  --final-trf-bed FILE        After Needle convergence, remove final numeric
                              bundle classes recurrently overlapping TRF and
                              remap the filtered full library. Disabled when
                              omitted; empty/comment-only input is a no-op.
  --final-trf-min-length BP   Minimum final-filter TRF interval length
                              (default: 0; no length filter).
  --final-trf-min-identity PCT
                              Minimum final-filter TRF identity (default: 80).
  --final-trf-full-coverage F
                              A single merged TRF interval must cover at least
                              this fraction of one mapped occurrence footprint
                              before its sequence qualifies (default: 0.95).
  --final-trf-removal-frequency F
                              Remove a class when at least this fraction of
                              its carrying sequences have a qualifying TRF-
                              covered occurrence (default: 0.50). Copies within
                              one sequence count once.
  --filter-final-bundles      After locus dedup, prune classes shortest first
                              with 95% sequence coverage / 5% population gate
                              and presence/absence-only haplotype merges; retain
                              candidates that do not change haplotype grouping.
                              Off by default; enabled by the complete pipeline.
  --no-final-bundle-filter    Disable final class pruning (initial remap default).
  --final-min-sequence-coverage F
                              Final sequence coverage threshold (default: 0.95).
  --final-max-low-coverage-fraction F
                              Allowed low-coverage sequence fraction (default: 0.05).
  --final-max-bundle-bp BP    Only bundle classes with median consensus length
                              <=BP are deletion candidates (default: 500).
  --keep-intermediates        Retain redundant mapping-split FASTA copies
                              (default: remove newly created duplicates only
                              after successful final post-processing).
  -h, --help                  Show this help.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bundle_fasta="bundle.fa"
genome_fasta="COMPLEX_17.rename.fa"
merged_bed=""
bundle_metadata_bed=""
output_dir="refine_v2_results"
minimap2_cmd="minimap2"
needle_cmd="needle"
mafft_cmd="mafft"
mafft_extra=""
mafft_params="--auto"
sample_map=""
threads=8
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
novel_fraction=0.25
mapping_split=1
split_min_truncated_count=5
split_truncated_fraction=0.60
split_max_unmapped_bp=10000
split_cut_merge_bp=100
split_min_cut_support=1
split_min_terminal_bp=1000
split_gap_min_mapping_bp=10000
split_gap_max_overlap_bp=1000
common_deletion_split=0
common_deletion_min_bp=1000
common_deletion_min_count=5
common_deletion_min_frequency=0.01
common_deletion_event_overlap=0.80
common_deletion_merge_identity=0.90
final_self_min_overlap_fraction=0.20
final_self_min_overlap_bp=10000
self_overlap_reconciliation=1
post_filter_needle_refinement=0
needle_gapopen=10.0
needle_gapextend=0.5
post_filter_needle_max_rounds=100
final_trf_bed=""
final_trf_min_length=0
final_trf_min_identity=80.0
final_trf_full_coverage=0.95
final_trf_removal_frequency=0.50
keep_intermediates=0
filter_final_bundles=0
final_min_sequence_coverage=0.95
final_max_low_coverage_fraction=0.05
final_max_bundle_bp=500

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bundle-fasta)
            bundle_fasta="${2:?missing value for --bundle-fasta}"
            shift 2
            ;;
        --genome-fasta)
            genome_fasta="${2:?missing value for --genome-fasta}"
            shift 2
            ;;
        --merged-bed)
            merged_bed="${2:?missing value for --merged-bed}"
            shift 2
            ;;
        --bundle-metadata-bed)
            bundle_metadata_bed="${2:?missing value for --bundle-metadata-bed}"
            shift 2
            ;;
        --output-dir)
            output_dir="${2:?missing value for --output-dir}"
            shift 2
            ;;
        --minimap2)
            minimap2_cmd="${2:?missing value for --minimap2}"
            shift 2
            ;;
        --needle)
            needle_cmd="${2:?missing value for --needle}"
            shift 2
            ;;
        --mafft)
            mafft_cmd="${2:?missing value for --mafft}"
            shift 2
            ;;
        --mafft-extra)
            mafft_extra="${2:?missing value for --mafft-extra}"
            shift 2
            ;;
        --mafft-params)
            mafft_params="${2:?missing value for --mafft-params}"
            shift 2
            ;;
        --sample-map)
            sample_map="${2:?missing value for --sample-map}"
            shift 2
            ;;
        --threads)
            threads="${2:?missing value for --threads}"
            shift 2
            ;;
        --min-bundle-coverage|--min-query-coverage)
            min_bundle_coverage="${2:?missing value for --min-bundle-coverage}"
            shift 2
            ;;
        --min-identity)
            min_identity="${2:?missing value for --min-identity}"
            shift 2
            ;;
        --min-mapq)
            min_mapq="${2:?missing value for --min-mapq}"
            shift 2
            ;;
        --max-bundle-gap|--max-target-gap)
            max_bundle_gap="${2:?missing value for --max-bundle-gap}"
            shift 2
            ;;
        --max-genome-gap|--max-query-gap)
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
        --novel-fraction)
            novel_fraction="${2:?missing value for --novel-fraction}"
            shift 2
            ;;
        --mapping-split-min-truncated-count)
            split_min_truncated_count="${2:?missing value}"; shift 2 ;;
        --mapping-split-truncated-fraction)
            split_truncated_fraction="${2:?missing value}"; shift 2 ;;
        --mapping-split-max-unmapped-bp)
            split_max_unmapped_bp="${2:?missing value}"; shift 2 ;;
        --mapping-split-cut-merge-bp)
            split_cut_merge_bp="${2:?missing value}"; shift 2 ;;
        --mapping-split-min-cut-support)
            split_min_cut_support="${2:?missing value}"; shift 2 ;;
        --mapping-split-min-terminal-bp)
            split_min_terminal_bp="${2:?missing value}"; shift 2 ;;
        --mapping-split-gap-min-mapping-bp)
            split_gap_min_mapping_bp="${2:?missing value}"; shift 2 ;;
        --mapping-split-gap-max-overlap-bp)
            split_gap_max_overlap_bp="${2:?missing value}"; shift 2 ;;
        --no-mapping-split)
            mapping_split=0; shift ;;
        --common-deletion-split)
            common_deletion_split=1; shift ;;
        --common-deletion-min-bp)
            common_deletion_min_bp="${2:?missing value}"; shift 2 ;;
        --common-deletion-min-count)
            common_deletion_min_count="${2:?missing value}"; shift 2 ;;
        --common-deletion-min-frequency)
            common_deletion_min_frequency="${2:?missing value}"; shift 2 ;;
        --common-deletion-event-overlap)
            common_deletion_event_overlap="${2:?missing value}"; shift 2 ;;
        --common-deletion-merge-identity)
            common_deletion_merge_identity="${2:?missing value}"; shift 2 ;;
        --final-self-min-overlap-fraction)
            final_self_min_overlap_fraction="${2:?missing value}"; shift 2 ;;
        --final-self-min-overlap-bp)
            final_self_min_overlap_bp="${2:?missing value}"; shift 2 ;;
        --no-self-overlap-reconciliation)
            self_overlap_reconciliation=0; shift ;;
        --post-filter-needle-refinement)
            post_filter_needle_refinement=1; shift ;;
        --no-post-filter-needle-refinement)
            post_filter_needle_refinement=0; shift ;;
        --needle-gapopen)
            needle_gapopen="${2:?missing value for --needle-gapopen}"; shift 2 ;;
        --needle-gapextend)
            needle_gapextend="${2:?missing value for --needle-gapextend}"; shift 2 ;;
        --post-filter-needle-max-rounds)
            post_filter_needle_max_rounds="${2:?missing value for --post-filter-needle-max-rounds}"; shift 2 ;;
        --final-trf-bed)
            final_trf_bed="${2:?missing value for --final-trf-bed}"; shift 2 ;;
        --final-trf-min-length)
            final_trf_min_length="${2:?missing value for --final-trf-min-length}"; shift 2 ;;
        --final-trf-min-identity)
            final_trf_min_identity="${2:?missing value for --final-trf-min-identity}"; shift 2 ;;
        --final-trf-full-coverage)
            final_trf_full_coverage="${2:?missing value for --final-trf-full-coverage}"; shift 2 ;;
        --final-trf-removal-frequency)
            final_trf_removal_frequency="${2:?missing value for --final-trf-removal-frequency}"; shift 2 ;;
        --filter-final-bundles)
            filter_final_bundles=1; shift ;;
        --no-final-bundle-filter)
            filter_final_bundles=0; shift ;;
        --final-min-sequence-coverage)
            final_min_sequence_coverage="${2:?missing value}"; shift 2 ;;
        --final-max-low-coverage-fraction)
            final_max_low_coverage_fraction="${2:?missing value}"; shift 2 ;;
        --final-max-bundle-bp)
            final_max_bundle_bp="${2:?missing value}"; shift 2 ;;
        --keep-intermediates)
            keep_intermediates=1; shift ;;
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

if [[ -z "$merged_bed" && -z "$bundle_metadata_bed" ]]; then
    echo "error: --bundle-metadata-bed is required without --merged-bed" >&2
    usage >&2
    exit 2
fi
if [[ -z "$bundle_metadata_bed" ]]; then
    bundle_metadata_bed="$merged_bed"
fi
input_paths=("$bundle_fasta" "$genome_fasta" "$bundle_metadata_bed")
if [[ -n "$merged_bed" ]]; then input_paths+=("$merged_bed"); fi
if [[ -n "$sample_map" ]]; then input_paths+=("$sample_map"); fi
if [[ -n "$final_trf_bed" ]]; then input_paths+=("$final_trf_bed"); fi
for input_path in "${input_paths[@]}"; do
    if [[ ! -f "$input_path" ]]; then
        echo "error: input file not found: $input_path" >&2
        exit 2
    fi
done
if [[ ! "$threads" =~ ^[0-9]+$ ]] || [[ "$threads" -le 0 ]]; then
    echo "error: --threads must be a positive integer" >&2
    exit 2
fi
if [[ "$minimap2_cmd" == */* ]]; then
    if [[ ! -x "$minimap2_cmd" ]]; then
        echo "error: minimap2 is not executable: $minimap2_cmd" >&2
        exit 2
    fi
elif ! command -v "$minimap2_cmd" >/dev/null 2>&1; then
    echo "error: minimap2 was not found; provide --minimap2 /path/to/minimap2" >&2
    exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "error: python3 was not found" >&2
    exit 2
fi
if [[ "$common_deletion_split" -eq 1 ]]; then
    if [[ "$mafft_cmd" == */* ]]; then
        if [[ ! -x "$mafft_cmd" ]]; then
            echo "error: MAFFT is not executable: $mafft_cmd" >&2
            exit 2
        fi
    elif ! command -v "$mafft_cmd" >/dev/null 2>&1; then
        echo "error: MAFFT was not found; provide --mafft /path/to/mafft" >&2
        exit 2
    fi
fi
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    if [[ "$needle_cmd" == */* ]]; then
        if [[ ! -x "$needle_cmd" ]]; then
            echo "error: EMBOSS needle is not executable: $needle_cmd" >&2
            exit 2
        fi
    elif ! command -v "$needle_cmd" >/dev/null 2>&1; then
        echo "error: EMBOSS needle was not found; provide --needle /path/to/needle" >&2
        exit 2
    fi
fi

postprocess_script="$script_dir/refine_bundle_mappings.py"
split_script="$script_dir/split_bundles_by_mapping.py"
cleanup_script="$script_dir/pipeline_intermediates.py"
common_deletion_script="$script_dir/split_bundles_by_common_deletions.py"
final_self_script="$script_dir/refine_self_aligned_library.py"
needle_refine_script="$script_dir/refine_overlapping_bundles_with_needle.py"
final_trf_filter_script="$script_dir/filter_final_library_by_trf.py"
if [[ ! -f "$postprocess_script" ]]; then
    echo "error: PAF post-processing script not found: $postprocess_script" >&2
    exit 2
fi
if [[ "$mapping_split" -eq 1 && ! -f "$split_script" ]]; then
    echo "error: mapping-split script not found: $split_script" >&2
    exit 2
fi
if [[ ! -f "$cleanup_script" ]]; then
    echo "error: intermediate-cleanup script not found: $cleanup_script" >&2
    exit 2
fi
if [[ "$common_deletion_split" -eq 1 && ! -f "$common_deletion_script" ]]; then
    echo "error: common-deletion splitter not found: $common_deletion_script" >&2
    exit 2
fi
if [[ ! -f "$final_self_script" ]]; then
    echo "error: final self-alignment script not found: $final_self_script" >&2
    exit 2
fi
if [[ "$post_filter_needle_refinement" -eq 1 && ! -f "$needle_refine_script" ]]; then
    echo "error: post-filter Needle refinement script not found: $needle_refine_script" >&2
    exit 2
fi
if [[ -n "$final_trf_bed" && ! -f "$final_trf_filter_script" ]]; then
    echo "error: final-library TRF filter not found: $final_trf_filter_script" >&2
    exit 2
fi
python3 - "$split_min_truncated_count" "$split_truncated_fraction" \
    "$split_max_unmapped_bp" "$split_cut_merge_bp" "$split_min_cut_support" \
    "$max_bundle_gap" "$max_genome_gap" "$max_overlap_fraction" \
    "$fixed_bundle_min_coverage" "$locus_dedup_min_bundle_coverage" \
    "$max_fixed_overlap_fraction" \
    "$max_fixed_overlap_bp" "$min_bundle_coverage" \
    "$split_gap_min_mapping_bp" "$split_gap_max_overlap_bp" \
    "$split_min_terminal_bp" "$final_min_sequence_coverage" \
    "$final_max_low_coverage_fraction" "$final_max_bundle_bp" \
    "$common_deletion_min_bp" "$common_deletion_min_count" \
    "$common_deletion_min_frequency" "$common_deletion_event_overlap" \
    "$common_deletion_merge_identity" "$final_self_min_overlap_fraction" \
    "$final_self_min_overlap_bp" "$needle_gapopen" "$needle_gapextend" \
    "$post_filter_needle_max_rounds" "$final_trf_min_length" \
    "$final_trf_min_identity" "$final_trf_full_coverage" \
    "$final_trf_removal_frequency" <<'PY'
import math
import sys
try:
    if int(sys.argv[1]) < 1 or not 0 <= float(sys.argv[2]) <= 1:
        raise ValueError("invalid mapping-split count or fraction")
    if int(sys.argv[3]) < 0 or int(sys.argv[4]) < 0:
        raise ValueError("mapping-split bp thresholds must be non-negative")
    if int(sys.argv[5]) < 1:
        raise ValueError("mapping-split cut support must be positive")
    for name, value in zip(("max-bundle-gap", "max-genome-gap"), sys.argv[6:8]):
        if int(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    for name, value in zip(("max-overlap-fraction", "fixed-bundle-min-coverage",
                            "locus-dedup-min-bundle-coverage",
                            "max-fixed-overlap-fraction"), sys.argv[8:12]):
        if not 0 <= float(value) <= 1:
            raise ValueError(f"--{name} must be in [0, 1]")
    if int(sys.argv[12]) < 0:
        raise ValueError("--max-fixed-overlap-bp must be non-negative")
    if not 0 <= float(sys.argv[13]) <= 1:
        raise ValueError("--min-bundle-coverage must be in [0, 1]")
    for name, value in zip(("mapping-split-gap-min-mapping-bp",
                            "mapping-split-gap-max-overlap-bp"), sys.argv[14:16]):
        if int(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    if int(sys.argv[15]) > int(sys.argv[14]):
        raise ValueError("--mapping-split-gap-max-overlap-bp must not exceed "
                         "--mapping-split-gap-min-mapping-bp")
    if int(sys.argv[16]) < 0:
        raise ValueError("--mapping-split-min-terminal-bp must be non-negative")
    for name, value in zip(("final-min-sequence-coverage", "final-max-low-coverage-fraction"), sys.argv[17:19]):
        if not 0 <= float(value) <= 1:
            raise ValueError(f"--{name} must be in [0, 1]")
    if int(sys.argv[19]) < 0:
        raise ValueError("--final-max-bundle-bp must be non-negative")
    if int(sys.argv[20]) < 0:
        raise ValueError("--common-deletion-min-bp must be non-negative")
    if int(sys.argv[21]) < 1:
        raise ValueError("--common-deletion-min-count must be positive")
    for name, value in zip(("common-deletion-min-frequency",
                            "common-deletion-event-overlap"), sys.argv[22:24]):
        if not 0 < float(value) <= 1:
            raise ValueError(f"--{name} must be in (0, 1]")
    if not 0 <= float(sys.argv[24]) <= 1:
        raise ValueError("--common-deletion-merge-identity must be in [0, 1]")
    if not 0 <= float(sys.argv[25]) <= 1:
        raise ValueError("--final-self-min-overlap-fraction must be in [0, 1]")
    if int(sys.argv[26]) < 0:
        raise ValueError("--final-self-min-overlap-bp must be non-negative")
    for name, value in (("needle-gapopen", sys.argv[27]),
                        ("needle-gapextend", sys.argv[28])):
        if not math.isfinite(float(value)) or float(value) < 0:
            raise ValueError(f"--{name} must be non-negative")
    if int(sys.argv[29]) < 1:
        raise ValueError("--post-filter-needle-max-rounds must be positive")
    if int(sys.argv[30]) < 0:
        raise ValueError("--final-trf-min-length must be non-negative")
    if not math.isfinite(float(sys.argv[31])) or not 0 <= float(sys.argv[31]) <= 100:
        raise ValueError("--final-trf-min-identity must be in [0, 100]")
    if not math.isfinite(float(sys.argv[32])) or not 0 < float(sys.argv[32]) <= 1:
        raise ValueError("--final-trf-full-coverage must be in (0, 1]")
    if not math.isfinite(float(sys.argv[33])) or not 0 < float(sys.argv[33]) <= 1:
        raise ValueError("--final-trf-removal-frequency must be in (0, 1]")
except ValueError as error:
    sys.exit(f"error: {error}")
PY

# Ownership is recorded before this invocation creates or overwrites anything.
# Existing files, including broken symlinks, must never become cleanup targets.
split_dir="$output_dir/mapping_split"
split_fasta_generated=0
split_summary_generated=0
if [[ "$mapping_split" -eq 1 && ! -e "$split_dir/bundle.fa" && ! -L "$split_dir/bundle.fa" ]]; then
    split_fasta_generated=1
fi
if [[ "$mapping_split" -eq 1 && ! -e "$split_dir/mapping_split_summary.json" && ! -L "$split_dir/mapping_split_summary.json" ]]; then
    split_summary_generated=1
fi
mkdir -p "$output_dir"
raw_paf="$output_dir/bundle_to_genome.raw.paf"
final_bundle_fasta="$output_dir/final_bundle.fa"
effective_paf="$raw_paf"
effective_metadata_bed="$bundle_metadata_bed"

map_consensuses() {
    local mapping_gap="$3"
    "$minimap2_cmd" -x asm20 -p 0.2 -N 10 -g "$mapping_gap" -c -t "$threads" \
        "$1" "$genome_fasta" > "$2"
}

if [[ "$common_deletion_split" -eq 1 ]]; then refine_stage_total=5; else refine_stage_total=4; fi
echo "[refine 1/$refine_stage_total] Mapping population haplotypes to rank-priority bundle consensuses"
echo "Bundle FASTA headers remain unique; aliases identify classes, never shared chaining coordinates."
map_consensuses "$bundle_fasta" "$raw_paf" 10000

split_count=0
if [[ "$mapping_split" -eq 1 ]]; then
    echo "[refine 2/$refine_stage_total] Selecting recurrent truncation evidence in preliminary-BED uncovered regions"
    python3 "$split_script" "$bundle_fasta" "$raw_paf" \
        --bundle-metadata-bed "$bundle_metadata_bed" --output-dir "$split_dir" \
        --min-truncated-count "$split_min_truncated_count" \
        --truncated-fraction "$split_truncated_fraction" \
        --max-unmapped-bp "$split_max_unmapped_bp" --cut-merge-bp "$split_cut_merge_bp" \
        --min-cut-support "$split_min_cut_support" \
        --min-terminal-fragment-bp "$split_min_terminal_bp" \
        --preliminary-min-bundle-coverage "$min_bundle_coverage" \
        --gap-min-mapping-bp "$split_gap_min_mapping_bp" \
        --gap-max-overlap-bp "$split_gap_max_overlap_bp" \
        --min-identity "$min_identity" --min-mapq "$min_mapq" \
        --max-bundle-gap "$max_bundle_gap" --max-genome-gap "$max_genome_gap" \
        --max-overlap-fraction "$max_overlap_fraction" \
        --fixed-bundle-min-coverage "$fixed_bundle_min_coverage" \
        --max-fixed-overlap-fraction "$max_fixed_overlap_fraction" \
        --max-fixed-overlap-bp "$max_fixed_overlap_bp"
    cp "$split_dir/bundle.fa" "$final_bundle_fasta"
    split_count="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["split_consensus_count"])' \
        "$split_dir/mapping_split_summary.json")"
    if [[ "$split_count" -gt 0 ]]; then
        echo "Split $split_count consensus bundle(s); remapping their globally numbered fragments"
        effective_paf="$output_dir/bundle_to_genome.split.raw.paf"
        effective_metadata_bed="$split_dir/bundle_metadata.bed"
        map_consensuses "$final_bundle_fasta" "$effective_paf" 500
    else
        echo "No eligible mapping-driven cuts; reusing the initial PAF"
    fi
else
    echo "[refine 2/$refine_stage_total] Skipped; mapping-driven splitting disabled"
    cp "$bundle_fasta" "$final_bundle_fasta"
fi

postprocess_options=(
    --min-bundle-coverage "$min_bundle_coverage" \
    --min-identity "$min_identity" \
    --min-mapq "$min_mapq" \
    --max-bundle-gap "$max_bundle_gap" \
    --max-genome-gap "$max_genome_gap" \
    --max-overlap-fraction "$max_overlap_fraction" \
    --fixed-bundle-min-coverage "$fixed_bundle_min_coverage" \
    --max-fixed-overlap-fraction "$max_fixed_overlap_fraction" \
    --max-fixed-overlap-bp "$max_fixed_overlap_bp" \
    --novel-fraction "$novel_fraction" \
    --bundle-alias dup1_1=dup1 \
    --bundle-alias dup1_2=dup1
)

common_split_count=0
common_deletion_dir="$output_dir/common_deletion_split"
if [[ "$common_deletion_split" -eq 1 ]]; then
    preliminary_dir="$output_dir/pre_common_deletion_mapping"
    echo "[refine 3/5] Building a preliminary class-local mapping set for common-indel MSA"
    python3 "$postprocess_script" \
        "$effective_paf" \
        --output-dir "$preliminary_dir" \
        --bundle-metadata-bed "$effective_metadata_bed" \
        "${postprocess_options[@]}"

    common_deletion_args=(
        "$final_bundle_fasta" "$genome_fasta"
        "$preliminary_dir/refined_bundle_mappings.tsv"
        --bundle-metadata-bed "$effective_metadata_bed"
        --output-dir "$common_deletion_dir"
        --mafft "$mafft_cmd" --threads "$threads"
        --min-deletion-bp "$common_deletion_min_bp"
        --min-deletion-count "$common_deletion_min_count"
        --min-deletion-frequency "$common_deletion_min_frequency"
        --event-overlap-fraction "$common_deletion_event_overlap"
        --merge-identity "$common_deletion_merge_identity"
    )
    if [[ -n "$sample_map" ]]; then
        common_deletion_args+=(--sample-map "$sample_map")
    fi
    common_deletion_args+=("--mafft-params=$mafft_params")
    if [[ -n "$mafft_extra" ]]; then
        common_deletion_args+=("--mafft-extra=$mafft_extra")
    fi
    if [[ "$keep_intermediates" -eq 1 ]]; then
        common_deletion_args+=(--keep-intermediates)
    fi
    python3 "$common_deletion_script" "${common_deletion_args[@]}"
    common_split_count="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["split_consensus_count"])' \
        "$common_deletion_dir/common_deletion_summary.json")"
    cp "$common_deletion_dir/bundle.fa" "$final_bundle_fasta"
    if [[ "$common_split_count" -gt 0 ]]; then
        effective_metadata_bed="$common_deletion_dir/bundle_metadata.bed"
        echo "Common indels split $common_split_count consensus bundle(s); deferring remap until final self-alignment"
    else
        echo "No recurrent internal indel exceeded the common-indel thresholds; reusing the current PAF"
    fi
fi

final_self_dir="$output_dir/post_split_self_alignment"
final_self_performed=0
final_self_operation_count=0
if [[ "$split_count" -gt 0 || "$common_split_count" -gt 0 ]]; then
    if [[ "$common_deletion_split" -eq 1 ]]; then
        stage_label="[refine 4/5]"
    else
        stage_label="[refine 3/4]"
    fi
    echo "$stage_label Iteratively reconciling the final post-split bundle library by self-alignment"
    final_self_args=(
        "$final_bundle_fasta" "$effective_metadata_bed"
        --output-dir "$final_self_dir"
        --minimap2 "$minimap2_cmd" --threads "$threads"
        --min-overlap-fraction "$final_self_min_overlap_fraction"
        --min-overlap-bp "$final_self_min_overlap_bp"
        --min-identity "$min_identity" --min-mapq "$min_mapq"
        --merge-identity "$common_deletion_merge_identity"
    )
    if [[ "$self_overlap_reconciliation" -eq 0 ]]; then
        final_self_args+=(--no-overlap-reconciliation)
    fi
    if [[ "$keep_intermediates" -eq 1 ]]; then
        final_self_args+=(--keep-intermediates)
    fi
    python3 "$final_self_script" "${final_self_args[@]}"
    cp "$final_self_dir/bundle.fa" "$final_bundle_fasta"
    effective_metadata_bed="$final_self_dir/bundle_metadata.bed"
    final_self_operation_count="$(python3 -c 'import json, sys; counts=json.load(open(sys.argv[1]))["counts"]; print(counts.get("reconciliation_operation_count", counts.get("split_operation_count", 0)))' \
        "$final_self_dir/self_alignment_summary.json")"
    final_self_performed=1
    effective_paf="$output_dir/bundle_to_genome.post_split_self.raw.paf"
    echo "Final self-alignment completed $final_self_operation_count correction operation(s); performing final remap"
    map_consensuses "$final_bundle_fasta" "$effective_paf" 10000
else
    if [[ "$common_deletion_split" -eq 1 ]]; then
        stage_label="[refine 4/5]"
    else
        stage_label="[refine 3/4]"
    fi
    echo "$stage_label Skipped final self-alignment; no consensus bundle was split"
fi

postprocess_inputs=("$effective_paf")
final_filter_options=()
if [[ "$filter_final_bundles" -eq 1 ]]; then
    final_filter_options=(--filter-final-bundles --genome-fasta "$genome_fasta"
        --final-min-sequence-coverage "$final_min_sequence_coverage"
        --final-max-low-coverage-fraction "$final_max_low_coverage_fraction"
        --final-max-bundle-bp "$final_max_bundle_bp")
    echo "Final class pruning enabled: sequence coverage and signed-haplotype protection"
fi
if [[ -n "$merged_bed" ]]; then
    postprocess_inputs+=("$merged_bed")
    if [[ "$common_deletion_split" -eq 1 ]]; then stage_label="[refine 5/5]"; else stage_label="[refine 4/4]"; fi
    echo "$stage_label Joining adjacent PAF blocks, normalizing classes, and comparing merged-BED coverage"
else
    if [[ "$common_deletion_split" -eq 1 ]]; then stage_label="[refine 5/5]"; else stage_label="[refine 4/4]"; fi
    echo "$stage_label Joining adjacent PAF blocks and normalizing classes; coverage comparison disabled"
fi
run_final_postprocess() {
    local inputs=("$effective_paf")
    if [[ -n "$merged_bed" ]]; then inputs+=("$merged_bed"); fi
    python3 "$postprocess_script" \
        "${inputs[@]}" \
        --output-dir "$output_dir" \
        --bundle-metadata-bed "$effective_metadata_bed" \
        --cross-bundle-locus-dedup \
        --locus-dedup-min-bundle-coverage "$locus_dedup_min_bundle_coverage" \
        "${postprocess_options[@]}" \
        ${final_filter_options[@]+"${final_filter_options[@]}"}
}

run_final_postprocess

needle_round_count=0
needle_operation_count=0
needle_root="$output_dir/post_filter_needle"
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    mkdir -p "$needle_root"
    while true; do
        needle_round_count=$((needle_round_count + 1))
        round_name="$(printf 'round_%04d' "$needle_round_count")"
        needle_round_dir="$needle_root/$round_name"
        mkdir -p "$needle_round_dir"
        cp "$output_dir/coverage_filtered_bundle_mappings.tsv" \
            "$needle_round_dir/coverage_only_mappings.tsv"
        echo "[post-filter Needle $needle_round_count] Checking genomic-overlap candidates with EMBOSS Needle"
        needle_args=(
            "$final_bundle_fasta" "$effective_metadata_bed"
            "$needle_round_dir/coverage_only_mappings.tsv"
            --output-dir "$needle_round_dir"
            --needle "$needle_cmd" --threads "$threads"
            --min-overlap-fraction "$final_self_min_overlap_fraction"
            --min-overlap-bp "$final_self_min_overlap_bp"
            --min-identity "$min_identity"
            --merge-identity "$common_deletion_merge_identity"
            --gapopen "$needle_gapopen" --gapextend "$needle_gapextend"
        )
        if [[ "$keep_intermediates" -eq 1 ]]; then
            needle_args+=(--keep-intermediates)
        fi
        python3 "$needle_refine_script" "${needle_args[@]}" \
            > "$needle_round_dir/run_result.json"
        round_operations="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["reconciliation_operation_count"])' \
            "$needle_round_dir/post_filter_needle_summary.json")"
        if [[ "$round_operations" -eq 0 ]]; then
            echo "No Needle-confirmed actionable overlap remains after coverage filtering"
            break
        fi
        needle_operation_count=$((needle_operation_count + round_operations))
        if [[ "$needle_round_count" -ge "$post_filter_needle_max_rounds" ]]; then
            echo "error: post-filter Needle refinement reached --post-filter-needle-max-rounds=$post_filter_needle_max_rounds before convergence" >&2
            exit 2
        fi
        cp "$needle_round_dir/bundle.fa" "$final_bundle_fasta"
        effective_metadata_bed="$needle_round_dir/bundle_metadata.bed"
        effective_paf="$output_dir/bundle_to_genome.post_filter_needle.${round_name}.raw.paf"
        echo "Needle reconciled one overlap; remapping the updated full bundle library"
        map_consensuses "$final_bundle_fasta" "$effective_paf" 10000
        run_final_postprocess
    done
fi

final_trf_filter_dir="$output_dir/final_trf_filter"
final_trf_filter_performed=0
final_trf_removed_class_count=0
final_trf_removed_consensus_count=0
if [[ -n "$final_trf_bed" ]]; then
    mkdir -p "$final_trf_filter_dir"
    cp "$final_bundle_fasta" "$final_trf_filter_dir/bundle.before.fa"
    cp "$effective_metadata_bed" "$final_trf_filter_dir/bundle_metadata.before.bed"
    echo "[final TRF] Requiring one merged interval to cover >=$final_trf_full_coverage of an occurrence; counting each class once per carrying sequence"
    python3 "$final_trf_filter_script" \
        "$final_trf_filter_dir/bundle.before.fa" \
        "$final_trf_filter_dir/bundle_metadata.before.bed" \
        "$output_dir/coverage_filtered_bundle_mappings.tsv" \
        "$output_dir/refined_bundle_class_ids.tsv" \
        "$final_trf_bed" \
        --output-dir "$final_trf_filter_dir" \
        --trf-min-length "$final_trf_min_length" \
        --trf-min-identity "$final_trf_min_identity" \
        --full-coverage "$final_trf_full_coverage" \
        --removal-frequency "$final_trf_removal_frequency" \
        > "$final_trf_filter_dir/run_result.json"
    final_trf_removed_class_count="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["removed_bundle_class_count"])' \
        "$final_trf_filter_dir/final_trf_filter_summary.json")"
    final_trf_removed_consensus_count="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["removed_consensus_count"])' \
        "$final_trf_filter_dir/final_trf_filter_summary.json")"
    if [[ "$final_trf_removed_class_count" -gt 0 ]]; then
        final_trf_retained_consensus_count="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["retained_consensus_count"])' \
            "$final_trf_filter_dir/final_trf_filter_summary.json")"
        if [[ "$final_trf_retained_consensus_count" -eq 0 ]]; then
            echo "error: final TRF filtering removed every bundle consensus; no final remap is possible" >&2
            exit 2
        fi
        cp "$final_trf_filter_dir/bundle.fa" "$final_bundle_fasta"
        effective_metadata_bed="$final_trf_filter_dir/bundle_metadata.bed"
        effective_paf="$output_dir/bundle_to_genome.final_trf_filtered.raw.paf"
        echo "Final TRF filter removed $final_trf_removed_class_count class(es); remapping the filtered full library"
        map_consensuses "$final_bundle_fasta" "$effective_paf" 10000
        run_final_postprocess
        final_trf_filter_performed=1
    else
        echo "Final TRF filter removed no bundle classes"
    fi
fi

python3 - "$output_dir/refinement_summary.json" "$post_filter_needle_refinement" \
    "$needle_round_count" "$needle_operation_count" "$needle_root" \
    "$final_self_min_overlap_fraction" "$final_self_min_overlap_bp" \
    "$needle_gapopen" "$needle_gapextend" <<'PY'
import json
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
enabled = bool(int(sys.argv[2]))
root = Path(sys.argv[5])
summary = json.loads(summary_path.read_text())
round_count = int(sys.argv[3])
round_summaries = ([str((root / f"round_{number:04d}" / "post_filter_needle_summary.json").resolve())
                    for number in range(1, round_count + 1)] if enabled else [])
summary["post_filter_needle_refinement"] = {
    "enabled": enabled,
    "alignment_backend": "EMBOSS needle" if enabled else None,
    "candidate_source": "coverage_filtered_bundle_mappings.tsv before any locus deduplication",
    "one_operation_then_remap": True,
    "round_count": round_count,
    "reconciliation_operation_count": int(sys.argv[4]),
    "min_overlap_fraction": float(sys.argv[6]),
    "min_overlap_bp": int(sys.argv[7]),
    "gapopen": float(sys.argv[8]),
    "gapextend": float(sys.argv[9]),
    "round_summaries": round_summaries,
}
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
PY

python3 - "$output_dir/refinement_summary.json" "$final_trf_filter_performed" \
    "$final_trf_removed_class_count" "$final_trf_removed_consensus_count" \
    "$final_trf_filter_dir/final_trf_filter_summary.json" "$final_trf_bed" <<'PY'
import json
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
summary = json.loads(summary_path.read_text())
enabled = bool(sys.argv[6])
summary["final_trf_filter"] = {
    "enabled": enabled,
    "performed_remap": bool(int(sys.argv[2])),
    "removed_bundle_class_count": int(sys.argv[3]),
    "removed_consensus_count": int(sys.argv[4]),
    "summary": str(Path(sys.argv[5]).resolve()) if enabled else None,
    "counting_unit": "distinct carrying population sequence per numeric bundle class",
    "overlap_rule": (
        "one selected merged TRF interval must cover at least full_coverage "
        "of one mapped occurrence footprint; disconnected intervals are not summed"
    ),
}
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
PY

if [[ "$common_deletion_split" -eq 1 ]]; then
    python3 - "$output_dir/refinement_summary.json" \
        "$common_deletion_dir/common_deletion_summary.json" "$common_split_count" \
        "$final_self_performed" "$needle_operation_count" \
        "$final_trf_filter_performed" <<'PY'
import json
import sys
from pathlib import Path

summary_path, deletion_path = map(Path, sys.argv[1:3])
summary = json.loads(summary_path.read_text())
deletion = json.loads(deletion_path.read_text())
summary["common_deletion_split"] = {
    "enabled": True,
    "split_consensus_count": int(sys.argv[3]),
    "parameters": deletion["parameters"],
    "counts": deletion["counts"],
    "summary": str(deletion_path.resolve()),
    "final_pgr_source": (
        "post-final-TRF-filter remap" if int(sys.argv[6]) else
        "post-filter Needle reconciliation remap" if int(sys.argv[5]) else
        "post-final-split self-alignment remap" if int(sys.argv[4]) else
        "unchanged PAF"
    ),
}
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
PY
fi

python3 - "$output_dir/refinement_summary.json" "$final_self_performed" \
    "$final_self_operation_count" "$final_self_dir/self_alignment_summary.json" \
    "$split_count" "$common_split_count" <<'PY'
import json
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
performed = bool(int(sys.argv[2]))
self_summary_path = Path(sys.argv[4])
summary = json.loads(summary_path.read_text())
summary["post_split_self_alignment"] = {
    "enabled": True,
    "performed": performed,
    "trigger": "at least one mapping/common-indel consensus was split",
    "mapping_split_consensus_count": int(sys.argv[5]),
    "common_deletion_split_consensus_count": int(sys.argv[6]),
    "correction_operation_count": int(sys.argv[3]),
    "summary": str(self_summary_path.resolve()) if performed else None,
    "final_remap_performed": performed,
    "final_remap_minimap_gap_bp": 10000 if performed else None,
}
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
PY

# set -e ensures failed mapping, splitting, or final post-processing never
# reaches cleanup. Source inputs and every non-FASTA diagnostic stay intact.
cleanup_options=(--audit-path "$output_dir/intermediate_cleanup.json")
if [[ "$split_fasta_generated" -eq 1 ]]; then
    cleanup_options+=(--generated-path "$split_dir/bundle.fa")
fi
if [[ "$split_summary_generated" -eq 1 ]]; then
    cleanup_options+=(--generated-path "$split_dir/mapping_split_summary.json")
fi
for input_path in "${input_paths[@]}"; do
    cleanup_options+=(--protected-path "$input_path")
done
if [[ "$keep_intermediates" -eq 1 ]]; then
    cleanup_options+=(--keep-intermediates)
fi
python3 "$cleanup_script" final-refinement "$output_dir" \
    "${cleanup_options[@]}"

echo "Refinement complete"
echo "Raw PAF: $raw_paf"
echo "PAF used for final refinement: $effective_paf"
echo "Final bundle consensus FASTA: $final_bundle_fasta"
if [[ "$mapping_split" -eq 1 ]]; then
    echo "Preliminary bundle BED: $split_dir/preliminary_bundles.bed"
    echo "Uncovered genome intervals: $split_dir/uncovered_genome_regions.bed"
    echo "Mapping-driven split statistics: $split_dir/bundle_split_stats.tsv"
    echo "Parent/fragment ID map: $split_dir/bundle_fragments.tsv"
fi
if [[ "$common_deletion_split" -eq 1 ]]; then
    echo "Common-indel summary: $common_deletion_dir/common_deletion_summary.json"
    echo "Common-deletion events: $common_deletion_dir/deletion_events.tsv"
    echo "Common-insertion events: $common_deletion_dir/insertion_events.tsv"
    echo "Common-indel fragments: $common_deletion_dir/bundle_fragments.tsv"
fi
if [[ "$final_self_performed" -eq 1 ]]; then
    echo "Post-split self-alignment summary: $final_self_dir/self_alignment_summary.json"
    echo "Post-split self-alignment events: $final_self_dir/self_split_events.tsv"
fi
if [[ "$post_filter_needle_refinement" -eq 1 ]]; then
    echo "Post-filter Needle rounds: $needle_round_count"
    echo "Post-filter Needle operations: $needle_operation_count"
    echo "Post-filter Needle audit root: $needle_root"
fi
if [[ -n "$final_trf_bed" ]]; then
    echo "Final TRF-filter removed classes: $final_trf_removed_class_count"
    echo "Final TRF-filter audit: $final_trf_filter_dir/final_trf_filter_summary.json"
    if [[ "$final_trf_filter_performed" -eq 1 ]]; then
        echo "Final TRF-filter remap: $effective_paf"
    fi
fi
echo "Refined mapping table: $output_dir/refined_bundle_mappings.tsv"
echo "Coverage-only pre-dedup mapping table: $output_dir/coverage_filtered_bundle_mappings.tsv"
echo "Refined mappings: $output_dir/refined_bundle_mappings.bed"
echo "PGR-compatible mappings: $output_dir/refined_bundle_mappings.pgr.bed"
if [[ "$filter_final_bundles" -eq 1 ]]; then
    echo "Before-filter PGR BED: $output_dir/refined_bundle_mappings.before_filter.pgr.bed"
    echo "Final class deletion decisions: $output_dir/final_bundle_filter.tsv"
    echo "Signed haplotypes before/after: $output_dir/final_haplotypes.before.tsv / final_haplotypes.after.tsv"
fi
echo "Bundle class/ID map: $output_dir/refined_bundle_class_ids.tsv"
if [[ -n "$merged_bed" ]]; then
    echo "Novel-coverage report: $output_dir/novel_bundle_coverage.tsv"
    echo "Novel mappings: $output_dir/novel_bundle_mappings.bed"
fi
echo "Same-locus alias decisions: $output_dir/duplicate_same_locus_mappings.tsv"
echo "Run summary: $output_dir/refinement_summary.json"
echo "Intermediate cleanup report: $output_dir/intermediate_cleanup.json"
