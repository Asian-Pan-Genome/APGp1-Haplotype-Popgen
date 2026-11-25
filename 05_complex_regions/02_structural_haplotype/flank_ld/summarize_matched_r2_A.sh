#!/usr/bin/env bash
# Equal-weight matched-background superiority A for *_flanking_r2.*.tsv files.
set -euo pipefail

usage() {
  echo 'Usage: bash summarize_matched_r2_A.sh INPUT_DIR OUTPUT_PREFIX [--include-incomplete]' >&2
}
if (( $# < 2 )); then usage; exit 2; fi
input_dir=$1
output_prefix=$2
shift 2
include_incomplete=0
while (( $# )); do
  case "$1" in
    --include-incomplete) include_incomplete=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -d "$input_dir" ]] || { echo "Input directory not found: $input_dir" >&2; exit 2; }
mkdir -p "$(dirname "$output_prefix")"
per_region="${output_prefix}.per_region.tsv"
overall="${output_prefix}.overall.tsv"
printf 'id\tstatus\tstatistic\ttarget_r2\tm_i\tbackground_less\tbackground_equal\tbackground_greater\tA_i\n' > "$per_region"

shopt -s nullglob
files=("$input_dir"/*_flanking_r2.summary.tsv)
if (( ${#files[@]} == 0 )); then
  echo "No *_flanking_r2.summary.tsv files in $input_dir" >&2
  exit 2
fi
included=0
skipped=0
for summary in "${files[@]}"; do
  base=${summary%.summary.tsv}
  id=${base##*/}
  id=${id%_flanking_r2}
  backgrounds="${base}.backgrounds.tsv"
  if [[ ! -f "$backgrounds" ]]; then
    echo "WARNING: missing backgrounds for $id" >&2
    ((skipped+=1))
    continue
  fi
  # Parse named columns, tolerating CRLF, UTF-8 BOM, and blank lines.
  meta=$(awk -F '\t' '
    {sub(/\r$/, "", $0); if($0 ~ /^[[:space:]]*$/) next}
    !header {
      sub(/^[^A-Za-z_]+/, "", $1)
      for(i=1;i<=NF;i++) col[$i]=i
      if(!("median_r2" in col) || !("status" in col)) {
        print "Missing median_r2 or status column in " FILENAME > "/dev/stderr"; exit 3
      }
      header=1; next
    }
    {print $col["status"] "\t" $col["median_r2"]; found=1; exit}
    END {if(!found) {print "No data row in " FILENAME > "/dev/stderr"; exit 3}}
  ' "$summary") || { echo "Invalid summary: $summary" >&2; exit 2; }
  IFS=$'\t' read -r status target <<< "$meta"
  if [[ "$include_incomplete" == 0 && "$status" != ok ]]; then
    ((skipped+=1)); continue
  fi
  if [[ "$target" == NA || -z "$target" ]]; then
    ((skipped+=1)); continue
  fi
  counts=$(awk -F '\t' -v target="$target" -v allow_incomplete="$include_incomplete" '
    {sub(/\r$/, "", $0); if($0 ~ /^[[:space:]]*$/) next}
    !header {
      sub(/^[^A-Za-z_]+/, "", $1)
      for(i=1;i<=NF;i++) col[$i]=i
      if(!("median_r2" in col)) {
        print "Missing median_r2 column in " FILENAME > "/dev/stderr"; exit 3
      }
      header=1; next
    }
    {
      value=$col["median_r2"]
      if(value=="NA" || value=="") next
      if(allow_incomplete==0 && ("full_flanks" in col) && $col["full_flanks"]!="True") next
      if(value+0 < target+0) less++
      else if(value+0 > target+0) greater++
      else equal++
    }
    END {if(!header) {print "No header in " FILENAME > "/dev/stderr"; exit 3}
         print less+0 "\t" equal+0 "\t" greater+0}
  ' "$backgrounds") || { echo "Invalid backgrounds: $backgrounds" >&2; exit 2; }
  IFS=$'\t' read -r less equal greater <<< "$counts"
  m=$((less + equal + greater))
  if (( m == 0 )); then ((skipped+=1)); continue; fi
  a=$(awk -v less="$less" -v equal="$equal" -v m="$m" 'BEGIN {printf "%.12g", (less+0.5*equal)/m}')
  printf '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n' \
    "$id" "$status" "median_r2" "$target" "$m" "$less" "$equal" "$greater" "$a" >> "$per_region"
  ((included+=1))
done
if (( included == 0 )); then
  echo "No usable target/background pairs; check status or use --include-incomplete for descriptive analysis" >&2
  exit 2
fi
awk -F '\t' -v skipped="$skipped" 'BEGIN {OFS="\t"; print "n_regions","skipped_regions","total_background_values","A","fraction_A_above_0.5"}
  NR>1 {n++; m+=$5; a+=$9; if($9>0.5) above++}
  END {print n,skipped,m,a/n,above/n}' "$per_region" > "$overall"
echo "Wrote $per_region and $overall; included=$included skipped=$skipped" >&2
