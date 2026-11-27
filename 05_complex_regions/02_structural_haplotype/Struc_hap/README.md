# Structural haplotype refinement from PGR-TK bundles

Refine population PGR-TK bundle decompositions into a shared consensus library
and a new PGR-compatible BED for structural haplotype analysis.

The main entry point is **[`run_pgrtk_remap_pipeline.sh`](run_pgrtk_remap_pipeline.sh)**.
It starts from an original PGR-TK BED and the corresponding population sequences,
builds a library in haplotype-frequency order, and refines bundle boundaries by
remapping, repeat-aware processing, multiple alignment, and overlap reconciliation.

## Contents

- [Requirements](#requirements)
- [Input files](#input-files)
- [Quick start](#quick-start)
- [Workflow](#workflow)
- [External software commands](#external-software-commands)
- [Pipeline parameter defaults](#pipeline-parameter-defaults)
- [MAFFT configuration and quick mode](#mafft-configuration-and-quick-mode)
- [Outputs and logs](#outputs-and-logs)
- [Restarting and resuming](#restarting-and-resuming)
- [Optional short-bundle filtering](#optional-short-bundle-filtering)
- [Other workflows](#other-workflows)
- [Validation](#validation)

## Requirements

- Python **3.9 or later** (`python3`)
- [minimap2](https://github.com/lh3/minimap2)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [EMBOSS Needle](https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html)

[PGR-TK](https://github.com/GeneDx/pgr-tk) and Tandem Repeats Finder ([TRF](https://tandem.bu.edu/trf/trf.html)) are **upstream data-generation tools**. TRF annotations
are optional.

Use `--mafft /path/to/mafft`, `--minimap2 /path/to/minimap2`, or
`--needle /path/to/needle` to select executables.

## Input files

- PGR-TK BED: `--pgr-bed`

Use the native four-column bundle decomposition. 

```text
# path_id  start  end  bundle_id:bundle_size:direction:bundle_pos_start:bundle_pos_end:R|U
sampleA_h1	0	12000	1:100:0:0:99:U
sampleA_h1	12000	18000	2:50:1:49:0:U
sampleB_h1	0	12500	1:100:0:0:99:U
```

- Population FASTA: `--fasta`

Provide the population sequence of the same region. 

- TRF annotations: `--trf-bed` (optional)

Supply a BED conversion of TRF annotations. At least seven columns are required:

```text
# path_id  start  end  repeat_id  period  copy_number  identity
sampleA_h1	2000	4500	tr1	25	100	92
```
  Omitting `--trf-bed` skips TRF processing.

## Quick start

Run from the repository root, replacing the example paths with your files:

```bash
bash run_pgrtk_remap_pipeline.sh \
  --pgr-bed data/original.pgr.bed \
  --fasta data/population.fa \
  --output-dir results/region1 \
  --threads 8
```

To enable TRF processing and count diploid paths by sample:

```bash
bash run_pgrtk_remap_pipeline.sh \
  --pgr-bed data/original.pgr.bed \
  --fasta data/population.fa \
  --trf-bed data/population.trf.bed \
  --output-dir results/region1 \
  --threads 8 \
  > region1.pipeline.log 2>&1
```

Use a separate output directory for each region or parameter configuration.
The primary results are:

- `results/region1/pgrtk_remapped.bed`: final PGR-compatible bundle mappings.
- `results/region1/03_refined_bundle_mapping/final_bundle.fa`: final refined
  consensus library corresponding to those mappings.

The FASTA in `02_bundle_consensus/` is the **pre-final-refinement** library.
Use the `03_` FASTA for the final result.

## Workflow

**Step 1 — `01_bundle_alignments/`: group and rank haplotypes.** Optional TRF
filtering removes repetitive bundle classes before grouping. Paths with the
same ordered list of signed bundle IDs form a preliminary haplotype; orientation
and repeated occurrences are retained. Despite the directory name, this step
writes ranking information rather than performing source MSAs.

**Step 2 — `02_bundle_consensus/`: build the library iteratively.** For rank 1,
align occurrences at each bundle position, refine candidates with TRF when
available, and generate consensus alleles. Map the library to rank-1 sequences
for splitting only, then reconcile overlaps between bundles by self-alignment.
For each subsequent rank, map the existing library, evaluate supported truncations
using the current and previously processed ranks, split as needed, and remap to
the current rank. Subtract this current-rank coverage from the original bundle
occurrences. Refine only the uncovered fragments as candidate additions, then
self-align the updated library after additions. Inherited cuts prevent later
ranks from silently recreating an unsplit parent.

**Step 3 — `03_refined_bundle_mapping/`: refine population-wide mappings.**
[`scripts/refine_v2.sh`](scripts/refine_v2.sh) maps the library to all population
sequences, chains compatible partial PAF blocks, evaluates mapping-supported
cuts, and aligns retained occurrences per bundle class to detect common internal
deletions and insertions. A changed library undergoes self-alignment and remapping.
Coverage-qualified overlapping mappings nominate pairs for Needle confirmation;
actionable overlaps are reconciled, followed by remapping and another check.
A final TRF class filter runs when annotations are provided, and a changed
library is remapped before publishing the final BED.

Shared intervals with identity **>0.90** retain one representative by default;
intervals at or below that threshold remain distinct alleles of the **same
numeric bundle class**. Substantial nonshared flanks become separate bundles.

Final locus redundancy is resolved by two local maximum-selection passes at
fraction-either overlap `>=0.80`, first scoring matching bases, then identity.
Needle candidates are captured **before** this locus competition, so the
competition cannot hide an overlap from Needle. Mapping deduplication differs
from optional short-bundle class deletion, which this wrapper disables.

## External software commands

These are the arguments explicitly passed by the current pipeline, rather than
the programs' own defaults. Other options retain the installed program's
defaults. `THREADS` is `8` unless changed with `--threads`.

### minimap2

| Use | Explicit arguments |
| --- | --- |
| Rank-to-population mapping, including historical split evidence | `-x asm20 -p 0.2 -N 10 -g 500 -c -t THREADS` |
| Self-alignment at rank 1, after additions, and after final splits | `-x asm20 -X -g 500 -c -t THREADS` |
| Quick-mode source-to-reference placement | `-x asm20 -g 500 --secondary=no -c -t THREADS` |
| Initial step-3 population mapping and final population remaps | `-x asm20 -p 0.2 -N 10 -g 10000 -c -t THREADS` |
| Temporary step-3 remap immediately after mapping-driven splitting | `-x asm20 -p 0.2 -N 10 -g 500 -c -t THREADS` |

For
population remapping the actual file order is:

```bash
minimap2 -x asm20 -p 0.2 -N 10 -g 10000 -c -t 8 \
  bundle.fa population.fa > bundle_to_genome.raw.paf
```

The bundle library is the **target** (PAF columns 6–9); population sequences are
the **query** (columns 1–4). Downstream code converts these axes into genomic
BED coordinates. Preserve this order when supplying PAF to standalone helpers.
The `-g` settings are fixed by stage, not exposed as a wrapper option.

### MAFFT

Ordinary step-2 MSAs and step-3 common-indel MSAs both default to:

```bash
mafft --auto --thread 8 --inputorder input.fa > input.aln.fa
```

Quick-mode local additions also use `--addfragments`; valid partially mapped
records use `--add` against the assembled backbone MSA. `--mafft-02-params` and
`--mafft-03-params` replace `--auto` in their respective stages. `--mafft-extra`
appends arguments to both. See the
[MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html) for
alignment strategies; quick mode is a pipeline option described below.

### EMBOSS Needle

The default final overlap-confirmation command is:

```bash
needle -auto -asequence first.fa -bsequence second.fa \
  -gapopen 10.0 -gapextend 0.5 -endweight N \
  -aformat3 fasta -outfile alignment.fa
```

The second sequence is oriented according to the overlap evidence before
alignment. End gaps are unweighted. The wrapper does not expose Needle scoring
options. Mapping pairs are nominated when their genomic overlap exceeds 20%
of either mapping or 10,000 bp; Needle-confirmed sequence overlaps must also
pass the bundle-level 20%/10,000-bp gate. These gates follow the wrapper's
`--rank1-self-min-overlap-fraction` and `--rank1-self-min-overlap-bp` settings.
Direct `scripts/refine_v2.sh` calls support `--needle-gapopen` and
`--needle-gapextend`. Its convergence safeguard defaults to 100 rounds via
`--post-filter-needle-max-rounds`; reaching the limit is an error rather than a
successful final result. Needle is not controlled by `--threads`.

## Pipeline parameter defaults

Unless marked as percent, fractions range from `0` to `1`; lengths are in bp.
Use `bash run_pgrtk_remap_pipeline.sh --help` for the full interface.

### Consensus and self-overlap reconciliation

| Option | Default | Meaning |
| --- | --- | --- |
| `--threads` | `8` | Threads passed to MAFFT/minimap2 |
| `--min-column-occupancy` | `0.50` | Minimum occupied fraction for a consensus column |
| `--min-base-frequency` | `0.50` | Minimum majority-base frequency for a consensus call |
| `--merge-identity` | `0.90` | Shared-region identity must be strictly greater to collapse alleles |
| `--addition-min-bp` | `1` | Minimum uncovered source fragment considered for addition |
| `--rank1-self-min-overlap-fraction` | `0.20` | Self-overlap union must exceed this fraction of either bundle, or pass the bp gate |
| `--rank1-self-min-overlap-bp` | `10000` | Alternative strict self-overlap union-length gate |

Despite their names, `--rank1-self-*` thresholds also configure post-addition
and post-final-split self-alignment. An outer self-alignment cut is suppressed
only when the terminal fragment is **<200 bp AND <5%** of its current parent.
Internal pieces have no minimum length. This differs from the mapping-based
terminal guard below.

### PAF chaining and coverage

| Option | Default | Meaning |
| --- | --- | --- |
| `--min-bundle-coverage` | `0.60` | Minimum union bundle coverage for accepted ranked/preliminary chains |
| `--locus-dedup-min-bundle-coverage` | `0.70` | Final completed-chain bundle coverage threshold |
| `--min-identity` | `0.0` | Minimum PAF block identity; distinct from the allele-collapse threshold |
| `--min-mapq` | `0` | Minimum PAF block MAPQ |
| `--max-bundle-gap` | `1000` | Maximum gap/overlap on the bundle axis for chaining |
| `--max-genome-gap` | `1000` | Maximum gap/overlap on the genomic axis for chaining |
| `--max-overlap-fraction` | `0.10` | Maximum overlap on either axis, divided by bundle consensus length |
| `--fixed-bundle-min-coverage` | `0.80` | Fix sufficiently complete single PAF blocks and exclude them from chaining |
| `--max-fixed-overlap-fraction` | `0.20` | Maximum proposed-chain overlap divided by the other fixed class's consensus length |
| `--max-fixed-overlap-bp` | `10000` | Alternative maximum overlap with another fixed class |

Compatible partial blocks are **chained before final coverage filtering**.
They must share a consensus coordinate system, genomic path, strand, and
collinear order. A proposed chain is rejected when overlap with a different
fixed class exceeds either fixed-overlap limit. Coverage uses the **union of
aligned bundle intervals divided by full consensus length**, not the bounding
span or the sum of overlapping hits. The final 70% rule does not replace the
ranked 60% rule; raising the general threshold can also make final acceptance
stricter.

### Mapping-supported splitting

| Option | Default | Meaning |
| --- | --- | --- |
| `--mapping-split-min-truncated-count` | `5` | Required gap-supported short-occurrence count for a bundle |
| `--mapping-split-truncated-fraction` | `0.60` | Short occurrence length fraction must be strictly below this value |
| `--mapping-split-max-unmapped-bp` | `10000` | Report longer unaligned bundle intervals; missing coverage alone does not trigger a cut |
| `--mapping-split-cut-merge-bp` | `100` | Cluster nearby candidate cut positions |
| `--mapping-split-min-cut-support` | `1` | Per-cut support, in addition to the bundle-level count gate |
| `--mapping-split-min-terminal-bp` | `1000` | Suppress cuts producing shorter terminal pieces; internal pieces are retained |
| `--mapping-split-gap-min-mapping-bp` | `10000` | Each selected PAF must cover strictly more than this many bp inside a preliminary-BED gap |
| `--mapping-split-gap-max-overlap-bp` | `1000` | Selected alignments must overlap by strictly less than this; `0` requires no overlap |

Candidate evidence comes from gaps in the preliminary accepted BED. Within
each gap, select compatible alignments maximizing covered genomic union and
look for incomplete bundle mappings. Aggregate population evidence before
choosing cuts. Long uncovered bundle intervals **without PAF evidence are
reported, not automatically split**. Truncation decisions use alignment endpoints
rather than CIGAR deletion operations.

### TRF processing

| Option | Default | Meaning |
| --- | --- | --- |
| `--trf-min-length` | `0` | No minimum annotation-length filter |
| `--trf-min-identity` | `80` | Minimum TRF annotation identity, in percent |
| `--trf-full-coverage` | `0.90` | Initial/candidate occurrence coverage by one merged TRF interval |
| `--final-trf-full-coverage` | `0.95` | Final mapped occurrence coverage by one merged TRF interval |
| `--trf-bundle-removal-frequency` | `0.50` | Whole-class/candidate removal frequency; denominators differ by stage below |
| `--trf-split-min-bp` | `1000` | Minimum longest ungapped span of a population TRF-region MSA component |
| `--trf-min-sequence-coverage` | `0.50` | Required sequence span relative to the longest span to count as repeat-region presence |
| `--trf-min-population-presence` | `0.50` | Presence below this threshold indicates structural divergence |
| `--trf-min-length-divergence` | `0.05` | Minimum mean pairwise normalized repeat-length divergence |
| `--trf-min-nucleotide-divergence` | `0.05` | Minimum non-gap nucleotide disagreement fraction in the repeat-region MSA |
| `--trf-min-fragment-sequence-bp` | `1000` | Real bases required for a sequence to support a retained child |
| `--trf-min-fragment-count` | `2` | Minimum child support count |
| `--trf-min-fragment-frequency` | `0.02` | Minimum child support fraction |

Selected overlapping/touching TRF annotations are merged on each genomic path.
For **whole-bundle removal**, a single merged interval must meet the coverage
threshold; disconnected intervals are never summed to reach it.

- **Before grouping:** remove every occurrence of an original class when at
  least 50% of its occurrences individually reach 90% TRF coverage.
- **Rank-1 and later candidate refinement:** apply the corresponding candidate
  occurrence rule; otherwise project TRF intervals into the source MSA and
  union their positions across occurrences. A sufficiently long region is
  removed/split if its length divergence, nucleotide divergence, or low population
  presence meets the respective criterion. Child support must reach
  `max(2, ceil(0.02 × occurrence_count))`, with at least 1000 real bases per
  supporting sequence. Repeat-region scores are not whole-bundle scores.
- **Final generated library:** among distinct sequences carrying a
  coverage-qualified mapping of a class, a sequence is positive if **any one**
  occurrence reaches 95% coverage by one merged interval. Remove the entire
  class when at least 50% of these carrier sequences are positive. Multiple
  copies in one sequence contribute only one vote. Any overlap alone is
  insufficient. Remap the retained library after removal.

In the final TRF filter, occurrence coverage uses the union of its mapped
genomic footprint as the denominator; gaps between chained blocks are excluded.
Repeat-length divergence is the mean of `abs(a - b) / max(a, b)` across pairs
of ungapped sequence lengths in the projected TRF-region MSA, with two zero
lengths contributing zero.

### Final common deletions and insertions

| Option | Default | Meaning |
| --- | --- | --- |
| `--common-deletion-min-bp` | `1000` | Event length must be **strictly greater**; applies to insertions too |
| `--common-deletion-min-count` | `5` | Minimum distinct supporting population units |
| `--common-deletion-min-frequency` | `0.01` | Minimum supporting fraction among units mapped to the class |
| `--common-deletion-event-overlap` | `0.80` | Reciprocal event overlap used to cluster shifted endpoints |

Required support is `max(5, ceil(0.01 × mapped_population_units))`, using samples
with `--sample-map` and paths otherwise. Only internal events with sequence on
both sides qualify. A deletion splits off the variable middle sequence as a
separate child; it is not simply discarded. A recurrent insertion relative to
all consensus alleles creates a cut at its consensus anchor; inserted bases
absent from the library are not invented as a consensus.

## MAFFT configuration and quick mode

For example, select progressive strategies independently in the two MSA stages:

```bash
bash run_pgrtk_remap_pipeline.sh \
  --pgr-bed data/original.pgr.bed \
  --fasta data/population.fa \
  --output-dir results/region1_fast \
  --mafft-02-params '--retree 2 --maxiterate 0' \
  --mafft-03-params '--retree 1 --maxiterate 0' \
  --threads 8
```

Changing alignment strategy may change inferred consensuses and cuts. For a
different step-2 source-alignment approach, add `--mafft-02-quick` (off by
default). It applies to **all source MSAs**, without a 100-kb trigger:

1. Select the longest original occurrence as reference; occurrences at least
   70% of its length form the backbone candidate set.
2. Exclude whole truncated records with length **<=10% of the reference AND
   <1000 bp**, including inherited cut templates. Exactly 1000 bp is not excluded.
3. Map remaining records to the reference with minimap2. Discard records without
   a valid forward, coordinate-consistent, CIGAR-backed placement.
4. Records with at least 95% of their own length mapped are projected into
   reference windows. Align backbone windows, then locally add shorter records
   using MAFFT `--addfragments`. The default maximum reference window width is
   `--mafft-02-window-bp 10000`; placement boundaries can make windows smaller.
5. Join the windows, then use MAFFT `--add` for records with a valid placement
   covering less than 95% of their own length. These records are not dropped
   solely for low mapping coverage.

Quick mode changes step-2 source MSAs; it does **not** window step-2 pairwise
candidate/allele comparisons or step-3 common-indel MSAs. For step-2
source/comparison sets whose longest sequence is **<100 bp**, a lightweight
built-in global alignment replaces MAFFT and retains the bundle. The quick-mode
short-record exclusion still applies. The final common-indel pass skips classes
whose longest consensus is <100 bp.

## Optional short-bundle filtering

To apply
coverage- and haplotype-aware filtering afterwards, use
[`scripts/filter_pgrtk_bed.py`](scripts/filter_pgrtk_bed.py):

```bash
# Optional upstream indexing tool: samtools
samtools faidx data/population.fa

python3 scripts/filter_pgrtk_bed.py \
  results/region1/pgrtk_remapped.bed data/population.fa.fai \
  --output-bed results/region1/pgrtk_remapped.filtered.bed \
  --max-bundle-bp 500 \
  --min-sequence-coverage 0.95 \
  --max-low-coverage-fraction 0.05
```

This filter considers classes with median occurrence span **<=500 bp**, shortest
first. Removal must change the haplotype partition through allowed indel
differences while preserving protected signed structural distinctions and the
sequence-coverage criterion. A class is retained if its removal leaves the
partition unchanged. Coverage is measured against full FAI sequence lengths;
optional `--trf-bed` adds immutable TRF interval coverage to the retained bundle
union. The default permits at most 5% of sequences below 95% coverage. Filtering
is a sequential single pass, not automatic iteration to convergence. It writes
a BED and decision audits; it does not alter the consensus FASTA.
