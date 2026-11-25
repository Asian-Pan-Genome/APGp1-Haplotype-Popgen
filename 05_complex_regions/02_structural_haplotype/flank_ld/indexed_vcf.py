"""On-demand indexed VCF/BCF retrieval using the external bcftools CLI."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import threading
from itertools import chain
from pathlib import Path
from typing import Dict, Mapping, Sequence

try:
    from .flanking_block_ld import Region, VCFData, merge_candidate_intervals, normalize_chrom, read_vcf
except ImportError:
    from flanking_block_ld import Region, VCFData, merge_candidate_intervals, normalize_chrom, read_vcf  # type: ignore


class VCFQueryError(OSError):
    """Fatal I/O/parser error, never a rejectable background-region QC error."""


class BCFtoolsVCFSource:
    QC_COUNTERS = ("candidate_variant_count", "filtered_non_snp_count", "filtered_vcf_filter_count",
                   "filtered_maf_count", "filtered_call_rate_count", "filtered_unphased_count")

    def __init__(self, args: argparse.Namespace, logger):
        self.args, self.logger = args, logger
        self.path = args.vcf.resolve()
        self.executable = shutil.which(getattr(args, "bcftools", "bcftools"))
        if self.executable is None:
            raise VCFQueryError("bcftools executable not found; install bcftools or specify --bcftools PATH; "
                                "no automatic full-VCF preload fallback")
        if not self.path.is_file():
            raise VCFQueryError(f"VCF/BCF does not exist: {self.path}")
        indexes = [Path(str(self.path) + suffix) for suffix in (".csi", ".tbi")]
        if not any(index.is_file() for index in indexes):
            raise VCFQueryError("on-demand access requires a BGZF VCF/BCF with a .csi or .tbi index; "
                                "prepare it with bcftools index FILE; plain VCF/gzip cannot be randomly accessed; "
                                "use --vcf-access preload only if a full scan is acceptable")
        header = self._capture(["view", "--header-only", str(self.path)])
        columns = [line for line in header.splitlines() if line.startswith("#CHROM\t")]
        if len(columns) != 1 or len(columns[0].split("\t")) < 10:
            raise VCFQueryError("bcftools returned no valid sample-bearing #CHROM header")
        self.sample_header = columns[0] + "\n"
        samples = columns[0].split("\t")[9:]
        if len(set(samples)) != len(samples):
            raise VCFQueryError("VCF has duplicate sample names")
        # --stats reads the index, not the variant records; --all includes empty contigs.
        index_stats = self._capture(["index", "--stats", "--all", str(self.path)])
        self.contigs: Dict[str, str] = {}
        for line in index_stats.splitlines():
            if not line.strip():
                continue
            actual = line.split("\t")[0]
            normalized = normalize_chrom(actual)
            if normalized in self.contigs and self.contigs[normalized] != actual:
                raise VCFQueryError(f"ambiguous chromosome normalization: {self.contigs[normalized]!r}, {actual!r}")
            self.contigs[normalized] = actual
        self.prefiltered = getattr(args, "vcf_prefiltered", False)
        self.data = VCFData(len(samples), {}, {}, {}, 0, 0, 0, 0, 0, 0, self.prefiltered)
        self.query_batches = self.query_intervals = self.records_returned = self.retained_snp_observations = 0
        self.qc_counts = {field: 0 for field in self.QC_COUNTERS}
        self.missing_contigs = set()
        self._state_lock = threading.RLock()
        logger.log(f"VCF INDEX READY samples={len(samples)} contigs={len(self.contigs)} "
                   f"bcftools={self.executable} mode=on_demand; no background/whole-VCF preload")

    def _capture(self, arguments: Sequence[str]) -> str:
        result = subprocess.run([self.executable, *arguments], capture_output=True, text=True, check=False)
        if result.returncode:
            raise VCFQueryError(f"bcftools {' '.join(arguments[:3])} failed (exit {result.returncode}): "
                                f"{result.stderr[:8000].strip()}")
        if result.stderr.strip():
            self.logger.log("bcftools: " + result.stderr[:4000].strip(), "WARNING")
        return result.stdout

    def fetch(self, requested: Sequence[Region]) -> VCFData:
        """Return only requested positions; merge overlaps before querying."""
        merged = merge_candidate_intervals(requested, 0, include_region_interiors=True)
        intervals = []
        for chrom, spans in sorted(merged.items()):
            if chrom not in self.contigs:
                with self._state_lock:
                    first_warning = chrom not in self.missing_contigs
                    self.missing_contigs.add(chrom)
                if first_warning:
                    self.logger.log(f"VCF requested chromosome {chrom!r} is absent from the index; "
                                    "its windows contain no available SNPs", "WARNING")
                continue
            intervals.extend(Region(chrom, self.contigs[chrom], start, end, "query", 0)
                             for start, end in spans)
        if not intervals:
            return self.data
        with self._state_lock:
            self.query_batches += 1
            self.query_intervals += len(intervals)
            batch = self.query_batches
        self.logger.periodic("vcf_query", f"VCF QUERY START batch={batch} intervals={len(intervals)} "
                             f"requested_bp={sum(item.end - item.start for item in intervals)}")
        final_counts = {}

        def report(counts: Mapping[str, int]) -> None:
            final_counts.update(counts)
            missing = "NA" if self.prefiltered else str(counts["filtered_missing"])
            with self._state_lock:
                records_returned = self.records_returned
            self.logger.periodic("vcf_query", f"VCF QUERY PROGRESS batch={batch} intervals={len(intervals)} "
                                 f"records_returned={records_returned + counts['records_scanned']} "
                                 f"batch_records={counts['records_scanned']} "
                                 f"batch_retained_snps={counts['retained_snps']} filtered_missing={missing}" +
                                 (" maf_missing_qc=skipped" if self.prefiltered else ""))

        # .bed explicitly selects 0-based half-open coordinates and supports
        # contig names that cannot safely be represented in a -r string.
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", encoding="utf-8") as bed, \
                tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stderr:
            for region in intervals:
                bed.write(f"{region.display_chrom}\t{region.start}\t{region.end}\n")
            bed.flush()
            command = [self.executable, "view", "--no-header", "--output-type", "v",
                       "--regions-overlap", "0", "--regions-file", bed.name, str(self.path)]
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr, text=True)
            try:
                try:
                    data = read_vcf(
                        self.path, intervals, 0, self.args.min_maf, 1 - self.args.max_site_missing_rate,
                        self.args.include_filtered, self.args.allow_unphased_missing,
                        include_region_interiors=True, progress_callback=report,
                        skip_maf_missing_qc=self.prefiltered,
                        vcf_lines=chain([self.sample_header], process.stdout),
                    )
                except (ValueError, UnicodeError) as error:
                    raise VCFQueryError(f"VCF query batch {batch}: invalid queried genotype/record: {error}") from error
                code = process.wait()
                stderr.seek(0)
                message = stderr.read(8000).strip()
                if code:
                    raise VCFQueryError(f"bcftools view query batch {batch} failed (exit {code}): {message}")
                if message:
                    self.logger.log("bcftools: " + message[:4000], "WARNING")
            finally:
                if process.poll() is None:
                    process.kill()
                process.wait()
                process.stdout.close()
        with self._state_lock:
            self.records_returned += final_counts.get("records_scanned", 0)
            self.retained_snp_observations += sum(len(values) for values in data.variants_by_chrom.values())
            for field in self.QC_COUNTERS:
                self.qc_counts[field] += getattr(data, field)
        return data

    def statistics(self) -> Dict[str, object]:
        with self._state_lock:
            return dict(mode="bcftools", bcftools=self.executable, query_batches=self.query_batches,
                        query_intervals=self.query_intervals, records_returned=self.records_returned,
                        retained_snp_observations=self.retained_snp_observations,
                        qc_counts=dict(self.qc_counts),
                        count_scope="queried record observations; overlapping queries/cache eviction may count a site again",
                        missing_requested_contigs=sorted(self.missing_contigs))
