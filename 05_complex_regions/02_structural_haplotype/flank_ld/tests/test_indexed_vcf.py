import contextlib
import io
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld import indexed_vcf
from flanking_block_ld import flanking_block_ld as reader
from flanking_block_ld.flanking_block_ld import Region
from flanking_block_ld.tests import test_background_block_ld as existing_tests
from flanking_block_ld.tests import test_background_progress as progress_tests
from flanking_block_ld.tests import test_background_region_pair_p as skip_tests


class FakeBCFtools:
    """Record the external CLI contract; never used by production code."""

    def __init__(self, content, index_stats=None):
        lines = content.splitlines(keepends=True)
        self.header = "".join(line for line in lines if line.startswith("#"))
        self.records = [line for line in lines if not line.startswith("#") and line.strip()]
        contigs = sorted({line.split("\t")[0] for line in self.records})
        self.index_stats = index_stats if index_stats is not None else "".join(f"{chrom}\t1000\t100\n" for chrom in contigs)
        self.calls, self.queries = [], []
        self.fail_next = False

    def run(self, command, **kwargs):
        self.calls.append(command)
        if command[1:3] == ["view", "--header-only"]:
            return subprocess.CompletedProcess(command, 0, self.header, "")
        if command[1:4] == ["index", "--stats", "--all"]:
            return subprocess.CompletedProcess(command, 0, self.index_stats, "")
        raise AssertionError(command)

    def popen(self, command, **kwargs):
        self.calls.append(command)
        assert command[1:6] == ["view", "--no-header", "--output-type", "v", "--regions-overlap"]
        assert command[6] == "0"
        bed = Path(command[command.index("--regions-file") + 1])
        assert bed.suffix == ".bed"
        spans = [(chrom, int(start), int(end))
                 for chrom, start, end in (line.split("\t") for line in bed.read_text().splitlines())]
        self.queries.append(spans)
        content = "".join(line for line in self.records if any(
            line.split("\t")[0] == chrom and start <= int(line.split("\t")[1]) - 1 < end
            for chrom, start, end in spans))
        code = 7 if self.fail_next else 0
        self.fail_next = False
        if code:
            kwargs["stderr"].write("broken index/read failure")
            kwargs["stderr"].flush()
        process = mock.Mock()
        process.stdout = io.StringIO(content)
        process.poll.return_value = code
        process.wait.return_value = code
        return process

    @contextlib.contextmanager
    def patches(self):
        with mock.patch.object(indexed_vcf.shutil, "which", return_value="/fake/bcftools"), \
                mock.patch.object(indexed_vcf.subprocess, "run", self.run), \
                mock.patch.object(indexed_vcf.subprocess, "Popen", self.popen):
            yield


class IndexedVCFTests(unittest.TestCase):
    def test_cli_default_and_explicit_preload_and_cache_validation(self):
        argv = ["target.bed", "background.bed", "cohort.vcf.gz", "-o", "output"]
        args = workflow.build_parser().parse_args(argv)
        self.assertEqual(args.vcf_access, "bcftools")
        self.assertEqual(args.bcftools, "bcftools")
        self.assertEqual(args.vcf_window_cache_size, 8192)
        self.assertEqual(args.pair_batch_size, 1000)
        self.assertEqual(args.region_batch_size, 1000)
        self.assertEqual(workflow.build_parser().parse_args(argv + ["--vcf-access", "preload"]).vcf_access, "preload")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            workflow.main(argv + ["--vcf-window-cache-size", "1"])
        for flag in ("--pair-batch-size", "--region-batch-size"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                workflow.main(argv + [flag, "0"])

    @staticmethod
    def boundary_content():
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        return header + "".join(f"chr1\t{pos}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
                                for pos in (1, 10, 11, 20, 21, 95, 105, 195, 205, 215, 305))

    def test_indexed_fetch_normalizes_names_merges_overlaps_and_never_opens_vcf(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf.write_text(self.boundary_content())
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(self.boundary_content())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with fake.patches(), mock.patch.object(reader, "open_text", side_effect=AssertionError("full scan")):
                    source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                    self.assertEqual(source.query_batches, 0)  # Header/index only at startup.
                    data = source.fetch([Region("1", "1", 0, 15, "a", 0),
                                         Region("1", "1", 10, 20, "b", 0)])
                self.assertEqual(fake.queries, [[("chr1", 0, 20)]])
                self.assertEqual([variant.pos for variant in data.variants_by_chrom["1"]], [1, 10, 11, 20])
                self.assertEqual(source.records_returned, 4)
                self.assertEqual(source.query_intervals, 1)
            finally:
                logger.close()

    def test_region_batches_flanks_caches_selected_blocks_and_releases_raw_data(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(args.vcf.read_text())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with fake.patches():
                    source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                    store = workflow.BlockStore(source.data, args, logger, source)
                    region = Region("1", "1", 300, 400, "test", 0)
                    rows = store.region_pairs(region)
                    self.assertEqual(len(fake.queries), 1)
                    self.assertEqual(fake.queries[0], [("1", 280, 300), ("1", 400, 420)])
                    store.region_entropy_profile(region, rows)
                    store.region_pairs(region)
                    self.assertEqual(len(fake.queries), 1)
                    self.assertEqual(store.data.variants_by_chrom, {})
                    self.assertEqual(len(store.window_cache), 4)
                    self.assertTrue(all(
                        block.selection_method == "deferred_pair_specific_fixed_cohort"
                        for block in store.window_cache.values()
                    ))
                    self.assertTrue(any(len(block.variants) > 1
                                        for block in store.window_cache.values()))
            finally:
                logger.close()

    def test_bounded_cache_eviction_triggers_only_needed_window_query(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf_window_cache_size = 4
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(args.vcf.read_text())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with fake.patches():
                    source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                    store = workflow.BlockStore(source.data, args, logger, source)
                    store.prefetch_windows([("1", start, "right") for start in (0, 10, 20, 30)])
                    self.assertEqual(len(fake.queries), 1)
                    store.window("1", 40, "right")
                    store.window("1", 0, "right")
                    self.assertEqual(len(store.window_cache), 4)
                    self.assertEqual(fake.queries[1:], [[("1", 40, 50)], [("1", 0, 10)]])
            finally:
                logger.close()

    def test_large_window_batch_uses_one_query_without_flooding_lru(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf_window_cache_size = 4
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(args.vcf.read_text())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with fake.patches():
                    source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                    store = workflow.BlockStore(source.data, args, logger, source)
                    requested = [("1", start, "right") for start in range(0, 100, 10)]
                    blocks = store.load_window_batch(requested)
                    self.assertEqual(len(fake.queries), 1)
                    self.assertEqual(set(blocks), set(requested))
                    self.assertEqual(len(store.window_cache), 0)
                    with store._cache_lock:
                        store._active_window_batch = blocks
                    try:
                        for key in requested:
                            self.assertIs(store.window(*key), blocks[key])
                    finally:
                        with store._cache_lock:
                            store._active_window_batch = None
                    self.assertEqual(len(fake.queries), 1)
            finally:
                logger.close()

    def test_missing_index_binary_or_ambiguous_contigs_fail_without_scan_fallback(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf.write_text(self.boundary_content())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with mock.patch.object(indexed_vcf.shutil, "which", return_value=None), self.assertRaisesRegex(OSError, "not found"):
                    indexed_vcf.BCFtoolsVCFSource(args, logger)
                with FakeBCFtools(self.boundary_content()).patches(), self.assertRaisesRegex(OSError, "requires a BGZF"):
                    indexed_vcf.BCFtoolsVCFSource(args, logger)
                Path(str(args.vcf) + ".csi").touch()
                fake = FakeBCFtools(self.boundary_content(), "1\t1000\t1\nchr1\t1000\t1\n")
                with fake.patches(), self.assertRaisesRegex(OSError, "ambiguous"):
                    indexed_vcf.BCFtoolsVCFSource(args, logger)
                self.assertEqual(fake.queries, [])
            finally:
                logger.close()

    def test_query_failure_and_invalid_gt_are_fatal_not_qc_rejections(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf.write_text(self.boundary_content())
            Path(str(args.vcf) + ".csi").touch()
            logger = workflow.RunLogger(directory / "test.log")
            try:
                for malformed in (False, True):
                    fake = FakeBCFtools(self.boundary_content().replace("0|1", "0|2") if malformed else self.boundary_content())
                    with fake.patches():
                        source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                        fake.fail_next = not malformed
                        with self.assertRaises(indexed_vcf.VCFQueryError) as error:
                            source.fetch([Region("1", "1", 0, 20, "a", 0)])
                        self.assertNotIsInstance(error.exception, ValueError)
                        self.assertEqual(source.records_returned, 0)
            finally:
                logger.close()

    def test_absent_chromosome_is_logged_once_and_returns_empty_window_data(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            directory = Path(temp_dir)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "test")
            args.vcf.write_text(self.boundary_content())
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(self.boundary_content())
            logger = workflow.RunLogger(directory / "test.log")
            try:
                with fake.patches():
                    source = indexed_vcf.BCFtoolsVCFSource(args, logger)
                    for _ in range(2):
                        data = source.fetch([Region("2", "2", 0, 20, "absent", 0)])
                        self.assertEqual(data.variants_by_chrom, {})
                self.assertEqual(fake.queries, [])
                self.assertEqual((directory / "test.log").read_text().count("absent from the index"), 1)
            finally:
                logger.close()

    def test_background_region_query_error_stops_run_instead_of_rejecting_candidate(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "fatal")
            args.vcf_access = "bcftools"
            args.background_regions = args.background_pairs_per_distance = 2
            args.skip_background_region_pair_p = True
            args.skip_pair_test = True
            Path(str(args.vcf) + ".csi").touch()
            fake = FakeBCFtools(args.vcf.read_text())
            original = workflow.BlockStore.load_window_batch

            def load_window_batch(store, requested):
                fake.fail_next = True
                return original(store, requested)

            with fake.patches(), mock.patch.object(
                    workflow.BlockStore, "load_window_batch", load_window_batch), \
                    contextlib.redirect_stderr(io.StringIO()), self.assertRaises(indexed_vcf.VCFQueryError):
                workflow.run(args)
            log = (directory / "fatal.log").read_text()
            self.assertIn("RUN FAILED VCFQueryError", log)
            self.assertNotIn("REGION BACKGROUND TOTAL", log)
            self.assertFalse((directory / "fatal.run.json").exists())

    def test_full_workflow_indexed_and_preload_results_match_with_snp_range(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            progress_tests.BackgroundProgressTests.range_inputs(directory)
            Path(str(directory / "cohort.vcf") + ".csi").touch()
            fake = FakeBCFtools((directory / "cohort.vcf").read_text())
            for prefix, access in (("preload", "preload"), ("indexed", "bcftools")):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.vcf_access = access
                args.snps_per_block = None
                args.min_snps_per_block, args.max_snps_per_block = 1, 2
                args.max_attempt_factor = 100
                args.skip_region_entropy_matching = False
                args.region_entropy_max_wasserstein = 1.0
                with fake.patches(), contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv", "background_regions.tsv", "background_regions.bed"):
                self.assertEqual((directory / f"preload.{suffix}").read_text(),
                                 (directory / f"indexed.{suffix}").read_text(), suffix)
            metadata = json.loads((directory / "indexed.run.json").read_text())
            self.assertEqual(metadata["vcf_access"]["mode"], "bcftools")
            self.assertGreater(metadata["vcf_access"]["query_batches"], 0)
            self.assertIsNone(metadata["retained_variant_count"])
            self.assertIn("VCF QUERY DONE", (directory / "indexed.log").read_text())
            self.assertNotIn("VCF START preload", (directory / "indexed.log").read_text())

    def test_candidate_batching_preserves_results_and_reduces_index_queries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            Path(str(directory / "cohort.vcf") + ".csi").touch()
            query_counts = {}
            for prefix, batch_size in (("single", 1), ("batched", 1000)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.vcf_access = "bcftools"
                args.background_pairs_per_distance = 4
                args.background_regions = 3
                args.pair_batch_size = batch_size
                args.region_batch_size = batch_size
                args.skip_background_region_pair_p = True
                fake = FakeBCFtools(args.vcf.read_text())
                with fake.patches(), contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
                query_counts[prefix] = len(fake.queries)
            for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv",
                           "background_regions.tsv", "background_regions.bed"):
                self.assertEqual(
                    (directory / f"single.{suffix}").read_text(),
                    (directory / f"batched.{suffix}").read_text(),
                    suffix,
                )
            self.assertLess(query_counts["batched"], query_counts["single"])
            metadata = json.loads((directory / "batched.run.json").read_text())
            self.assertTrue(all(info["candidate_batch_size"] == 1000
                                for info in metadata["pair_sampling_counts"].values()))
            self.assertEqual(metadata["region_sampling_counts"]["target"]["candidate_batch_size"], 1000)


BCFTOOLS = os.environ.get("BLOCK_LD_TEST_BCFTOOLS") or shutil.which("bcftools")


@unittest.skipUnless(BCFTOOLS, "real bcftools CLI unavailable; set BLOCK_LD_TEST_BCFTOOLS to enable integration")
class RealBCFtoolsTests(unittest.TestCase):
    def test_real_soft_distance_and_unmatched_snp_counts_match_preload(self):
        from flanking_block_ld.tests import test_background_pair_soft_controls as soft_tests

        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            soft_tests.BackgroundPairSoftControlTests.write_inputs(directory)
            raw = directory / "cohort.vcf"
            raw.write_text(raw.read_text().replace(
                "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"))
            args = soft_tests.BackgroundPairSoftControlTests.args(directory, "soft_preload")
            args.background_pairs_per_distance = args.background_regions = 2
            args.region_length_tolerance_bp = 0
            args.vcf_prefiltered = True
            with contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)
            for fmt in ("csi", "tbi", "bcf"):
                path = directory / ("soft.bcf" if fmt == "bcf" else f"soft_{fmt}.vcf.gz")
                subprocess.run([BCFTOOLS, "view", "-Ob" if fmt == "bcf" else "-Oz",
                                "-o", str(path), str(raw)], check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocess.run([BCFTOOLS, "index", "--tbi" if fmt == "tbi" else "--csi", str(path)],
                               check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                args.output_prefix = directory / f"soft_{fmt}"
                args.vcf, args.vcf_access, args.bcftools = path, "bcftools", BCFTOOLS
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
                for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv",
                               "background_regions.tsv", "background_regions.bed"):
                    self.assertEqual((directory / f"soft_preload.{suffix}").read_text(),
                                     Path(f"{args.output_prefix}.{suffix}").read_text(), (fmt, suffix))

    def test_real_bgzf_csi_tbi_and_bcf_match_preload_with_phasing_and_prefiltered_qc(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            # Complete format/contig headers are needed when generating BCF.
            raw = directory / "cohort.vcf"
            content = raw.read_text().replace("##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n"
                                             "##contig=<ID=1,length=1000>\n"
                                             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            raw.write_text(content)
            for prefiltered in (False, True):
                prefix = f"baseline_{prefiltered}"
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.background_pairs_per_distance = args.background_regions = 2
                args.vcf_prefiltered = prefiltered
                args.skip_region_entropy_matching = False
                args.window_step_bp = 5
                args.region_length_tolerance_bp = 2
                args.skip_background_region_pair_p = True
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
                for fmt in ("csi", "tbi", "bcf"):
                    path = directory / (f"{fmt}_{prefiltered}.bcf" if fmt == "bcf" else f"{fmt}_{prefiltered}.vcf.gz")
                    subprocess.run([BCFTOOLS, "view", "-Ob" if fmt == "bcf" else "-Oz", "-o", str(path), str(raw)], check=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    subprocess.run([BCFTOOLS, "index", "--tbi" if fmt == "tbi" else "--csi", str(path)], check=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    args.output_prefix = directory / f"indexed_{fmt}_{prefiltered}"
                    args.vcf, args.vcf_access, args.bcftools = path, "bcftools", BCFTOOLS
                    with contextlib.redirect_stderr(io.StringIO()):
                        workflow.run(args)
                    for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv", "background_regions.tsv", "background_regions.bed"):
                        self.assertEqual((directory / f"{prefix}.{suffix}").read_text(),
                                         Path(f"{args.output_prefix}.{suffix}").read_text(), (fmt, prefiltered, suffix))


if __name__ == "__main__":
    unittest.main()
