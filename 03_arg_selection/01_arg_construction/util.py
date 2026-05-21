# Common utilities for processing long-read/short-read data and associated IGD files,
# ARG tree-sequences, etc.
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from typing import Any, Tuple, Optional, Union, List, TextIO
import json
import os
import subprocess
import sys

THISDIR = os.path.dirname(os.path.realpath(__file__))


# You can optionally pass in a coordinate map that maps the ARG file's positions
# to another coordinate system.
def parse_map(filename: str) -> Tuple[str, Any]:
    parts = filename.split(":")
    mapped = None
    if len(parts) > 1:
        assert len(parts) == 2, f"Unexpected ':' in filename: {filename}"
        with open(parts[1]) as f:
            mapped = {int(k): int(v) for k, v in json.load(f).items()}
        filename = parts[0]
        print(f"Using coordinate map {parts[1]} for {filename}", file=sys.stderr)

    def mapit(x):
        if mapped is None:
            return x
        if x in mapped:
            return mapped[x]
        return None

    return filename, mapit


@dataclass_json
@dataclass
class ExperimentConfig:
    window_size: int  # In base-pairs
    grid_size: int  # For grid-based plots/analysis of paired coalescences.
    ne: float  # For any scaling that uses effective population size.
    mut_rate: float
    fake_mutmap: str
    chain_to_chm13: str  # Chain file from GRCH38 -> CHM13
    chain_to_grch38: str  # Chain file from CHM13 -> GRCH38
    max_missingness: float  # Maximum proportion of missing alleles allowed per sample
    grch38_ancestral: str
    chm13_ancestral: str
    grch38_ratemaps: str
    chm13_ratemaps: str
    mcmc_samples: int
    mcmc_thin: int

    # These optional arguments are loaded from the environment, not the config file.
    data_dir: Optional[str] = None
    output_dir: Optional[str] = None


def load_config() -> ExperimentConfig:
    with open(os.path.join(THISDIR, "config", "config.json")) as f:
        result = ExperimentConfig.from_dict(json.load(f))
    assert "DATA_DIR" in os.environ, "Please supply the DATA_DIR environment variable"
    assert (
        "OUTPUT_DIR" in os.environ
    ), "Please supply the OUTPUT_DIR environment variable"
    result.data_dir = os.environ["DATA_DIR"]
    if not os.path.exists(result.data_dir):
        os.mkdir(result.data_dir)
    result.output_dir = os.environ["OUTPUT_DIR"]
    if not os.path.exists(result.output_dir):
        os.mkdir(result.output_dir)

    def resolve(path: str) -> str:
        return path.format(
            **{
                "DATA_DIR": result.data_dir,
                "CONFIG_DIR": os.path.join(THISDIR, "config"),
            }
        )

    # Resolve paths that might contain format string variables.
    result.chain_to_chm13 = resolve(result.chain_to_chm13)
    result.chain_to_grch38 = resolve(result.chain_to_grch38)
    result.grch38_ancestral = resolve(result.grch38_ancestral)
    result.chm13_ancestral = resolve(result.chm13_ancestral)
    result.grch38_ratemaps = resolve(result.grch38_ratemaps)
    result.chm13_ratemaps = resolve(result.chm13_ratemaps)
    result.fake_mutmap = resolve(result.fake_mutmap)
    return result


def which(exe: str, required=False) -> Optional[str]:
    """
    Find the named executable, first via system PATH and then via the Python PATH.

        :param exe: The executable name.
    :param required: If True, throw an exception when not found instead of returning None.
    :return: None if the executable is not found.
    """
    try:
        result = (
            subprocess.check_output(["which", exe], stderr=subprocess.STDOUT)
            .decode("utf-8")
            .strip()
        )
    except subprocess.CalledProcessError:
        result = None
    if result is None:
        for p in sys.path + [os.path.realpath(os.path.dirname(__file__))]:
            p = os.path.join(p, exe)
            if os.path.isfile(p):
                result = p
                break
    if required and result is None:
        raise RuntimeError(f"Could not find executable {exe}")
    return result


def run(cmd: Union[str, List[str]], shell: bool = False, verbose: bool = False):
    if verbose:
        print(f"Running: {cmd}")
    if shell:
        subprocess.check_call(cmd, shell=True)
    else:
        subprocess.check_call([str(c) for c in cmd])


def remove_ext(filename: str, ext: Optional[str] = None) -> str:
    file_ext = filename.split(".")[-1]
    removed = ".".join(filename.split(".")[:-1])
    assert len(file_ext) < len(filename), "Filename has no extension"
    if ext is not None:
        assert ext == file_ext, f"Unexpected file extension on {filename}"
    return removed


def traverse_vcf(file_obj, callback, print_meta=False, context: Any = None) -> bool:
    traversed_entire = True
    for line in file_obj:
        line = line.strip()
        if line.startswith("##"):
            if print_meta:
                print(line, file=print_meta)
            continue
        elif line.startswith("#"):
            header = line.split("\t")
            individuals = header[9:]
        else:
            data = line.split("\t")
            chrom, position, var_id, ref, alt, qual, filt, info, fmt = data[0:9]
            genotype = data[9:]
            if not callback(
                individuals, int(position), ref, alt, genotype, context=context
            ):
                traversed_entire = False
                break
    return traversed_entire


def filter_vcf(file_obj: TextIO, callback, out_file: TextIO):
    for line in file_obj:
        line = line.strip()
        if line.startswith("##"):
            print(line, file=out_file)
            continue
        elif line.startswith("#"):
            header = line.split("\t")
            individuals = header[9:]
            print(line, file=out_file)
        else:
            data = line.split("\t")
            chrom, position, var_id, ref, alt, qual, filt, info, fmt = data[0:9]
            genotype = data[9:]
            if callback(individuals, int(position), ref, alt, genotype):
                print(line, file=out_file)
