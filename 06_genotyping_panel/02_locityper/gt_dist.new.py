#!/usr/bin/env python3

import argparse
from collections import defaultdict
import itertools
import sys
import operator
import os
import re
import numpy as np

import common


# Cigar operations, only small subset is allowed.
op_pattern = re.compile(r'[ID=X]')

def calc_directed_distance(line):
    cigar = None
    for entry in line[12:]:
        if entry.startswith('cg:'):
            cigar = entry[5:]
            break
    assert cigar

    # 1 = query, 2 = reference.
    # dist12 = number of basepairs, unique to ref.
    # dist21 = -//-, unique to query.
    dist12 = 0
    dist21 = 0
    n = len(cigar)
    i = 0
    for m in re.finditer(op_pattern, cigar):
        j = m.start()
        length = int(cigar[i:j])
        op = cigar[j]
        if op == 'X':
            dist12 += length
            dist21 += length
        elif op == 'D':
            dist12 += length
        elif op == 'I':
            dist21 += length
        i = j + 1
    assert i == n
    return dist12, dist21


class Distances:
    def __init__(self, discarded_path, paf_path, dir_dist=False):
        self.paf_path = paf_path
        self.discarded = {}
        if discarded_path is not None and os.path.exists(discarded_path):
            with common.open(discarded_path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    hap, haps2 = line.strip().split('=')
                    self.discarded[hap.strip()] = tuple(map(str.strip, haps2.split(',')))

        self.lengths = {}
        self.distances = defaultdict(dict)
        # dir_distances[hap1][hap2] = number of basepairs, unique to hap2 compared to hap1.
        self.dir_distances = defaultdict(dict) if dir_dist else None
        with common.open(paf_path) as paf:
            for line in paf:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                hap1 = line[0]
                hap2 = line[5]

                if hap1 not in self.lengths:
                    len1 = int(line[1])
                    for hap1a in self.group(hap1):
                        self.lengths[hap1a] = len1
                if hap2 not in self.lengths:
                    len2 = int(line[6])
                    for hap2a in self.group(hap2):
                        self.lengths[hap2a] = len2

                # Due to incorrect output PAF files in some version of `locityper align`, need to shift column indices.
                ix_shift = int(line[9] == '+')
                nmatches = int(line[9 + ix_shift])
                aln_size = int(line[10 + ix_shift])
                assert aln_size != 0, f'Missing alignment between {hap1} and {hap2}'
                assert nmatches <= aln_size
                dist = (aln_size - nmatches, aln_size)
                
                # Store all equivalent haplotype combinations in the distance dictionary.
                # Therefore, even if the PAF only contains APGp1_A versus APGp1_B,
                # entries such as HPRC_1 versus HPRC_2 are also created when the discarded file links them.
                for hap1a, hap2a in itertools.product(self.group(hap1), self.group(hap2)):
                    self.distances[hap1a][hap2a] = dist
                    self.distances[hap2a][hap1a] = dist

                if dir_dist:
                    dist12, dist21 = calc_directed_distance(line)
                    for hap1a, hap2a in itertools.product(self.group(hap1), self.group(hap2)):
                        self.dir_distances[hap1a][hap2a] = dist12
                        self.dir_distances[hap2a][hap1a] = dist21

        for hap, length in self.lengths.items():
            for hap1, hap2 in itertools.product(self.group(hap), repeat=2):
                self.distances[hap1][hap2] = (0, length)
                self.distances[hap2][hap1] = (0, length)
                if dir_dist:
                    self.dir_distances[hap1][hap2] = 0
                    self.dir_distances[hap2][hap1] = 0

        pattern = re.compile(r'[._][1-9]$')
        self.sample_haps = defaultdict(list)
        for hap in self.distances:
            m = re.search(pattern, hap)
            if m:
                self.sample_haps[hap[:m.start()]].append(hap)

    def group(self, hap):
        return (hap,) + self.discarded.get(hap, ())

    def group_size(self, hap):
        assert hap in self.lengths
        return len(self.group(hap))

    def get_sample_haplotypes(self, sample):
        return self.sample_haps.get(sample, ())

    def all_distances(self, genotype, allowed_samples=None):
        """
        Args:
            genotype: Target genotype (hap1, hap2, ...).
            allowed_samples: Optional set of sample names allowed as query haplotypes.
                             This is used to estimate availability for a smaller panel from the full PAF.
        """
        # Check whether all genotype haplotypes are present in the distance database.
        for hap in genotype:
            if hap not in self.distances:
                return None

        pred_dists = {}
        
        # Helper function: check whether a haplotype belongs to the allowed sample set.
        def is_allowed(hap_name):
            if allowed_samples is None:
                return True
            # Assume the naming format is Sample.Haplotype, for example HG002.1 or C003-CHA-E03.2.
            # Extract the sample name by removing the terminal haplotype suffix.
            sample_name = hap_name.rsplit('.', 1)[0]
            return sample_name in allowed_samples

        # ================== Haploid genotype handling ==================
        if len(genotype) == 1:
            hap1 = genotype[0]
            # Iterate over distances between this haplotype and all other haplotypes.
            for hap2, (edit, size) in self.distances[hap1].items():
                if not is_allowed(hap2):
                    continue
                if size == 0:
                    continue
                div = edit / size
                query = (hap2,)  # The query genotype is also haploid.
                pred_dists[query] = (div, edit, size)

        # ================== Diploid genotype handling ==================
        elif len(genotype) == 2:
            hap_dists = [self.distances[genotype[0]].items(), self.distances[genotype[1]].items()]
            for (h1, (e1, s1)), (h2, (e2, s2)) in itertools.product(*hap_dists):
                # Filter out query genotypes unless both haplotypes are allowed.
                if not (is_allowed(h1) and is_allowed(h2)):
                    continue
                
                edit = e1 + e2
                size = s1 + s2
                if size == 0:
                    continue
                div = edit / size
                query = (h1, h2) if h1 <= h2 else (h2, h1)
                last_pred = pred_dists.get(query)
                if last_pred is None or last_pred[0] > div:
                    pred_dists[query] = (div, edit, size)
        
        else:
            raise ValueError(f"Unsupported ploidy: {len(genotype)}. Genotype must be haploid or diploid.")
            
        return pred_dists

    def calc_distance(self, gt1, gt2):
        assert len(gt1) == len(gt2)
        best_div = np.inf
        best_distances = None

        for perm2 in itertools.permutations(gt2):
            distances = []
            sum_edit = 0
            sum_size = 0
            for hap1, hap2 in zip(gt1, perm2):
                if hap1 is None:
                    distances.append((None, None))
                    continue
                try:
                    edit, size = self.distances[hap1][hap2]
                except KeyError:
                    sys.stderr.write(f'Cannot calculate distance between {",".join(gt1)} and {",".join(gt2)}'
                        f' (missing distance {hap1} - {hap2}) (see {self.paf_path})\n')
                    distances.append((None, None))
                    continue
                sum_edit += edit
                sum_size += size
                distances.append((edit, size))

            div = sum_edit / sum_size if sum_size else np.inf
            if div <= best_div:
                best_div = div
                best_distances = distances
        return GtDist(best_distances)

    def find_closest(self, gt, loo=True, excl_haps=()):
        """
        Find closest genotype to `gt`.
        loo: bool - this is leave-one-out experiment, do not include haplotypes from the same sample.
        excl_haps: set of excluded haplotype names.
        """
        loo_gt = []
        distances = []
        for hap in gt:
            if hap is None:
                loo_gt.append(None)
                distances.append((None, None))
                continue

            best_hap = None
            best_div = np.inf
            best_edit = None
            for hap2, (edit, size) in self.distances[hap].items():
                if (loo and hap2 in gt) or hap2 in excl_haps:
                    continue
                if edit / size < best_div:
                    best_div = edit / size
                    best_edit = (edit, size)
                    best_hap = hap2
            loo_gt.append(best_hap)
            distances.append(best_edit)
        return loo_gt, GtDist(distances)

    def average_divergence(self):
        edits = []
        divs = []
        for hap1, hap_dists in self.distances.items():
            for hap2, (edit, size) in hap_dists.items():
                if hap1 < hap2:
                    edits.append(edit)
                    divs.append(edit / size)
        return np.mean(edits), np.mean(divs)


def edit_to_str(edit, size):
    div = edit / size
    qv = np.inf if div == 0 else -10 * np.log10(div)
    return f'{edit}\t{size}\t{div:.9f}\t{qv:.9f}'


class GtDist:
    def __init__(self, distances):
        self.distances = distances

    def iter_strs(self):
        sum_edit = 0
        sum_size = 0
        all_present = True
        for edit, size in self.distances:
            if edit is None:
                all_present = False
                yield 'NA\tNA\tNA\tNA'
            else:
                sum_edit += edit
                sum_size += size
                yield edit_to_str(edit, size)

        if all_present:
            yield edit_to_str(sum_edit, sum_size)
        else:
            yield 'NA\tNA\tNA\tNA'


def get_genotype(s, split, sep):
    tup = tuple(map(str.strip, s.strip().split(split)))
    assert tup
    gt_str = ','.join(tup)

    if len(tup) == 1:
        # A single input value is treated as a haploid genotype.
        # Users can provide "sample.1,sample.2" to specify a diploid genotype.
        return gt_str, (tup[0],)
        
    if len(tup) == 2:
        # Two comma-separated values are treated as a diploid genotype.
        return gt_str, tup

    # Other ploidies are not supported.
    raise ValueError(f'Genotype must be haploid (e.g., "hap1") or diploid (e.g., "hap1,hap2"), but got "{s}"')


def load_target_genotypes(args):
    genotypes = []
    if args.genotype:
        for val in args.genotype:
            genotypes.append(get_genotype(val, ',', args.sep))
    if args.gt_file:
        with common.open(args.gt_file) as f:
            for line in f:
                genotypes.append(get_genotype(line, ',', args.sep))
    return genotypes

# Load the list of samples allowed as query haplotypes.
def load_allowed_query_samples(path):
    if not path:
        return None
    allowed = set()
    with common.open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Assume each line contains one sample ID, for example HG002.
                allowed.add(line.split()[0])
    return allowed

def is_loo(target, query):
    return all(h1 != h2 for h1, h2 in itertools.product(target, query))


def calc_gt_distances(genotypes, distances, out, max_entries, loo=False, allowed_samples=None):
    for gt_str, genotype in genotypes:
        # Pass allowed_samples to all_distances for query-panel filtering.
        pred_dists = distances.all_distances(genotype, allowed_samples=allowed_samples)
        
        if pred_dists is None:
            out.write(f'{gt_str}\t*\tNA\tNA\tNA\tNA\tNA\n')
            continue
        pred_dists = sorted(pred_dists.items(), key=operator.itemgetter(1))

        num_entries = 0
        for query, (div, edit, size) in pred_dists:
            loo_status = 'T' if is_loo(genotype, query) else 'F'
            if loo and loo_status == 'F':
                continue
            if num_entries >= max_entries:
                break
            qv = np.inf if div == 0 else -10 * np.log10(div)
            out.write(f'{gt_str}\t{",".join(query)}\t{loo_status}\t{edit}\t{size}\t{div:.9f}\t{qv:.4f}\n')
            num_entries += 1


def main():
    parser = argparse.ArgumentParser(
        description='Calculating distances between a target genotype and all other genotypes',
        usage='%(prog)s -i alns.paf -d discarded.txt (-g hap1,hap2 | -G genotypes.txt) -o out.csv')
    parser.add_argument('-i', '--input', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-d', '--discarded', metavar='FILE',
        help='File with discarded haplotypes.')
    parser.add_argument('-g', '--genotype', metavar='STR', nargs='+',
        help='One or more target genotype (haplotypes through comma) or sample name.')
    parser.add_argument('-G', '--gt-file', metavar='FILE',
        help='List of target genotypes/samples.')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output CSV file.')
    parser.add_argument('-s', '--sep', metavar='STR', default='.',
        help='Separator between sample and haplotype [default: %(default)s].')
    parser.add_argument('-n', '--max-entries', metavar='INT', type=int,
        help='Output at most INT entries per target genotype [default: all].')
    parser.add_argument('--loo', action='store_true',
        help='Experiment performed in the leave-one-out setting. [default: False]', default=False)
    # Optional query-sample whitelist.
    parser.add_argument('-Q', '--query-samples', metavar='FILE',
        help='File containing list of sample names allowed to be used as queries. '
             'If provided, only haplotypes belonging to these samples will be considered as matches.')
    
    args = parser.parse_args()

    genotypes = load_target_genotypes(args)
    distances = Distances(args.discarded, args.input)
    allowed_samples = load_allowed_query_samples(args.query_samples)

    max_entries = args.max_entries or sys.maxsize
    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('target\tquery\tloo\tedit_dist\taln_size\tdivergence\tqv\n')
        calc_gt_distances(genotypes, distances, out, max_entries, args.loo, allowed_samples)


if __name__ == '__main__':
    main()
    
