#!/usr/bin/env python
"""
Attempts to identify groups of clones in a dataset. The script
(1) conducts pairwise comparisons (genetic similarity) for all individuals in
a `.vcf` file,
(2) produces a histogram of genetic similarities,
(3) lists the highest matches to assess for a potential clonal threshold,
(4) clusters the groups of clones based on a particular threshold (supplied or
roughly inferred),
(5) lists the clonal individuals that can be removed from the dataset
(so that one individual with the least amount of missing data remains), and
(6) produces a multi-panel PDF report (NJ tree with clone groups highlighted +
per-sample % genotyped bars + within/between-population similarity histograms).

If optional popfile is given, then clonal groups are sorted by population and
the histogram contrasts within- vs between-population comparisons.

Note: Firstly, the script is run with a `.vcf` file and an optional popfile
to produce an output file (e.g. `python3 vcf_clone_detect.py --vcf
vcf_file.vcf --pop pop_file.txt --output compare_file.csv`). Secondly, it can
be rerun using the precalculated similarities under different thresholds
(e.g. `python3 vcf_clone_detect.py --input compare_file.csv --threshold 94.5`).

Several similarity measures are available via `--method` (default `ibs`, which
reproduces the original allelic-similarity behaviour exactly). `het-masked` and
`dosage` give a sharper clone/non-clone boundary; `single-read` (AD-weighted
when allelic depths are present) is robust to genotype-call artifacts and gives
cleaner deep topology, but is NOT recommended for clone thresholds.

Similarity measures and their PLINK / ANGSD correspondence
----------------------------------------------------------
All measures are reported as a 0-100 % similarity over sites genotyped in both
individuals (pairwise-complete); the genetic distance used for the tree is
`1 - similarity/100`. `da`, `db` are alt-allele dosages (0, 1, 2).

  ibs          Allele-sharing IBS: per-site score `1 - |da-db|/2` (het-vs-het
               counts as a full match). This is the L1 / allele-count family.
               Equivalent to PLINK `--distance 1-ibs flat-missing` (verified
               identical to rounding), and to PLINK's default
               `--distance` / `allele-ct` output (same metric, rescaled to
               `n_variants * mean|da-db|`).

  dosage       Squared-Euclidean (L2) distance on allele dosage: per-site
               `1 - (da-db)^2 / max_sq`. Corresponds to PLINK `--make-rel cov`
               (centred covariance): `sum (da-db)^2 == M * (Cjj + Ckk - 2*Cjk)`
               (verified identical to rounding). NOT the same as the default,
               allele-frequency-standardised GRM (`--make-rel`), and NOT a
               PLINK `--distance` flavour.

  single-read  Emulates ANGSD `-doIBS 1 -makeMatrix 1` (single-read / random
               haploid sampling). Expected per-site mismatch `p1 + p2 - 2*p1*p2`
               where p is the alt fraction (AD-weighted alt/(ref+alt) when the
               AD field is present, else the called-genotype dosage/2). Differs
               from `ibs` only in that het-vs-het has an expected distance of
               0.5 rather than 0. NB: true ANGSD single-read needs read-level
               data (BAM/genotype likelihoods); this approximates it from a VCF.

  het-masked   Homozygous-only IBS (heterozygous sites ignored). No direct
               PLINK/ANGSD equivalent.

Correspondences hold on matching SNP sets: this script uses pairwise-complete
sites, whereas PLINK/ANGSD apply their own missing-data scaling, so values can
drift slightly on data with missing genotypes (the underlying metric is the
same).
"""
import sys
import os
import argparse
import operator
import itertools
import math
import numpy as np

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
FIRST_GENOTYPE_COLUMN = 9
IND_STR_LEN = 50
POP_STR_LEN = 50
MATCHES_ADDITIONAL_ROWS = 5
HIST_RANGE = range(1, 101)
DEF_THRESHOLD = 85.0

METHODS = ('ibs', 'het-masked', 'dosage', 'single-read')
DEFAULT_METHOD = 'ibs'

C_IND1 = 'ind1'
C_IND2 = 'ind2'
C_IND1_SNPS = 'ind1_snps'
C_IND2_SNPS = 'ind2_snps'
C_BOTH_SNPS = 'both_snps'
C_MATCH = 'match'
C_MATCH_PERC = 'match_perc'
C_POP = 'pop'

COMPARISONS_DTYPES = [(C_IND1, np.str_, IND_STR_LEN),
                      (C_IND2, np.str_, IND_STR_LEN),
                      (C_IND1_SNPS, np.int64),
                      (C_IND2_SNPS, np.int64),
                      (C_BOTH_SNPS, np.int64),
                      (C_MATCH, np.float64),
                      (C_MATCH_PERC, np.float64),
                      (C_POP, np.str_, POP_STR_LEN)]

OUTPUT_FILE_DELIM = ','
OUTPUT_FILE_HEADER = OUTPUT_FILE_DELIM.join([C_IND1, C_IND2, C_IND1_SNPS,
                                             C_IND2_SNPS, C_BOTH_SNPS,
                                             C_MATCH, C_MATCH_PERC, C_POP])
OUTPUT_FILE_FORMAT = OUTPUT_FILE_DELIM.join(['%s', '%s', '%i', '%i', '%i',
                                             '%f', '%f', '%s'])
NLOCI_COMMENT = 'n_loci='


class CloneGroup(object):

    def __init__(self, row):
        self.indivs = set([str(row[C_IND1]), str(row[C_IND2])])
        self.pops = set(str(row[C_POP]).split('-'))
        self.min_sim_score = self.max_sim_score = float(row[C_MATCH_PERC])
        if int(row[C_IND1_SNPS]) >= int(row[C_IND2_SNPS]):
            self.best_indiv = str(row[C_IND1])
            self.best_indiv_snps = int(row[C_IND1_SNPS])
        else:
            self.best_indiv = str(row[C_IND2])
            self.best_indiv_snps = int(row[C_IND2_SNPS])

    def add_clone_from_row(self, row):
        self.indivs.update([str(row[C_IND1]), str(row[C_IND2])])
        self.pops.update(str(row[C_POP]).split('-'))
        if float(row[C_MATCH_PERC]) < self.min_sim_score:
            self.min_sim_score = float(row[C_MATCH_PERC])
        if row[C_MATCH_PERC] > self.max_sim_score:
            self.max_sim_score = float(row[C_MATCH_PERC])
        if row[C_IND1_SNPS] >= self.best_indiv_snps:
            self.best_indiv = str(row[C_IND1])
            self.best_indiv_snps = int(row[C_IND1_SNPS])
        if row[C_IND2_SNPS] >= self.best_indiv_snps:
            self.best_indiv = str(row[C_IND2])
            self.best_indiv_snps = int(row[C_IND2_SNPS])

    def get_formatted_clone_info(self):
        if self.min_sim_score == self.max_sim_score:
            score_range = '{0} %'.format(self.min_sim_score)
        else:
            score_range = '{0}-{1} %'.format(self.min_sim_score,
                                             self.max_sim_score)
        info = '{0}: {1} ({2})'.format('-'.join(sorted(self.pops)),
                                       ', '.join(sorted(self.indivs)),
                                       score_range)
        return info

    def merge_group(self, other, row):
        """ Merge another clone group into this one (single-linkage), also
        folding in the linking row's score/pop. Only reached for non-transitive
        clusters, which the original implementation rejected outright. """
        self.indivs.update(other.indivs)
        self.pops.update(other.pops)
        self.pops.update(str(row[C_POP]).split('-'))
        mp = float(row[C_MATCH_PERC])
        self.min_sim_score = min(self.min_sim_score, other.min_sim_score, mp)
        self.max_sim_score = max(self.max_sim_score, other.max_sim_score, mp)
        if other.best_indiv_snps > self.best_indiv_snps:
            self.best_indiv = other.best_indiv
            self.best_indiv_snps = other.best_indiv_snps

    def get_samples_to_remove(self):
        return sorted(self.indivs - set([self.best_indiv]))


def get_snp_match(genotype1, genotype2):
    """ Get match value for two genotypes (one SNP) [reference implementation] """
    if genotype1 == genotype2:
        match_score = 1
    elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2] or \
            genotype1[0] == genotype2[2] or genotype1[2] == genotype2[0]:
        match_score = 0.5
    else:
        match_score = 0
    return match_score


def get_pop_assignments_from_popfile(pop_filename):
    """ Initialise dict of pops with lists of indvs from popfile """
    indivs_pops = {}
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        cols = line.rstrip().replace(',', ' ').split()
        if not cols:
            continue
        indiv = cols[0]
        pop = cols[1]
        indivs_pops[indiv] = pop
    return indivs_pops


def get_pop_group(individual1, individual2, indivs_pops):
    """ Define whether a comparison is within or between pops (as in original) """
    if individual1 in indivs_pops and individual2 in indivs_pops:
        if indivs_pops[individual1] == indivs_pops[individual2]:
            return indivs_pops[individual1]
        pops = sorted([indivs_pops[individual1], indivs_pops[individual2]])
        return '-'.join(pops)
    return 'NA'


def load_vcf_arrays(vcf_filename, need_ad=False):
    """ Read a VCF once into vectorised arrays.

    Returns a dict with:
      names       : list of sample names (VCF column order)
      n_loci      : total number of SNP records
      state_code  : sites x samples int array; small code per ordered genotype
                    (e.g. '0/1' vs '1/0' are distinct codes), -1 for missing
      code_match  : (n_states x n_states) matrix of get_snp_match() scores
      dosage      : sites x samples int8 (allele1+allele2), -1 for missing
      max_dosage  : maximum observed alt-allele count (for dosage normalisation)
      multiallelic: True if any site has a multi-allele ALT
      ad_p        : sites x samples float alt-read fraction (or None)
      ad_present  : sites x samples float {0,1} where reads available (or None)
    """
    states = {}                  # genotype-string (gt[:3]) -> code
    state_tokens = []            # code -> gt[:3] string
    rows_state = []
    rows_dos = []
    rows_p = []
    rows_r = []
    names = None
    ad_seen = False
    multiallelic = False
    max_dosage = 1

    with open(vcf_filename, 'r') as vcf_file:
        for line in vcf_file:
            if line[:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
                names = line.rstrip().split('\t')[FIRST_GENOTYPE_COLUMN:]
                continue
            if line[0] == HEADER_CHAR:
                continue
            cols = line.rstrip().split('\t')
            if len(cols) <= FIRST_GENOTYPE_COLUMN:
                continue
            if ',' in cols[4]:
                multiallelic = True
            fmt = cols[8].split(':')
            ad_idx = fmt.index('AD') if (need_ad and 'AD' in fmt) else -1
            srow = []
            drow = []
            prow = []
            rrow = []
            for col in cols[FIRST_GENOTYPE_COLUMN:]:
                gt = col[:3]
                if gt[0] == '.':
                    srow.append(-1)
                    drow.append(-1)
                else:
                    code = states.get(gt)
                    if code is None:
                        code = len(state_tokens)
                        states[gt] = code
                        state_tokens.append(gt)
                    srow.append(code)
                    a0 = int(gt[0])
                    a2 = int(gt[2])
                    drow.append(a0 + a2)
                    if a0 > max_dosage:
                        max_dosage = a0
                    if a2 > max_dosage:
                        max_dosage = a2
                if ad_idx >= 0:
                    parts = col.split(':')
                    if len(parts) > ad_idx and parts[ad_idx] not in ('.', './.', ''):
                        ad = parts[ad_idx].split(',')
                        try:
                            ref = int(ad[0])
                            alt = int(ad[1])
                        except (ValueError, IndexError):
                            ref = alt = 0
                        tot = ref + alt
                        if tot > 0:
                            prow.append(alt / tot)
                            rrow.append(1.0)
                            ad_seen = True
                        else:
                            prow.append(0.0)
                            rrow.append(0.0)
                    else:
                        prow.append(0.0)
                        rrow.append(0.0)
                else:
                    prow.append(0.0)
                    rrow.append(0.0)
            rows_state.append(srow)
            rows_dos.append(drow)
            rows_p.append(prow)
            rows_r.append(rrow)

    if names is None:
        sys.exit('Error: no #CHROM header line found in VCF')

    state_code = np.array(rows_state, dtype=np.int16)
    dosage = np.array(rows_dos, dtype=np.int8)
    n_states = len(state_tokens)
    # Build the get_snp_match() score matrix between ordered genotype states
    code_match = np.zeros((n_states, n_states), dtype=np.float64)
    for a in range(n_states):
        for b in range(n_states):
            code_match[a, b] = get_snp_match(state_tokens[a], state_tokens[b])
    # Homozygous iff the two alleles are identical (exact, incl. multiallelic)
    state_is_hom = np.array([tok[0] == tok[2] for tok in state_tokens],
                            dtype=bool)

    result = {
        'names': names,
        'n_loci': state_code.shape[0],
        'state_code': state_code,
        'code_match': code_match,
        'n_states': n_states,
        'state_is_hom': state_is_hom,
        'dosage': dosage,
        'max_dosage': max(max_dosage, 1),
        'multiallelic': multiallelic,
        'ad_p': None,
        'ad_present': None,
    }
    if need_ad and ad_seen:
        result['ad_p'] = np.array(rows_p, dtype=np.float64)
        result['ad_present'] = np.array(rows_r, dtype=np.float64)
    return result


def _match_denom_matrices(data, method):
    """ Compute pairwise (match, denom) matrices for a method.

    Returns (match_mat, denom_mat) as samples x samples numpy arrays (diagonal
    is meaningless). The similarity percentage is then computed per-pair as
    round((match/denom)*100, 2) (matching the original operation order), so the
    `--threshold` semantics carry across all methods. For `dosage`/`single-read`
    `match` is set to similarity*denom so this identity still holds.
    """
    n = len(data['names'])
    dosage = data['dosage']
    present = (data['state_code'] >= 0).astype(np.float64)   # sites x samples
    both = present.T @ present                               # both-genotyped count

    if method == 'ibs':
        # Exact replica of original get_snp_match() via ordered-state one-hots
        state_code = data['state_code']
        W = data['code_match']
        ns = data['n_states']
        onehot = [((state_code == s).astype(np.float64)) for s in range(ns)]
        match = np.zeros((n, n), dtype=np.float64)
        for a in range(ns):
            Xa = onehot[a]
            for b in range(ns):
                w = W[a, b]
                if w != 0.0:
                    match += w * (Xa.T @ onehot[b])
        denom = both

    elif method == 'het-masked':
        # Only sites homozygous in both individuals; match if same hom state
        state_code = data['state_code']
        ns = data['n_states']
        hom_codes = [s for s in range(ns) if data['state_is_hom'][s]]
        hom_onehot = [((state_code == s).astype(np.float64)) for s in hom_codes]
        hom_present = np.zeros_like(present)
        for X in hom_onehot:
            hom_present += X
        denom = hom_present.T @ hom_present
        match = np.zeros((n, n), dtype=np.float64)
        for X in hom_onehot:
            match += X.T @ X            # same homozygous state in both

    elif method == 'dosage':
        # similarity = 100 * (1 - mean((da-db)^2) / max_sq)
        max_dose = data['max_dosage'] + data['max_dosage']  # max allele-count value
        max_sq = float(max_dose * max_dose) if max_dose > 0 else 1.0
        dvals = sorted(set(int(v) for v in np.unique(dosage) if v >= 0))
        onehot = {v: ((dosage == v).astype(np.float64)) for v in dvals}
        sumsq = np.zeros((n, n), dtype=np.float64)
        for u in dvals:
            for v in dvals:
                if u != v:
                    sumsq += ((u - v) ** 2) * (onehot[u].T @ onehot[v])
        denom = both
        with np.errstate(divide='ignore', invalid='ignore'):
            mean_sq = np.where(denom > 0, sumsq / denom, 0.0)
        sim = 1.0 - (mean_sq / max_sq)
        match = sim * denom            # so match/denom == sim
        return match, denom

    elif method == 'single-read':
        if data['ad_p'] is not None:
            p = data['ad_p']
            r = data['ad_present']
            denom = r.T @ r
            pp = p * r
            dist = (pp.T @ r + r.T @ pp - 2.0 * (pp.T @ pp))
        else:
            sys.stderr.write('Note: no AD field found; single-read falls back '
                             'to genotype-based (haploidised) sampling.\n')
            f = np.where(dosage >= 0, dosage / 2.0, 0.0)
            fp = f * present
            denom = both
            dist = (fp.T @ present + present.T @ fp - 2.0 * (fp.T @ fp))
        with np.errstate(divide='ignore', invalid='ignore'):
            dist = np.where(denom > 0, dist / denom, 0.0)
        sim = 1.0 - dist
        match = sim * denom
        return match, denom

    else:
        sys.exit('Error: unknown method `{0}`'.format(method))

    return match, denom


def build_comparisons(data, indivs_pops, method):
    """ Build the structured comparisons array (same schema as original) using
    vectorised, method-aware similarity computation. """
    names = data['names']
    n = len(names)
    match_mat, both_mat = _match_denom_matrices(data, method)
    # per-sample genotyped counts (== ind*_snps in the original)
    geno_count = (data['state_code'] >= 0).sum(axis=0).astype(np.int64)

    unique_pairs = int((math.pow(n, 2) - n) / 2)
    comparisons = np.zeros(unique_pairs, dtype=COMPARISONS_DTYPES)
    index = 0
    for i, j in itertools.combinations(range(n), 2):
        ind1, ind2 = names[i], names[j]
        match = float(match_mat[i, j])
        both = int(round(both_mat[i, j]))
        # Python round() and original (match/both)*100 order, for exact parity
        match_perc = round((match / both) * 100, 2) if both > 0 else 0.0
        comparisons[index] = (str(ind1), str(ind2),
                              int(geno_count[i]), int(geno_count[j]),
                              both, match, match_perc,
                              str(get_pop_group(ind1, ind2, indivs_pops)))
        index += 1
    comparisons[::-1].sort(order=C_MATCH_PERC)
    print('{0} comparisons completed'.format(comparisons.size))
    return comparisons


def get_pairwise_comparisons_from_input_file(input_filename):
    """ Load pairwise comparisons between all individuals from input_file """
    comparisons = np.genfromtxt(input_filename, dtype=COMPARISONS_DTYPES,
                                delimiter=',')
    comparisons[::-1].sort(order=C_MATCH_PERC)
    print('{0} comparisons loaded from `{1}`'.format(comparisons.size,
                                                     input_filename))
    return comparisons


def read_nloci_from_input_file(input_filename):
    """ Read the optional `# n_loci=N` comment from a comparison file """
    with open(input_filename, 'r') as f:
        for line in f:
            if not line.startswith(HEADER_CHAR):
                break
            if NLOCI_COMMENT in line:
                try:
                    return int(line.split(NLOCI_COMMENT)[1].split()[0])
                except (ValueError, IndexError):
                    pass
    return None


def save_pairwise_comparisons_to_input_file(output_filename, comparisons,
                                            n_loci=None):
    """ Save pairwise comparisons to output_file (with optional n_loci) """
    header = OUTPUT_FILE_HEADER
    if n_loci is not None:
        header = '{0} {1}\n{2}'.format(NLOCI_COMMENT, n_loci, OUTPUT_FILE_HEADER)
    np.savetxt(output_filename, comparisons, fmt=OUTPUT_FILE_FORMAT,
               header=header)
    print('Comparisons outputted to file: `{0}`'.format(output_filename))


def output_ascii_hist(raw_values, bin_values):
    """ Plot a text-based histogram """
    values, bins = np.histogram(raw_values, bins=bin_values)
    output_lines = []
    graph_multiplier = 1
    lower_bound_flag = False
    previous_value = breakpoint = display_lines = 0
    for index, value in enumerate(values):
        if value > 0:
            if not lower_bound_flag:
                lower_bound_flag = True
                lower_bound = bins[index]
            else:
                upper_bound = bins[index]
        if lower_bound_flag:
            graph_bar = '*' * int(value * graph_multiplier)
            if len(graph_bar) > 70:
                graph_bar = graph_bar[:69] + '#'

            output_lines.append(
                '{:3d} {:7d} {:70s}'.format(bins[index], value,
                                            graph_bar[:70]))
    print('\n'.join(output_lines[:(upper_bound - lower_bound + 2)]))


def output_highest_matches(comparisons, threshold):
    """ Output list of highest matches """
    extra_rows = last_value = diff = highest_diff = 0
    highest_diff_max_perc = highest_diff_min_perc = 0
    output_lines = []
    for row in np.nditer(comparisons):
        # Only output matches above threshold (+ several additional rows)
        if threshold == 0 and row[C_MATCH_PERC] < DEF_THRESHOLD:
            break
        if threshold > 0 and row[C_MATCH_PERC] < threshold:
            extra_rows += 1
            if extra_rows > MATCHES_ADDITIONAL_ROWS:
                break
        # Keep track of largest difference between sequential matches
        if last_value != 0:
            diff = round(last_value - row[C_MATCH_PERC], 2)
            if diff > highest_diff:
                highest_diff = diff
                highest_diff_min_perc = row[C_MATCH_PERC]
                highest_diff_max_perc = last_value

        output_lines.append(('{0}\t{1}\t[{2}]\t{3} vs {4}'
                             '\t{5}/{6}\t{7}\t{8}').format(row[C_MATCH_PERC],
                                                           diff,
                                                           row[C_POP],
                                                           row[C_IND1], row[
                                                               C_IND2],
                                                           row[C_MATCH],
                                                           row[C_BOTH_SNPS],
                                                           row[C_IND1_SNPS],
                                                           row[C_IND2_SNPS]))
        last_value = row[C_MATCH_PERC]

    # Determine threshold value
    if threshold > 0:
        threshold_msg = 'Manual threshold'
    elif threshold == 0:
        threshold_msg = 'Potential threshold'
        if int(highest_diff_max_perc) > highest_diff_min_perc:
            threshold = float(int(highest_diff_max_perc))
        else:
            threshold = (highest_diff_max_perc - highest_diff_min_perc) / 2

    # Output list of matches with a break at the threshold

    for line in output_lines:
        if float(line.split()[0]) < threshold and threshold_msg:
            print('{0}\t{1} {2} {1}'.format(round(threshold, 2), '-' * 20,
                                            threshold_msg))
            threshold_msg = ''  # Use as flag so threshold only occurs once
        print(line)
    return threshold


def cluster_clones(comparisons, threshold):
    """ Cluster groups of clones together """
    clone_groups = []       # list of CloneGroup instances
    clone_indexes = {}      # clone_indexes[indiv] = (index for clone_groups)
    clone_index_count = -1  # to keep track of index for clone_groups
    for row in np.nditer(comparisons):
        # Stop iterating when reaching simularity threshold
        if row[C_MATCH_PERC] < threshold:
            break
        # Determine clone-index if one of indivs is already in list
        ind1 = str(row[C_IND1])
        ind2 = str(row[C_IND2])
        if ind1 in clone_indexes.keys() and ind2 in clone_indexes.keys():
            gi1 = clone_indexes[ind1]
            gi2 = clone_indexes[ind2]
            if gi1 != gi2:
                # Non-transitive link: merge the two groups (single-linkage)
                clone_groups[gi1].merge_group(clone_groups[gi2], row)
                for ind in clone_groups[gi2].indivs:
                    clone_indexes[ind] = gi1
                clone_groups[gi2] = None
        elif ind1 in clone_indexes.keys():
            clone_indexes[ind2] = clone_indexes[ind1]
            clone_groups[clone_indexes[ind1]].add_clone_from_row(row)
        elif ind2 in clone_indexes.keys():
            clone_indexes[ind1] = clone_indexes[ind2]
            clone_groups[clone_indexes[ind2]].add_clone_from_row(row)
        else:
            clone_index_count += 1
            clone_indexes[ind1] = clone_index_count
            clone_indexes[ind2] = clone_index_count
            clone_groups.append(CloneGroup(row))

    # Drop groups emptied by merging (None); preserves order otherwise
    return [group for group in clone_groups if group is not None]


def get_perc_genotyped(comparisons, n_loci):
    """ Per-sample % genotyped from ind*_snps columns and total loci count """
    counts = {}
    for row in np.nditer(comparisons):
        counts[str(row[C_IND1])] = int(row[C_IND1_SNPS])
        counts[str(row[C_IND2])] = int(row[C_IND2_SNPS])
    if n_loci:
        return {ind: 100.0 * c / n_loci for ind, c in counts.items()}, True
    # Fall back to absolute genotyped count if total loci unknown
    return counts, False


def write_pdf_report(comparisons, clone_groups, threshold, perc_genotyped,
                     is_percentage, method, has_popfile, pdf_filename):
    """ Build the multi-panel PDF report (portrait):
        left  = NJ tree (clone groups highlighted) + per-sample bars
        right = within/between similarity histogram (top) + zoomed (bottom)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import Bio.Phylo
    import Bio.Phylo.TreeConstruction

    # --- Build distance matrix (1 - similarity) and NJ tree ---
    names = sorted(perc_genotyped.keys())
    name_idx = {nm: i for i, nm in enumerate(names)}
    n = len(names)
    dist = np.zeros((n, n), dtype=np.float64)
    for row in np.nditer(comparisons):
        i = name_idx.get(str(row[C_IND1]))
        j = name_idx.get(str(row[C_IND2]))
        if i is None or j is None:
            continue
        d = max(0.0, 1.0 - float(row[C_MATCH_PERC]) / 100.0)
        dist[i, j] = dist[j, i] = d

    lower = [[dist[i][j] for j in range(i + 1)] for i in range(n)]
    dm = Bio.Phylo.TreeConstruction._DistanceMatrix(names=list(names),
                                                    matrix=lower)
    sys.stderr.write('Building NJ tree ({0} tips)...\n'.format(n))
    tree = Bio.Phylo.TreeConstruction.DistanceTreeConstructor().nj(dm)
    tree.ladderize()
    for clade in tree.get_nonterminals():
        clade.name = None

    # --- Clone-group colours (one per group; non-clones grey) ---
    n_groups = len(clone_groups)
    if n_groups <= 10:
        cmap = plt.cm.tab10
    elif n_groups <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.gist_rainbow
    label_colors = {}
    for gi, cg in enumerate(clone_groups):
        col = cmap(gi / max(n_groups - 1, 1))
        for ind in cg.indivs:
            label_colors[ind] = col

    # --- Figure layout: fixed A4 portrait page ---
    # Left half  = NJ tree (+ per-sample bars).
    # Right half = two stacked histograms (top 50% = all comparisons,
    #              bottom 50% = zoomed to the clonal threshold).
    n_tips = n
    fig_w, fig_h = 8.27, 11.69            # A4 portrait (inches)
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Left half: tree (with a thin bar panel added later, both within x < 0.5)
    ax_tree = fig.add_axes([0.07, 0.05, 0.27, 0.90])
    # Right half: two histograms, top half and bottom half
    ax_hist = fig.add_axes([0.575, 0.565, 0.385, 0.385])   # top-right
    ax_zoom = fig.add_axes([0.575, 0.075, 0.385, 0.385])   # bottom-right

    font_size = max(1.0, min(9, 600.0 / n_tips))
    Bio.Phylo.draw(tree, axes=ax_tree,
                   label_colors=label_colors,
                   label_func=lambda c: '  ' + c.name if c.name else '',
                   do_show=False)
    ax_tree.set_ylabel('')
    ax_tree.set_xlabel('Genetic distance (1 - {0} similarity)'.format(method))

    # Recover tip y-positions from the drawn labels
    tip_y = {}
    for text in ax_tree.texts:
        label = text.get_text().strip()
        if label in name_idx:
            tip_y[label] = text.get_position()[1]
            text.set_fontsize(font_size)
            if label in label_colors:
                x = text.get_position()[0]
                ax_tree.plot(x, tip_y[label], 'o', color=label_colors[label],
                             markersize=max(1.5, font_size),
                             markeredgewidth=0.3, markeredgecolor='black',
                             zorder=5, clip_on=False)

    # --- Per-sample bar panel, aligned to tree tips (overlaid axis) ---
    tree_pos = ax_tree.get_position()
    bar_w = 0.12
    ax_bar = fig.add_axes([tree_pos.x1 + 0.005, tree_pos.y0,
                           bar_w, tree_pos.height])
    if tip_y:
        ys = np.array([tip_y[nm] for nm in names if nm in tip_y])
        vals = np.array([perc_genotyped[nm] for nm in names if nm in tip_y])
        bar_colors = [label_colors.get(nm, '0.5')
                      for nm in names if nm in tip_y]
        ax_bar.barh(ys, vals, height=0.8, color=bar_colors,
                    edgecolor='none')
        ax_bar.set_ylim(ax_tree.get_ylim())
        ax_bar.set_yticks([])
        if is_percentage:
            ax_bar.set_xlim(max(0, vals.min() - 5), 100)
            ax_bar.set_title('% Genotyped', fontsize=8)
        else:
            ax_bar.set_xlim(0, vals.max() * 1.05)
            ax_bar.set_title('# SNPs genotyped', fontsize=8)
        ax_bar.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax_bar.tick_params(labelsize=6)
        for spine in ('top', 'right'):
            ax_bar.spines[spine].set_visible(False)

    # --- Histograms: within vs between population ---
    perc = comparisons[C_MATCH_PERC]
    pops = comparisons[C_POP].astype(str)
    is_between = np.array(['-' in p for p in pops])
    is_within = np.array([(p != 'NA' and '-' not in p) for p in pops])

    def draw_hist(ax, bins):
        if has_popfile and is_within.any():
            ax.hist(perc[is_between], bins=bins, color='0.6', alpha=0.8,
                    label='Between populations')
            ax.hist(perc[is_within], bins=bins, color='#1f6f6f', alpha=0.7,
                    label='Within populations')
        else:
            ax.hist(perc, bins=bins, color='0.5', alpha=0.85,
                    label='All pairs')
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        ax.set_xlabel('Genetic similarity (%)', fontsize=8)
        ax.set_ylabel('Pairwise comparisons', fontsize=8)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.tick_params(labelsize=7)

    def add_threshold_line(ax):
        ax.axvline(threshold, color='red', linestyle='--', linewidth=1)
        ax.text(threshold, ax.get_ylim()[1] * 0.98,
                ' threshold = {0}%'.format(round(threshold, 2)),
                color='red', fontsize=7, ha='left', va='top', rotation=90)

    def series_peak(lo, hi, width=0.5):
        """ Tallest single-series bar count within [lo, hi] """
        edges = np.arange(lo, hi + width, width)
        peak = 0
        if has_popfile and is_within.any():
            groups = (perc[is_between], perc[is_within])
        else:
            groups = (perc,)
        for g in groups:
            c, _ = np.histogram(g[(g >= lo) & (g <= hi)], bins=edges)
            if c.size and c.max() > peak:
                peak = int(c.max())
        return max(peak, 1)

    valid = perc[np.isfinite(perc)]
    lo_all = math.floor(valid.min()) if valid.size else 0
    draw_hist(ax_hist, np.linspace(lo_all, 100, max(20, int(100 - lo_all))))
    ax_hist.legend(fontsize=7, loc='upper left', framealpha=0.9)
    ax_hist.set_title('All pairwise comparisons', fontsize=9)
    add_threshold_line(ax_hist)

    # Zoomed: a small margin below the threshold up to 100%, with the y-axis
    # capped to the clone peak (counts at/above threshold) so it is legible
    # even though the bulk just below threshold is far taller (it is clipped).
    zoom_lo = max(lo_all, threshold - 1)
    draw_hist(ax_zoom, np.arange(zoom_lo, 100.0 + 0.5, 0.5))
    ax_zoom.set_xlim(zoom_lo, 100)
    clone_peak = series_peak(threshold, 100.0)
    ax_zoom.set_ylim(0, clone_peak * 1.4)
    ax_zoom.set_title('Zoomed to clonal threshold', fontsize=9)
    add_threshold_line(ax_zoom)

    fig.savefig(pdf_filename, format='pdf')
    plt.close(fig)
    sys.stderr.write('PDF report written to {0}\n'.format(pdf_filename))


def derive_pdf_filename(vcf_filename, input_filename, output_filename,
                        pdf_output):
    """ Decide the PDF output filename """
    if pdf_output:
        return pdf_output
    base = None
    for candidate in (output_filename, vcf_filename, input_filename):
        if candidate:
            base = os.path.splitext(candidate)[0]
            break
    if base is None:
        base = 'clone_detect'
    return base + '_clones.pdf'


def main(vcf_filename, input_filename, output_filename, pop_filename,
         threshold, method=DEFAULT_METHOD, make_pdf=True, pdf_output=None):

    print('###1 - Pairwise comparisons of all individuals')

    n_loci = None
    perc_genotyped = None
    is_percentage = False

    # Input data (vcf_file or input_file)
    if vcf_filename:
        if pop_filename:
            indivs_pops = get_pop_assignments_from_popfile(pop_filename)
        else:
            indivs_pops = {}
        need_ad = (method == 'single-read')
        data = load_vcf_arrays(vcf_filename, need_ad=need_ad)
        if data['multiallelic']:
            sys.stderr.write('Warning: multiallelic sites present; `dosage` '
                             'and `single-read` treat alt-allele counts only.\n')
        comparisons = build_comparisons(data, indivs_pops, method)
        n_loci = data['n_loci']
    elif input_filename:
        comparisons = get_pairwise_comparisons_from_input_file(input_filename)
        n_loci = read_nloci_from_input_file(input_filename)
    else:
        sys.exit('Error: Please provide either a vcf_file or input_file.')

    # Output data (as input_file)
    if output_filename:
        save_pairwise_comparisons_to_input_file(output_filename, comparisons,
                                                n_loci=n_loci)

    print('\n###2 - Histogram (of pairwise genetic similarities)')
    output_ascii_hist(comparisons[C_MATCH_PERC], HIST_RANGE)

    print('\n###3 - List of highest matches')
    threshold = output_highest_matches(comparisons, float(threshold))

    print('\n###4 - Clonal groups (threshold: {0})'.format(threshold))
    clone_groups = cluster_clones(comparisons, float(threshold))
    for clone_group in clone_groups:
        print(clone_group.get_formatted_clone_info())

    print(('\n###5 - Individuals to remove from dataset (retaining indiv'
           ' with least amount of missing data)'))
    for clone_group in clone_groups:
        print('\n'.join(clone_group.get_samples_to_remove()))

    # PDF report (section 6)
    if make_pdf:
        perc_genotyped, is_percentage = get_perc_genotyped(comparisons, n_loci)
        pdf_filename = derive_pdf_filename(vcf_filename, input_filename,
                                           output_filename, pdf_output)
        print('\n###6 - PDF report: {0}'.format(pdf_filename))
        try:
            write_pdf_report(comparisons, clone_groups, float(threshold),
                             perc_genotyped, is_percentage, method,
                             bool(pop_filename) or _has_pop_info(comparisons),
                             pdf_filename)
        except Exception as e:
            sys.stderr.write('Warning: PDF report failed ({0})\n'.format(e))


def _has_pop_info(comparisons):
    """ Whether the comparisons contain population assignments (not all NA) """
    pops = comparisons[C_POP].astype(str)
    return bool(np.any(pops != 'NA'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-p', '--pop', dest='pop_filename',
                        metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations (to accompany `.vcf` file)')
    parser.add_argument('-i', '--input', dest='input_filename',
                        metavar='compare_file',
                        help='input file (csv) with previously \
                        calculated pairwise comparisons (using the \
                        `--outputfile` option)')
    parser.add_argument('-o', '--output', dest='output_filename',
                        metavar='compare_file',
                        help='output file (csv) for all pairwise comparisons \
                        (can later be used as input with `--inputfile`)')
    parser.add_argument('-t', '--threshold', dest='threshold',
                        metavar='threshold', default=0.0,
                        help='manual similarity threshold (e.g. `94.5` means \
                        at least 94.5 percent allelic similarity for \
                        individuals to be considered clones)')
    parser.add_argument('-m', '--method', dest='method', default=DEFAULT_METHOD,
                        choices=METHODS,
                        help='similarity measure (default: ibs, identical to \
                        the original allelic similarity). PLINK/ANGSD \
                        correspondence: ibs = PLINK `--distance 1-ibs \
                        flat-missing` (and default `--distance`/allele-ct); \
                        dosage = PLINK `--make-rel cov` (squared-Euclidean, \
                        NOT the standardised GRM); single-read = ANGSD \
                        `-doIBS 1 -makeMatrix 1` (AD-weighted if available); \
                        het-masked has no PLINK/ANGSD equivalent. `het-masked` \
                        and `dosage` give a sharper clone boundary; \
                        `single-read` is structure/topology-oriented and NOT \
                        recommended for clone thresholds (see module docstring)')
    parser.add_argument('--pdf-output', dest='pdf_output', default=None,
                        metavar='pdf_file',
                        help='filename for the PDF report (default: derived \
                        from the vcf/output/input basename)')
    parser.add_argument('--no-pdf', dest='no_pdf', action='store_true',
                        help='do not generate the PDF report (text output only)')
    args = parser.parse_args()
    main(args.vcf_filename, args.input_filename, args.output_filename,
         args.pop_filename, args.threshold, method=args.method,
         make_pdf=not args.no_pdf, pdf_output=args.pdf_output)
