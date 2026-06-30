#!/usr/bin/env python
"""
Identify groups of clones (near-identical samples) in a `.vcf` dataset. The
script:
 (1) computes a pairwise genetic similarity between every pair of individuals,
 (2) prints a text histogram of those similarities,
 (3) lists the highest-similarity pairs and, unless --threshold is given,
     infers a candidate clonal threshold (see "Threshold inference" below),
 (4) clusters individuals into clonal groups at the threshold,
 (5) lists, per group, the members to remove so that the member with the most
     genotyped loci is retained, and
 (6) writes an A4 PDF report: a neighbour-joining tree with clonal groups
     coloured and a per-sample % genotyped bar panel, plus similarity
     histograms (full, and zoomed to the threshold).

Usage is two-step. First compute and save the pairwise comparisons from a VCF:
    vcf_clone_detect.py --vcf data.vcf --pop pops.txt --output compare.csv
Then re-threshold without recomputing by reading that file back:
    vcf_clone_detect.py --input compare.csv --threshold 94.5
A popfile (sample, population per line) sets the matrix/tree order and splits
the histograms into within- vs between-population comparisons.

Similarity measures (--method)
------------------------------
Each measure is computed only over sites genotyped in both individuals of a
pair, and reported as a 0-100 % similarity (100 = identical). The tree distance
is 1 - similarity/100. da, db are alt-allele dosages (0, 1, 2).

  ibs          (default) Mean allele-sharing per site: score 1 if the two
               genotypes are identical, 0.5 if they share one allele, 0 if they
               share none (a heterozygote pair scores 1). Equals PLINK
               `--distance 1-ibs flat-missing` (and the default
               `--distance`/`allele-ct`, rescaled).

  het-masked   As ibs but using only sites where BOTH individuals are
               homozygous; the per-site score is 1 if the genotypes are
               identical, else 0. No PLINK/ANGSD equivalent.

  dosage       Mean squared dosage difference per site, 1 - (da-db)^2 / dmax^2
               (dmax = maximum dosage, 2 for biallelic data). Squared-Euclidean
               distance on dosage; equals PLINK `--make-rel cov` via
               sum((da-db)^2) = M*(Cjj+Ckk-2*Cjk). Not the standardised GRM.

  single-read  Mean of p1 + p2 - 2*p1*p2 per site, where p is the alt-allele
               fraction: alt/(ref+alt) from the AD field when present, otherwise
               the genotype dosage/2. This is the expected mismatch when one
               allele is drawn at random from each individual; it emulates ANGSD
               `-doIBS 1 -makeMatrix 1`. With called genotypes a heterozygote
               pair has an expected per-site mismatch of 0.5. (ANGSD itself
               needs read-level data; this approximates it from a VCF.)

With --method ibs the result matches the previous version of this script
exactly. PLINK/ANGSD equivalences are exact only on identical SNP sets: this
script uses pairwise-complete sites while PLINK/ANGSD apply their own
missing-data scaling, so values can differ slightly when data are missing.

Threshold inference (when --threshold is not given)
---------------------------------------------------
Pairs are sorted by descending similarity. Considering only pairs at or above
85 %, the largest drop in similarity between two consecutive pairs is taken as
the clone/non-clone break. The threshold is the integer part of the similarity
just above that drop if that integer lies within the drop, otherwise the
midpoint of the drop. This heuristic assumes clones form a tight cluster at the
top of the distribution and is calibrated to the ibs %-scale; for the other
measures, or whenever the inferred value looks wrong, set --threshold.

Multiple cryptic lineages
-------------------------
When a dataset contains several genetic lineages, a single clonal threshold is
unreliable: each lineage has a different similarity distribution (different
informative-SNP counts and per-pair denominators), so clones are best called
within each lineage on its own scale. Two ways to do this:
  - Supply a 3-column popfile (sample, population, lineage): clone detection
    runs independently within each lineage, each with its own histogram and
    inferred/supplied threshold.
  - Or give `-k/--n-lineages K`: the dataset is split into K lineages by
    average-linkage clustering of this script's own distance matrix (which is
    clone-robust, unlike STRUCTURE/SNAPCLUST), the inferred 3-column popfile is
    written to `<base>_lineages.txt` for review, and detection then runs per
    lineage.
Each lineage gets its own comparison CSV and PDF; a combined list of all
samples to remove across the dataset is written to `<base>_clones_remove.txt`.
`--no-lineage` forces a single global run. Only within-lineage pairs are ever
eligible to be clones. The comparison CSV is always written (default name from
the input basename if `--output` is omitted).

Output columns (--output / --input CSV)
---------------------------------------
ind1, ind2, ind1_snps, ind2_snps, both_snps, match, match_perc, pop.
ind1_snps/ind2_snps: each sample's genotyped-loci count. both_snps: number of
sites the chosen measure used for the pair (genotyped in both; for het-masked,
homozygous in both; for AD-based single-read, with reads in both). match: summed
per-site similarity score over those sites; match_perc = round(100 * match /
both_snps, 2). pop: the shared population, 'popA-popB' for a between-population
pair, or 'NA'.
"""
import sys
import os
import argparse
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

# Popfile columns (3rd column = genetic lineage is optional)
COL_INDIV = 0
COL_POP = 1
COL_LINEAGE = 2

# Default output suffixes (used when --output is not given)
COMPARE_SUFFIX = '_compare.csv'
PDF_SUFFIX = '_clones.pdf'
REMOVE_SUFFIX = '_clones_remove.txt'
LINEAGES_SUFFIX = '_lineages.txt'

# PDF report font sizes (standardised across the three panels)
FS_TITLE = 9
FS_LABEL = 8
FS_TICK = 7

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
    """ Allele-sharing score for two genotypes at one SNP: 1 if identical,
    0.5 if they share an allele, 0 otherwise. Used to build the per-state score
    matrix for the vectorised `ibs` computation. """
    if genotype1 == genotype2:
        match_score = 1
    elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2] or \
            genotype1[0] == genotype2[2] or genotype1[2] == genotype2[0]:
        match_score = 0.5
    else:
        match_score = 0
    return match_score


def get_assignments_from_popfile(pop_filename):
    """ Read a (tsv/csv) popfile and return (indivs_pops, indivs_lineages).
    Column 1 = sample, column 2 = population, optional column 3 = genetic
    lineage. indivs_lineages is empty when no 3rd column is present. """
    indivs_pops = {}
    indivs_lineages = {}
    with open(pop_filename, 'r') as pop_file:
        for line in pop_file:
            cols = line.rstrip().replace(',', ' ').split()
            if not cols:
                continue
            indiv = cols[COL_INDIV]
            indivs_pops[indiv] = cols[COL_POP] if len(cols) > COL_POP else 'NA'
            if len(cols) > COL_LINEAGE:
                indivs_lineages[indiv] = cols[COL_LINEAGE]
    return indivs_pops, indivs_lineages


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


def compute_pair_matrices(data, method):
    """ Compute the global NxN match/both matrices once, plus per-sample
    genotyped-loci counts (== ind*_snps). """
    match_mat, both_mat = _match_denom_matrices(data, method)
    geno_count = (data['state_code'] >= 0).sum(axis=0).astype(np.int64)
    return match_mat, both_mat, geno_count


def comparisons_for_indices(names, match_mat, both_mat, geno_count,
                            indivs_pops, indices):
    """ Build the structured comparisons array for all pairs within `indices`
    (column positions into the global matrices). With indices = range(n) this
    reproduces the original all-pairs output exactly. """
    indices = list(indices)
    npairs = len(indices) * (len(indices) - 1) // 2
    comparisons = np.zeros(npairs, dtype=COMPARISONS_DTYPES)
    index = 0
    for i, j in itertools.combinations(indices, 2):
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
    return comparisons


def delimit_lineages(match_mat, both_mat, names, k):
    """ Partition samples into k lineages by average-linkage (UPGMA) clustering
    of the pairwise distance (1 - similarity). Returns {sample: 'L<n>'}.
    Distance is clone-robust, so clones cluster within their lineage rather than
    distorting the deep splits. """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    n = len(names)
    with np.errstate(divide='ignore', invalid='ignore'):
        dist = np.where(both_mat > 0, 1.0 - match_mat / both_mat, 1.0)
    dist = np.clip(dist, 0.0, None)
    dist = (dist + dist.T) / 2.0            # enforce exact symmetry
    np.fill_diagonal(dist, 0.0)
    linkage_matrix = linkage(squareform(dist, checks=False), method='average')
    labels = fcluster(linkage_matrix, t=k, criterion='maxclust')
    return {names[i]: 'L{0}'.format(labels[i]) for i in range(n)}


def write_lineage_popfile(filename, names, indivs_pops, lineages):
    """ Write an inferred 3-column popfile (sample, population, lineage) for
    review / re-use with --pop. """
    with open(filename, 'w') as out:
        for nm in names:
            out.write('{0}\t{1}\t{2}\n'.format(nm, indivs_pops.get(nm, 'NA'),
                                               lineages[nm]))
    sys.stderr.write('Inferred lineage assignments written to {0}\n'.format(
        filename))


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


def output_ascii_hist(raw_values, bin_values):
    """ Plot a text-based histogram """
    values, bins = np.histogram(raw_values, bins=bin_values)
    output_lines = []
    graph_multiplier = 1
    lower_bound_flag = False
    lower_bound = upper_bound = 0
    for index, value in enumerate(values):
        if value > 0:
            if not lower_bound_flag:
                lower_bound_flag = True
                lower_bound = bins[index]
            upper_bound = bins[index]      # last populated bin (>=1 handled)
        if lower_bound_flag:
            graph_bar = '*' * int(value * graph_multiplier)
            if len(graph_bar) > 70:
                graph_bar = graph_bar[:69] + '#'

            output_lines.append(
                '{:3d} {:7d} {:70s}'.format(bins[index], value,
                                            graph_bar[:70]))
    if output_lines:
        print('\n'.join(output_lines[:(upper_bound - lower_bound + 2)]))


def output_highest_matches(comparisons, threshold):
    """ Print the highest-similarity pairs and return the threshold to use.

    If `threshold` is 0, it is inferred: scanning pairs from highest similarity
    down to DEF_THRESHOLD, the largest drop between two consecutive pairs is
    taken as the clone/non-clone break. The returned threshold is the integer
    part of the similarity just above that drop when that integer falls within
    the drop, otherwise the midpoint of the drop. A non-zero `threshold` is
    used as given (and reported as a manual threshold). """
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
            # Drop spans no integer: use its midpoint (a value inside the gap)
            threshold = (highest_diff_max_perc + highest_diff_min_perc) / 2

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
    tree_bottom, tree_top = 0.05, 0.95
    ax_tree = fig.add_axes([0.07, tree_bottom, 0.27, tree_top - tree_bottom])
    # Right half: two histograms. Top panel's top aligns with the tree top;
    # bottom panel's bottom (its x-axis) aligns with the tree bottom x-axis.
    hist_h = 0.385
    ax_hist = fig.add_axes([0.575, tree_top - hist_h, 0.385, hist_h])
    ax_zoom = fig.add_axes([0.575, tree_bottom, 0.385, hist_h])

    font_size = max(1.0, min(9, 600.0 / n_tips))
    Bio.Phylo.draw(tree, axes=ax_tree,
                   label_colors=label_colors,
                   label_func=lambda c: '  ' + c.name if c.name else '',
                   do_show=False)
    ax_tree.set_ylabel('')
    ax_tree.set_xlabel('Genetic distance (1 - {0} similarity)'.format(method),
                       fontsize=FS_LABEL)
    ax_tree.tick_params(labelsize=FS_TICK)

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
            ax_bar.set_title('% Genotyped', fontsize=FS_LABEL)
        else:
            ax_bar.set_xlim(0, vals.max() * 1.05)
            ax_bar.set_title('# SNPs genotyped', fontsize=FS_LABEL)
        ax_bar.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax_bar.tick_params(labelsize=FS_TICK)
        for spine in ('top', 'right'):
            ax_bar.spines[spine].set_visible(False)

    # --- Histograms: within vs between population ---
    perc = comparisons[C_MATCH_PERC]
    pops = comparisons[C_POP].astype(str)
    is_between = np.array(['-' in p for p in pops])
    is_within = np.array([(p != 'NA' and '-' not in p) for p in pops])

    def draw_hist(ax, bins, log=False):
        if has_popfile and is_within.any():
            ax.hist(perc[is_between], bins=bins, color='0.6', alpha=0.8,
                    log=log, label='Between populations')
            ax.hist(perc[is_within], bins=bins, color='#1f6f6f', alpha=0.7,
                    log=log, label='Within populations')
        else:
            ax.hist(perc, bins=bins, color='0.5', alpha=0.85, log=log,
                    label='All pairs')
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        ax.set_xlabel('Genetic similarity (%)', fontsize=FS_LABEL)
        ax.set_ylabel('Pairwise comparisons', fontsize=FS_LABEL)
        if not log:
            ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.tick_params(labelsize=FS_TICK)

    def add_threshold_line(ax):
        ax.axvline(threshold, color='red', linestyle='--', linewidth=1)
        ax.text(threshold, ax.get_ylim()[1] * 0.98,
                ' threshold = {0}%'.format(round(threshold, 2)),
                color='red', fontsize=FS_TICK, ha='left', va='top',
                rotation=90)

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
    # Top panel: log y-scale so the few clonal pairs on the right are visible
    # against the much larger non-clone bulk.
    draw_hist(ax_hist, np.linspace(lo_all, 100, max(20, int(100 - lo_all))),
              log=True)
    ax_hist.legend(fontsize=FS_TICK, loc='upper left', framealpha=0.9)
    ax_hist.set_title('All pairwise comparisons', fontsize=FS_TITLE)
    add_threshold_line(ax_hist)

    # Bottom panel: zoom from the upper tail of the sub-threshold peak up to
    # 100%, with the y-axis capped to the clone-region peak so clones are
    # legible (the much taller bulk tail clips at the top).
    nonclone = perc[(perc < threshold) & np.isfinite(perc)]
    tail_edge = np.percentile(nonclone, 95) if nonclone.size else threshold - 1
    zoom_lo = min(tail_edge, threshold - 0.5)        # show some sub-threshold tail
    zoom_lo = max(zoom_lo, threshold - 8, lo_all)    # but not the whole distribution
    draw_hist(ax_zoom, np.arange(zoom_lo, 100.0 + 0.5, 0.5))
    ax_zoom.set_xlim(zoom_lo, 100)
    clone_peak = series_peak(threshold, 100.0)
    ax_zoom.set_ylim(0, clone_peak * 1.4)
    ax_zoom.set_title('Zoomed to clonal threshold', fontsize=FS_TITLE)
    add_threshold_line(ax_zoom)

    fig.savefig(pdf_filename, format='pdf')
    plt.close(fig)
    sys.stderr.write('PDF report written to {0}\n'.format(pdf_filename))


def derive_base(vcf_filename, input_filename, output_filename):
    """ Basename (no extension) for deriving default output filenames """
    for candidate in (output_filename, vcf_filename, input_filename):
        if candidate:
            return os.path.splitext(candidate)[0]
    return 'clone_detect'


def analyze_set(comparisons, threshold, method, n_loci, make_pdf,
                pdf_filename, output_filename, has_popfile):
    """ Run sections 1-6 on one set of comparisons (whole dataset or one
    lineage). Returns (threshold_used, clone_groups, removed_samples). """
    # 1 - always save the comparisons (skip only when no path is available)
    if output_filename:
        save_pairwise_comparisons_to_input_file(output_filename, comparisons,
                                                n_loci=n_loci)
        print('###1 - Pairwise comparisons: written to {0}'.format(
            output_filename))
    else:
        print('###1 - Pairwise comparisons: {0} pairs'.format(comparisons.size))

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
    removed = []
    for clone_group in clone_groups:
        samples = clone_group.get_samples_to_remove()
        removed.extend(samples)
        print('\n'.join(samples))

    if make_pdf:
        perc_genotyped, is_percentage = get_perc_genotyped(comparisons, n_loci)
        print('\n###6 - PDF report: {0}'.format(pdf_filename))
        try:
            write_pdf_report(comparisons, clone_groups, float(threshold),
                             perc_genotyped, is_percentage, method,
                             has_popfile, pdf_filename)
        except Exception as e:
            sys.stderr.write('Warning: PDF report failed ({0})\n'.format(e))
    return threshold, clone_groups, removed


def main(vcf_filename, input_filename, output_filename, pop_filename,
         threshold, method=DEFAULT_METHOD, make_pdf=True, pdf_output=None,
         n_lineages=None, no_lineage=False):

    base = derive_base(vcf_filename, input_filename, output_filename)

    # --- Re-run from a precomputed comparison file (single global set) ---
    if input_filename and not vcf_filename:
        comparisons = get_pairwise_comparisons_from_input_file(input_filename)
        n_loci = read_nloci_from_input_file(input_filename)
        pdf = pdf_output or (base + PDF_SUFFIX)
        analyze_set(comparisons, threshold, method, n_loci, make_pdf, pdf,
                    output_filename, _has_pop_info(comparisons))
        return

    if not vcf_filename:
        sys.exit('Error: Please provide either a vcf_file or input_file.')

    # --- Compute pairwise matrices from the VCF (once) ---
    indivs_pops, indivs_lineages = ({}, {})
    if pop_filename:
        indivs_pops, indivs_lineages = get_assignments_from_popfile(pop_filename)
    need_ad = (method == 'single-read')
    data = load_vcf_arrays(vcf_filename, need_ad=need_ad)
    if data['multiallelic']:
        sys.stderr.write('Warning: multiallelic sites present; `dosage` '
                         'and `single-read` treat alt-allele counts only.\n')
    names = data['names']
    n_loci = data['n_loci']
    name_to_col = {nm: i for i, nm in enumerate(names)}
    match_mat, both_mat, geno_count = compute_pair_matrices(data, method)
    has_pop = bool(indivs_pops)

    # --- Decide lineage assignments (supplied, inferred, or none) ---
    lineages = None
    if not no_lineage:
        if indivs_lineages:
            lineages = {nm: indivs_lineages.get(nm, 'NA') for nm in names}
        elif n_lineages and n_lineages > 1:
            lineages = delimit_lineages(match_mat, both_mat, names, n_lineages)
            write_lineage_popfile(base + LINEAGES_SUFFIX, names, indivs_pops,
                                  lineages)

    # --- Single global run (no lineages) ---
    if not lineages:
        comparisons = comparisons_for_indices(names, match_mat, both_mat,
                                              geno_count, indivs_pops,
                                              range(len(names)))
        out_csv = output_filename or (base + COMPARE_SUFFIX)
        pdf = pdf_output or (base + PDF_SUFFIX)
        analyze_set(comparisons, threshold, method, n_loci, make_pdf, pdf,
                    out_csv, has_pop)
        return

    # --- Per-lineage runs ---
    groups = {}
    for nm in names:
        groups.setdefault(lineages[nm], []).append(nm)
    sys.stderr.write('Lineage mode: {0} lineages ({1})\n'.format(
        len(groups), ', '.join('{0}={1}'.format(k, len(v))
                               for k, v in sorted(groups.items()))))
    all_removed = []
    for lineage in sorted(groups):
        members = groups[lineage]
        print('\n========== Lineage {0} ({1} samples) =========='.format(
            lineage, len(members)))
        if len(members) < 2:
            print('(skipped: fewer than 2 samples)')
            continue
        idx = [name_to_col[nm] for nm in members]
        comparisons = comparisons_for_indices(names, match_mat, both_mat,
                                              geno_count, indivs_pops, idx)
        out_csv = '{0}_{1}{2}'.format(base, lineage, COMPARE_SUFFIX)
        pdf = '{0}_{1}{2}'.format(base, lineage, PDF_SUFFIX)
        _, _, removed = analyze_set(comparisons, threshold, method, n_loci,
                                    make_pdf, pdf, out_csv, has_pop)
        all_removed.extend(removed)

    # --- Combined removal list across lineages ---
    print('\n========== Combined samples to remove (all lineages) ==========')
    all_removed = sorted(set(all_removed))
    print('\n'.join(all_removed))
    remove_file = base + REMOVE_SUFFIX
    with open(remove_file, 'w') as out:
        out.write('\n'.join(all_removed) + ('\n' if all_removed else ''))
    sys.stderr.write('Combined removal list written to {0}\n'.format(
        remove_file))


def _has_pop_info(comparisons):
    """ Whether the comparisons contain population assignments (not all NA) """
    pops = comparisons[C_POP].astype(str)
    return bool(np.any(pops != 'NA'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='vcf_file',
                        help='input VCF file with SNP data')
    parser.add_argument('-p', '--pop', dest='pop_filename',
                        metavar='pop_file',
                        help='population file (sample, population[, lineage] per '
                             'line); a 3rd lineage column enables per-lineage '
                             'detection')
    parser.add_argument('-i', '--input', dest='input_filename',
                        metavar='compare_file',
                        help='read pairwise comparisons from a prior --output '
                             'CSV instead of a VCF')
    parser.add_argument('-o', '--output', dest='output_filename',
                        metavar='compare_file',
                        help='comparisons CSV (default: derived from input '
                             'basename; always written)')
    parser.add_argument('-t', '--threshold', dest='threshold',
                        metavar='threshold', default=0.0,
                        help='minimum %% similarity to call two samples clones '
                             '(default: inferred per set)')
    parser.add_argument('-m', '--method', dest='method', default=DEFAULT_METHOD,
                        choices=METHODS,
                        help='similarity measure (default: ibs); see above')
    parser.add_argument('-k', '--n-lineages', dest='n_lineages', type=int,
                        default=None, metavar='K',
                        help='delimit K lineages from the distance matrix and '
                             'run clone detection within each (writes an '
                             'inferred 3-col popfile); ignored if --pop already '
                             'has a lineage column')
    parser.add_argument('--no-lineage', dest='no_lineage', action='store_true',
                        help='force a single global run even with a 3-col '
                             'popfile or --n-lineages')
    parser.add_argument('--pdf-output', dest='pdf_output', default=None,
                        metavar='pdf_file',
                        help='PDF report filename (default: from input '
                             'basename; per-lineage names are auto-derived)')
    parser.add_argument('--no-pdf', dest='no_pdf', action='store_true',
                        help='skip the PDF report (text output only)')
    args = parser.parse_args()
    main(args.vcf_filename, args.input_filename, args.output_filename,
         args.pop_filename, args.threshold, method=args.method,
         make_pdf=not args.no_pdf, pdf_output=args.pdf_output,
         n_lineages=args.n_lineages, no_lineage=args.no_lineage)
