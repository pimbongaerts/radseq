#!/usr/bin/env python
"""
Attempts to identify groups of clones in a dataset. The script
(1) conducts pairwise comparisons (allelic similarity) for all individuals in
a `.vcf` file,
(2) produces a histogram of genetic similarities,
(3) lists the highest matches to assess for a potential clonal threshold,
(4) clusters the groups of clones based on a particular threshold (supplied or
roughly inferred), and
(5) lists the clonal individuals that can be removed from the dataset
(so that one individual with the least amount of missing data remains).
If optional popfile is given, then clonal groups are sorted by population.
Note: Firstly, the script is run with a `.vcf` file and an optional popfile
to produce an output file (e.g. `python3 vcf_clone_detect.py.py --vcf
vcf_file.vcf --pop pop_file.txt --output compare_file.csv`). Secondly, it can
be rerun using the precalculated similarities under different thresholds
(e.g. `python3 vcf_clone_detect.py.py --input compare_file.csv
--threshold 94.5`)
"""
import sys
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
        info = '{0}: {1} ({2})'.format('-'.join(self.pops),
                                       ', '.join(self.indivs),
                                       score_range)
        return info

    def get_samples_to_remove(self):
        return self.indivs - set([self.best_indiv])


def get_snp_match(genotype1, genotype2):
    """ Get match value for two genotypes (one SNP) """
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
        indiv = cols[0]
        pop = cols[1]
        indivs_pops[indiv] = pop
    return indivs_pops


def get_genotypes_from_vcf(vcf_filename):
    """ Read genotypes from vcf and store in dictionary """
    vcf_file = open(vcf_filename, 'r')
    individuals = {}
    genotypes = {}

    # Iterate through vcf file and store individual names and genotypes
    for line in vcf_file:
        cols = line.split()
        # Store name of individuals in temporary dictionary
        if line[:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
            for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
                individuals[x] = cols[x]
        # Store genotype for each individual in dictionary
        elif line[0] != HEADER_CHAR:
            for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
                genotypes.setdefault(individuals[x], []).append(cols[x])
    vcf_file.close()
    return genotypes


def get_match_info_as_row(individual1, genotypes1,
                          individual2, genotypes2, indivs_pops):
    """ Get match value for two genotypes (one SNP) """
    ind1_snps = ind2_snps = both_snps = match = 0
    for x in range(0, len(genotypes1)):
        genotype1 = genotypes1[x][:3]
        genotype2 = genotypes2[x][:3]
        if genotype1[0] != '.' and genotype2[0] != '.':  # both genotyped
            ind1_snps += 1
            ind2_snps += 1
            both_snps += 1
            match += get_snp_match(genotype1, genotype2)
        elif genotype1[0] != '.' and genotype2[0] == '.':  # genotype of ind2 missing
            ind1_snps += 1
        elif genotype1[0] == '.' and genotype2[0] != '.':  # genotype of ind1 missing
            ind2_snps += 1
        elif genotype1[0] == '.' and genotype2[0] == '.':   # both missing
            pass # do nothing
        else:
            sys.exit('Unexpected genotype: {0} {1}'.format(genotype1, genotype2))

    match_perc = round((match / both_snps) * 100, 2)

    if individual1 in indivs_pops and individual2 in indivs_pops:
        # Define whether comparison is within or between pops
        if indivs_pops[individual1] == indivs_pops[individual2]:
            pop_group = indivs_pops[individual1]
        else:
            pops = sorted([indivs_pops[individual1],
                           indivs_pops[individual2]])
            pop_group = '-'.join(pops)
    else:
        # Not applicable as indiv(s) not in popfile or popfile not provided
        pop_group = 'NA'

    return (str(individual1), str(individual2), int(ind1_snps),
            int(ind2_snps), int(both_snps), float(match),
            float(match_perc), str(pop_group))


def get_all_pairwise_comparisons(genotypes, indivs_pops):
    """ Conduct pairwise comparisons between all individuals """
    unique_pairs = int((math.pow(len(genotypes), 2) - len(genotypes)) / 2)
    comparisons = np.zeros(unique_pairs, dtype=COMPARISONS_DTYPES)

    index = 0
    for individual1, individual2 in itertools.combinations(genotypes, 2):

        comparisons[index] = get_match_info_as_row(individual1,
                                                   genotypes[individual1],
                                                   individual2,
                                                   genotypes[individual2],
                                                   indivs_pops)
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


def save_pairwise_comparisons_to_input_file(output_filename, comparisons):
    """ Load pairwise comparisons between all individuals to output_file """
    np.savetxt(output_filename, comparisons, fmt=OUTPUT_FILE_FORMAT,
               header=OUTPUT_FILE_HEADER)
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
            if clone_indexes[ind1] != clone_indexes[ind2]:
                sys.exit('Error: clone matrix not properly sorted')
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

    return clone_groups


def main(vcf_filename, input_filename, output_filename, pop_filename,
         threshold):

    print('###1 - Pairwise comparisons of all individuals')

    # Input data (vcf_file or input_file)
    if vcf_filename:
        # Population assignments
        if pop_filename:
            indivs_pops = get_pop_assignments_from_popfile(pop_filename)
        else:
            indivs_pops = []
        # Genotypes
        genotypes = get_genotypes_from_vcf(vcf_filename)
        # Comparisons
        comparisons = get_all_pairwise_comparisons(genotypes, indivs_pops)
    elif input_filename:
        # Comparisons from file
        comparisons = get_pairwise_comparisons_from_input_file(input_filename)
    else:
        sys.exit('Error: Please provide either a vcf_file or input_file.')

    # Output data (as input_file)
    if output_filename:
        save_pairwise_comparisons_to_input_file(output_filename, comparisons)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
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
    args = parser.parse_args()
    main(args.vcf_filename, args.input_filename, args.output_filename,
         args.pop_filename, args.threshold)
