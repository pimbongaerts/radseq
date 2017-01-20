#!/usr/bin/env python
"""
Calculates Genetic Distance (Hamming / p-distance) for each pair of individuals
in a `.vcf` file and outputs as matrix. Popfile is supplied to indicate order
in matrix.
"""
import sys
import argparse
import operator
import itertools
from collections import Counter

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


VCF_HEADER_CHAR = '#'
VCF_INDIV_HEADER = '#CHROM'
ERROR_CODE = -1
SELF_MATCH = 0

SEPARATOR = '\t'
MIN_THRESHOLD = 100         # Minimum number of SNPs without missing data


def get_genotypes_from_vcf(vcf_filename):
    """ Read genotypes from VCF into dict of lists (indivs: genotypes) """
    indivs_gts = {}
    colnrs_indivs = {}
    # Iterate through VCF and store genotypes
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        cols = line.rstrip().split()
        if line[0:len(VCF_INDIV_HEADER)] == VCF_INDIV_HEADER:
            # Extract indiv names from header
            colnrs_indivs = {i: indiv for (i, indiv) in
                             enumerate(cols[9:], start=9)}
            # Create dict item with list for indiv
            indivs_gts = {indiv: [] for indiv in cols[9:]}
        elif line[0] != VCF_HEADER_CHAR:
            # Extract genotypes for all individuals (current SNP)
            for index, genotype in enumerate(cols[9:], start=9):
                indivs_gts[colnrs_indivs[index]].append(genotype)
    vcf_file.close()
    return indivs_gts


def compare_indivs(indiv1_gts, indiv2_gts):
    """ Calculate genetic distance between 2 individuals """
    matches_count = total_count = 0

    # Evaluate match across each individual SNP
    for x in range(0, len(indiv1_gts)):
        genotype1 = indiv1_gts[x][:3]
        genotype2 = indiv2_gts[x][:3]

        # Only consider SNPs that are genotyped for both individuals
        if genotype1[0] != '.' and genotype2[0] != '.':
            if genotype1 == genotype2:
                matches_count += 1
            elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2]:
                matches_count += 0.5
            total_count += 1

    # Evaluate overall match if above threshold
    if total_count > MIN_THRESHOLD:
        # percentage_match = round((matches_count/total_count) * 100, 2)
        return round(1 - (matches_count / total_count), 4)
    else:
        return ERROR_CODE


def get_comparison_key(pairs_gds, indiv1, indiv2):
    """ Obtain comparison key from dict for pair of individuals """
    comparison_key1 = '{0}:{1}'.format(indiv1, indiv2)
    comparison_key2 = '{0}:{1}'.format(indiv2, indiv1)
    if comparison_key1 in pairs_gds:
        return comparison_key1
    elif comparison_key2 in pairs_gds:
        return comparison_key2
    else:
        return sys.exit('Unexpected error')


def get_match_score(pairs_gds, indiv1, indiv2):
    """ Obtain comparison key from dict for pair of individuals """
    if indiv1 == indiv2:
        return str(SELF_MATCH)
    else:
        comparison_key = get_comparison_key(pairs_gds, indiv1, indiv2)
        return str(pairs_gds[comparison_key])


def main(vcf_filename, pop_filename):
    # Read order of individuals from popfile
    with open(pop_filename) as file:
        indiv_order = [line.replace(',', ' ').split()[0] for line in file]

    # Read genotypes for all individual and loci into a dict
    indivs_gts = get_genotypes_from_vcf(vcf_filename)

    # Conduct a pairwise comparison between all individuals
    pairs_gds = {}
    for indiv1, indiv2 in itertools.combinations(indivs_gts, 2):
        comparison_key = '{0}:{1}'.format(indiv1, indiv2)
        pairs_gds[comparison_key] = compare_indivs(indivs_gts[indiv1],
                                                   indivs_gts[indiv2])

        # Exit when comparison is below threshold
        if pairs_gds[comparison_key] == ERROR_CODE:
            sys.exit(('Error: less than {0} SNPs shared between {1} '
                      'and {2}').format(MIN_THRESHOLD, indiv1, indiv2))

    # Create header for pairwise distance matrix
    print('\t{0}'.format(SEPARATOR.join(indiv_order)))

    # Output distance matrix row for each individual
    for indiv1 in indiv_order:
        matrix_row = []
        matrix_row.append(indiv1)
        for indiv2 in indiv_order:
            matrix_row.append(get_match_score(pairs_gds, indiv1, indiv2))

        # Output distance matrix row to STDOUT
        print(SEPARATOR.join(matrix_row))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename)
