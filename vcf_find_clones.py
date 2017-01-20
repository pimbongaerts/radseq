#!/usr/bin/env python
"""
Script compares the allelic similarity of individuals in a VCF, and outputs
all pairwise comparisons. This can be used to detect potential clones
based on percentage match.

Note: highest matches can be assessed in the output file by using `$ sort -rn
--key=5 output_file.txt | head -n 50` in the terminal.
"""
import sys
import argparse
import operator
import itertools
import math

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
FIRST_GENOTYPE_COLUMN = 9
OUTPUT_HEADER = 'Ind1,Ind2,SNPs,Match_score,perc_match\n'


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


def get_match_score(genotype1, genotype2):
    """ Get matching score of two individuals """
    if genotype1 == genotype2:
        match_score = 1
    elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2] or \
            genotype1[0] == genotype2[2] or genotype1[2] == genotype2[0]:
        match_score = 0.5
    else:
        match_score = 0
    return match_score


def main(vcf_filename):
    # Get genotypes from text file
    genotypes = get_genotypes_from_vcf(vcf_filename)

    # Pairwise comparison of each combination of individuals
    progress_count = last_progress_update = 0
    unique_combinations = (math.pow(len(genotypes), 2) - len(genotypes)) / 2
    for individual1, individual2 in itertools.combinations(genotypes, 2):
        genotypes1 = genotypes[individual1]
        genotypes2 = genotypes[individual2]

        # Compare genotypes of both individuals for each SNP
        match_count = total_count = 0
        for x in range(0, len(genotypes1)):
            genotype1 = genotypes1[x][:3]
            genotype2 = genotypes2[x][:3]
            if genotype1[0] != '.' and genotype2[0] != '.':
                match_score = get_match_score(genotype1, genotype2)
                match_count += match_score
                total_count += 1

        # Output allelic identity match
        if total_count > 0:
            perc_match = round((match_count / total_count) * 100, 2)
            print('{0}\t{1}\t{2}\t{3}\t{4}'.format(individual1,
                                                   individual2,
                                                   total_count,
                                                   match_count,
                                                   perc_match))
        else:
            print('{0}\t{1}\t{2}\t{3}\t{4}'.format(individual1,
                                                   individual2,
                                                   total_count,
                                                   match_count,
                                                   "NA"))
        progress_count += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    args = parser.parse_args()
    main(args.vcf_filename)
