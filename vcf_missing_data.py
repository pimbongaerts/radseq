#!/usr/bin/env python
"""
Outputs list of missing data (# and % of SNPs) for each sample in VCF, to
identify poor-performing samples to eliminate prior to SNP filtering.

Takes vcf_filename as argument. Outputs to STDOUT (no output file).

"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
MISSING_CHAR = "."
OUTPUT_HEADER = 'INDIVIDUAL\tMISS\tGENO\tTOTAL\t% GENOTYPED'


def main(vcf_filename):
    # Read in genotypes for all individuals
    individuals = {}
    genotypes = {}
    with open(vcf_filename, 'r') as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line[0:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
                cols = line.split('\t')
                # Store individual names with col_index as key
                for col_index, col in enumerate(cols):
                    if col_index >= 9:
                        individual_name = col
                        individuals[col_index] = individual_name
                        genotypes[individual_name] = []
            elif not line[0:len(HEADER_CHAR)] == HEADER_CHAR:
                cols = line.split('\t')
                # Store genotypes for each individual
                for col_index, col in enumerate(cols):
                    if col_index >= 9:
                        genotype = col
                        individual = individuals[col_index]
                        genotypes[individual].append(genotype[0:3])

    # Assess missing data for each individual
    print(OUTPUT_HEADER)
    for individual in sorted(genotypes.keys()):
        # A genotype is counted as missing if either allele is missing
        # (e.g. `./.` or partial calls such as `./0`)
        missing_count = sum(1 for genotype in genotypes[individual]
                            if MISSING_CHAR in genotype)
        total_count = len(genotypes[individual])
        genotyped_count = total_count - missing_count
        if total_count > 0:
            perc_count = round((genotyped_count / total_count) * 100, 2)
        else:
            perc_count = 'NA'
        print('{0}\t{1}\t{2}\t{3}\t{4}'.format(individual, missing_count,
                                               genotyped_count, total_count,
                                               perc_count))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    args = parser.parse_args()
    main(args.vcf_filename)
