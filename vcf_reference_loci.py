#!/usr/bin/env python
"""
Lists all loci (using CHROM column) in `.vcf` that are genotyped for at
least one of the indicated samples/individuals. This can be used to reduce the
dataset to loci matching an included reference (e.g. aposymbiotic) samples.

Note: `vcf` can subsequently be filtered by using the output as inclusion_file
for `vcf_include_chrom.py`.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
CHROM_COLUMN = 0
FIRST_GENOTYPE_COLUMN = 9


def main(vcf_filename, reference_samples):
    vcf_file = open(vcf_filename, 'r')
    genotyped_loci = set()
    all_loci = set()
    individuals = {}

    # Iterate through vcf file and store individual names and genotypes
    for line in vcf_file:
        cols = line.split()
        # Store name of individuals in temporary dictionary
        if line[:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
            for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
                individuals[x] = cols[x]
        # Store genotyped loci for each individual in dictionary set
        elif line[0] != HEADER_CHAR:
            chrom = cols[CHROM_COLUMN]
            all_loci.add(chrom)
            for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
                if individuals[x] in reference_samples and cols[x][0] != '.':
                    genotyped_loci.add(chrom)
    vcf_file.close()

    # Output the list of loci genotyped for the reference samples
    for locus in sorted(genotyped_loci):
        print('{0}'.format(locus))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('reference_samples', metavar='reference_samples',
                        nargs='*', help='sample(s) against which the \
                        remainder of the dataset will be compared')
    args = parser.parse_args()
    main(args.vcf_filename, args.reference_samples)
