#!/usr/bin/env python
"""
Cleans up popfile by eliminating any individuals that are not in `.vcf` file.
"""
import argparse
import sys

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

HEADER_STR = '#CHROM'
SAMPLE_START_COL = 9


def main(vcf_filename, pop_filename):
    # Open vcf file and extract samples
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        if line[0:len(HEADER_STR)] == HEADER_STR:
            samples = line.split()[SAMPLE_START_COL:]
            break
    vcf_file.close()

    # Iterate over popfile and only output those that are extracted from vcf
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        if line.replace(',', ' ').split()[0] in samples:
            print(line, end='')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename)
