#!/usr/bin/env python
"""
Translates sequential marker numbers back to CHROM/POS from original `.vcf`
file. Several programs only allow for integers to identify markers, this
script is to restore the original CHROM/POS for markers that were identified.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = '#'


def get_numbers_list(numbers_filename):
    """Read list of numbers from text file"""
    numbers_file = open(numbers_filename, 'r')
    numbers_list = []
    for line in numbers_file:
        number = line.strip()
        if number not in numbers_list:
            numbers_list.append(int(number))
    return numbers_list


def main(vcf_filename, numbers_filename):
    # Read list of numbers from text file
    numbers_list = get_numbers_list(numbers_filename)

    # Look up CHROM and POS for each SNP number in list
    vcf_file = open(vcf_filename, 'r')
    snp_counter = 0
    for line in vcf_file:
        # Output headers to output file
        if not line[0] == HEADER_CHAR:
            snp_counter += 1
            if snp_counter in numbers_list:
                cols = line.split()
                chrom = cols[0].strip()
                pos = cols[1].strip()
                print('{0}\t{1}\n'.format(chrom, pos), end='')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('numbers_filename', metavar='markernumbers_file',
                        help='text file with SNP numbers that were identified')
    args = parser.parse_args()
    main(args.vcf_filename, args.numbers_filename)
