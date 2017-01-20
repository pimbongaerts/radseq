#!/usr/bin/env python
"""
Removes SNPs from `.vcf` that are above a certain POS value.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = '#'


def main(vcf_filename, highest_pos):
    # Iterate through VCF and output to file
    vcf_file = open(vcf_filename, 'r')
    snp_count = total_count = 0
    for line in vcf_file:
        # Output headers to output file
        if line[0] == HEADER_CHAR:
            print(line, end='')
        else:
            cols = line.split()
            pos = int(cols[1])
            if pos <= int(highest_pos):
                # Output SNP if below threshold
                print(line, end='')
                snp_count += 1
            total_count += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('highest_pos', type=int,
                        help='max. POS value allowed in `.vcf`')
    args = parser.parse_args()
    main(args.vcf_filename, args.highest_pos)
