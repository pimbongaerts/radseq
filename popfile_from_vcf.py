#!/usr/bin/env python
"""
Creates tab-separated popfile from `.vcf`, using a subset of the sample name
as population. For example, to use the substring `MGD` from `AFMGD6804H` as
population designation, run script as `python3 popfile_from_vcf vcf_file 3 5`.
"""
import argparse
import sys

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

HEADER_STR = '#CHROM'
SAMPLE_START_COL = 9


def main(vcf_filename, start_pos, end_pos):
    # Open input and output files
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        if line[0:len(HEADER_STR)] == HEADER_STR:
            samples = line.split()[SAMPLE_START_COL:]
            break

    # Output samples and pop designations
    for sample in samples:
            print('{0}\t{1}'.format(sample, sample[start_pos-1:end_pos]))
    vcf_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('start_pos', metavar='start_pos', type=int,
                        choices=range(1, 100),
                        help='character start position in sample name to \
                        be used for population name (one-based)')
    parser.add_argument('end_pos', metavar='end_pos', type=int,
                        choices=range(1, 100),
                        help='character end position in sample name to \
                        be used for population name (one-based)')
    args = parser.parse_args()
    main(args.vcf_filename, args.start_pos, args.end_pos)
