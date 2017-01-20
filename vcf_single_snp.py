#!/usr/bin/env python
"""
Reduces `.vcf` file to a single 'random' SNP per locus/chrom. Use for analyses
that require SNPs that are not physically linked. (although note that they of
course still may be - particularly so when dealing with short loci)
"""
import sys
import argparse
import random

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"


def main(vcf_filename):
    input_file = open(vcf_filename, 'r')

    total_count = chrom_count = 0
    previous_chrom = ''
    lines_for_chrom = []

    for line in input_file:
        # NOTE: if only wanting to use a substring of the chrom name, split on
        # different delimiting character, e.g.: cols = line.split('_')
        cols = line.split()
        chrom = cols[0]

        # Output headers to output file
        if line[0] == HEADER_CHAR:
            print('{0}'.format(line), end='')
        else:
            # Output random SNP when reaching new chrom
            if len(lines_for_chrom) > 0 and chrom != previous_chrom:
                random_line_no = random.randint(0, len(lines_for_chrom) - 1)
                random_line = lines_for_chrom[random_line_no]
                print('{0}'.format(random_line), end='')
                chrom_count += 1
                lines_for_chrom = []

            lines_for_chrom.append(line)
            total_count += 1
        previous_chrom = chrom

    input_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    args = parser.parse_args()
    main(args.vcf_filename)
