#!/usr/bin/env python
"""
Counts SNP occurrence frequency for each POS in `.vcf` file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = '#'
MAX_READ_LENGTH = 1000
MAX_STATUS_BAR_LENGTH = 60


def main(vcf_filename):
    # Initialise dict with SNP positions
    pos_occurrence = dict.fromkeys(range(1, MAX_READ_LENGTH), 0)

    # Open VCF and iterate over SNPs
    vcf_file = open(vcf_filename, 'r')
    highest_count = highest_pos = 0
    for line in vcf_file:
        if line[0] != HEADER_CHAR:
            pos = int(line.split()[1])
            if pos in pos_occurrence:
                pos_occurrence[pos] += 1
                if pos_occurrence[pos] > highest_count:
                    highest_count = pos_occurrence[pos]
                if pos > highest_pos:
                    highest_pos = pos
            else:
                print('Error - unexpected pos: {0}'.format(pos))

    # Initialise file and variables for outputting
    output_file = open(vcf_filename.replace('.vcf', '') + '_occur.txt', 'w')
    group_number = group_count = 0
    last_group_number = 1
    graph_multiplier = MAX_STATUS_BAR_LENGTH / (highest_count * 4)

    # Iterate over occurrence dict and output result
    for pos in pos_occurrence:
        # Output to file
        output_file.write('{0}\t{1}\n'.format(pos, pos_occurrence[pos]))
        # Output to screen every 5 bp
        group_count += pos_occurrence[pos]
        group_number += 1
        graph_bar = '*' * int(group_count * graph_multiplier)
        if group_number == 5:
            print('{0}-{1}\t{2}\t{3}'.format(last_group_number,
                                             pos, group_count, graph_bar))
            last_group_number = pos + 1
            group_number = group_count = 0
        elif pos >= highest_pos:
            print('{0}-{1}\t{2}\t{3}'.format(last_group_number,
                                             pos, group_count, graph_bar))
            break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    args = parser.parse_args()
    main(args.vcf_filename)
