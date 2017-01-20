#!/usr/bin/env python
"""
Create FASTA file with a representative sequence (using first sample) for
each locus in pyRAD/ipyrad `.loci` or `.allele` file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(input_filename):
    input_file = open(input_filename, 'r')

    current_seqlength = 0
    current_sequence = ''

    for line in input_file:
        if line[0] == '/':
            # Retrieve locus number and write to output file
            cols = line.split('|')
            locus_number = cols[1].strip()
            print('>{0}'.format(locus_number))
            print('{0}'.format(current_sequence))
            current_sequence = ''
        elif current_sequence == '':
            # Store first sequence in locus
            cols = line.split()
            current_sequence = cols[1].strip()

    input_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_filename', metavar='pyrad_file',
                        help='PyRAD allele file (`.loci` or `.allele`)')
    args = parser.parse_args()
    main(args.input_filename)
