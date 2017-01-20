#!/usr/bin/env python
"""
Reduces FASTA file to only those loci listed in supplied text file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(fasta_filename, list_filename):
    # Read list from file and convert to set (unique values only)
    with open(list_filename) as file:
        lines = [line.strip() for line in file]
    loci = set(lines)

    # Parse FASTA and only output those in set
    fasta_file = open(fasta_filename, 'r')
    output_sequence = False
    locuscount = 0

    for line in fasta_file:
        if line[0] == '>':
            # Evaluate if sequence in list of loci
            locusname = line[1:].strip()
            if locusname in loci:
                output_sequence = True
                print(line.strip())
            else:
                output_sequence = False
        elif output_sequence:
            # Output sequence if in list of loci
            print(line.strip())

    fasta_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta_filename', metavar='fasta_file',
                        help='FASTA input file (`.fasta`/ `.fa`)')
    parser.add_argument('inclusion_filename', metavar='inclusion_file',
                        help='text file with names of loci to be included')
    args = parser.parse_args()
    main(args.fasta_filename, args.inclusion_filename)
