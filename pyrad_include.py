#!/usr/bin/env python
"""
Reduces (i)pyrad file to only those samples listed in supplied text file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(loci_filename, list_filename, min_samples):
    # Read list from file and convert to set (unique values only)
    with open(list_filename) as file:
        lines = [line.strip() for line in file]
    samples = set(lines)

    # Parse FASTA and only output those in set
    loci_file = open(loci_filename, 'r')
    seqs_to_output = []

    for line in loci_file:
        if line[0] == '>':
            # Evaluate if sample in list of samples
            if line[1:].split()[0] in samples:
                seqs_to_output.append(line.strip())
        elif line[0] == '/' and len(seqs_to_output) >= int(min_samples):
            # Output locus info if sequences were outputted
            for seq in seqs_to_output:
                print(seq)
            print(line.strip())
            seqs_to_output = []
        else:
            seqs_to_output = []
    loci_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('loci_filename', metavar='loci_file',
                        help='samples input file')
    parser.add_argument('inclusion_filename', metavar='inclusion_file',
                        help='text file with names of samples to be included')
    parser.add_argument('min_samples', metavar='min_samples',
                        help='minimum number of samples for locus to be included')
    args = parser.parse_args()
    main(args.loci_filename, args.inclusion_filename, args.min_samples)
