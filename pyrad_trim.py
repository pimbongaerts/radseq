#!/usr/bin/env python
"""
Trims sequence length in PyRAD/ipyrad `.alleles` or `.loci` file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def determine_seq_start(line):
    """ Determine the start position of the sequence """
    sample, seq = line.split()
    seq_start = line.find(seq)
    if seq_start < 0:
        sys.exit('Unexpected error')
    else:
        return seq_start


def main(pyrad_filename, seq_length):
    # Open input file
    pyrad_file = open(pyrad_filename, 'r')
    # Iterate through file and trim sequences
    seq_start = 0
    for line in pyrad_file:
        if line[0] == '/':
            # Extract locus info and number
            cols = line.strip().split('|')
            locus_info = cols[0]
            locus_nr = cols[1]
            # Trim locus info
            locus_info_trim = locus_info[0:output_line_length]
            print('{0}|{1}'.format(locus_info_trim, locus_nr))
        else:
            # Determine delimiter if not yet set (first line)
            if seq_start == 0:
                seq_start = determine_seq_start(line)
            # Extract sample name and sequence
            sample, seq = line.strip().split()
            # Trim sequence and determine new line length
            trim_seq = seq[0:int(seq_length)]
            # Output trimmed sequence
            output_line = '{0}{1}'.format(sample.ljust(seq_start), trim_seq)
            print(output_line)
            output_line_length = len(output_line)
    pyrad_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pyrad_filename', metavar='pyrad_file',
                        help='PyRAD allele file (`.loci` or `.allele`)')
    parser.add_argument('seq_length', type=int,
                        help='length to which sequences are trimmed')
    args = parser.parse_args()
    main(args.pyrad_filename, args.seq_length)
