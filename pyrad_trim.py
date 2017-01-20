#!/usr/bin/env python
"""
Trims sequence length in PyRAD/ipyrad `.alleles` or `.loci` file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(pyrad_filename, seq_length):
    # Open input file
    pyrad_file = open(pyrad_filename, 'r')
    # Iterate through file and trim sequences
    output_line_length = 0
    for line in pyrad_file:
        if line[0] == '/':
            # Extract locus info and number
            cols = line.split('|')
            locus_info = cols[0]
            locus_nr = cols[1].strip()
            # Trim locus info
            locus_info_trim = locus_info[0:output_line_length + 1]
            print('{0}|{1}'.format(locus_info_trim, locus_nr))
        else:
            # Extract sample name and sequence
            cols = line.split()
            sample = cols[0]
            seq = cols[1].strip()
            # Trim sequence and determine new line length
            trim_seq = cols[1][0:int(seq_length)]
            # Output trimmed sequence
            output_line = '{0}\t{1}'.format(sample, trim_seq)
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
