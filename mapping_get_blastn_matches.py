#!/usr/bin/env python
"""
Extracts a list of loci that have a blastn e-value below a certain threshold,
and outputs the (first) matching reference locus, as well as the alignment
length, nident, and pident. Results are outputted to file with the chosen
e-value as post-fix, and STDOUT gives minimum alignment stats for filtered
loci.

Note: fields in input file should be (in this order): query id, subject id,
alignment length, identity, perc. identity, evalue, bitscore (additional fields
beyond that are fine). This can be achieved by using `blastn` with the
following argument: `-outfmt 7 qseqid sseqid length nident pident evalue
bitscore`.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

COMMENT_CHARACTER = '#'
COL_QSEQID = 0
COL_SSEQID = 1
COL_LENGTH = 2
COL_NIDENT = 3
COL_PIDENT = 4
COL_EVALUE = 5


def is_comment(line):
    """ Return True if line is a comment """
    return line[:len(COMMENT_CHARACTER)] == COMMENT_CHARACTER


def main(input_filename, max_evalue):
    # Open output file
    output_postfix = '_match{0}.txt'.format(max_evalue)
    output_filename = '{0}{1}'.format(input_filename.replace('.txt', ''),
                                      output_postfix)
    output_filename = open(output_filename, 'w')

    # Initialise variables
    min_length = min_nident = min_pident = match_count = 0

    # Iterate through blastn file and extract matches
    input_file = open(input_filename, 'r')
    for line in input_file:
        if is_comment(line):
            match_flag = False
        elif not match_flag:
            cols = line.split('\t')
            if float(cols[COL_EVALUE]) <= float(max_evalue):
                # Write to output file
                output_filename.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    cols[COL_QSEQID], cols[COL_SSEQID], cols[COL_LENGTH],
                    cols[COL_NIDENT], cols[COL_EVALUE]))
                # Keep track of min. length, nident and pident values
                if min_length == 0 or float(cols[COL_LENGTH]) < min_length:
                    min_length = float(cols[COL_LENGTH])
                if min_nident == 0 or float(cols[COL_NIDENT]) < min_nident:
                    min_nident = float(cols[COL_NIDENT])
                if min_pident == 0 or float(cols[COL_PIDENT]) < min_pident:
                    min_pident = float(cols[COL_PIDENT])
                #
                match_count += 1
                match_flag = True
    # Close input and output files
    input_file.close()
    output_filename.close()

    # Output match summary
    print(('Matches: {0} | Min.length: {1} bp | Min. nident: {2} bp | '
           'Min. pident: {3} %').format(match_count, min_length,
                                        min_nident, min_pident))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('blastn_filename', metavar='blastn_file',
                        help='blastn output file with the following fields \
                        (in that order): query id, subject id, alignment \
                        length, identity, perc. identity, evalue, bitscore')
    parser.add_argument('evalue_cut_off', metavar='evalue_cut_off', type=float,
                        help='maximum e-value for match to be included')
    args = parser.parse_args()
    main(args.blastn_filename, args.evalue_cut_off)
