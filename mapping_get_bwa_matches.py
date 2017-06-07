#!/usr/bin/env python
"""
Extracts a list of succesfully mapped loci from `.sam` file (produced with
`bwa mem`). Successfully mapped loci are identified by default identified as
those with flags 0 and 16 (can be adjusted in MATCH_FLAGS constant), and a 
mapping quality of >=20. Configured for use with single-end reads.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_ID = '@'
COL_QNAME = 0
COL_FLAG = 1
COL_RNAME = 2
COL_POS = 3
COL_MAPQ = 4
COL_SEQ = 9
MATCH_FLAGS = [0, 16]
MIN_MAPQ = 20


def is_header(line):
    """ Return True if line is a header """
    return line[:len(HEADER_ID)] == HEADER_ID


def main(sam_filename):
    # Iterate trough SAM file and identify matching loci
    match_count = ref_count = 0
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if not is_header(line):
            cols = line.split('\t')
            if (int(cols[COL_FLAG]) in MATCH_FLAGS and 
                int(cols[COL_MAPQ]) >= MIN_MAPQ):
                # Output match to file
                print('{0}\t{1}\t{2}\t{3}'.format(cols[COL_QNAME],
                                                  cols[COL_RNAME],
                                                  cols[COL_POS],
                                                  cols[COL_FLAG]))
                match_count += 1
    sam_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sam_filename', metavar='sam_file',
                        help='`bwa mem` output file (`.sam`)')
    args = parser.parse_args()
    main(args.sam_filename)
