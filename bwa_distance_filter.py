#!/usr/bin/env python
"""
Filters list of loci mapped to genome scaffolds, so that they are spaced at
least a certain distance. Input file should be tab-separated with columns in
the following order (no header): rad_locus, ref_scaffold, ref_start_pos, flag.
Required spacing between POS will be spacing + max_locus_length.
"""
import sys
import argparse
import numpy as np

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

INPUT_FILE_DELIM = '\t'
STR_LENGTH = 255
C_RADLOCUS = 'rad_locus'
C_REFSCAFFOLD = 'ref_scaffold'
C_REFSTARTPOS = 'ref_start_pos'
C_FLAG = "flag"

MAPPING_RESULTS_DTYPES = [(C_RADLOCUS, np.str_, STR_LENGTH),
                          (C_REFSCAFFOLD, np.str_, STR_LENGTH),
                          (C_REFSTARTPOS, np.int),
                          (C_FLAG, np.int)]


def load_mapping_results(filename):
    """ Load mapping results into numpy array and sort """
    mapping_results = np.genfromtxt(filename,
                                    delimiter=INPUT_FILE_DELIM,
                                    dtype=MAPPING_RESULTS_DTYPES)
    mapping_results.sort(order=(C_REFSCAFFOLD, C_REFSTARTPOS))
    return mapping_results

def load_locus_names(loci_filename):
    """ Read list with locus names to be considered (others skipped) """
    loci_to_consider = []
    if loci_filename:
        loci_file = open(loci_filename, 'r')
        for line in loci_file:
            cols = line.rstrip().replace(',', ' ').split()
            loci_to_consider.append(cols[0])
    return loci_to_consider

def output_spaced_loci(mapping_results, required_spacing, loci_to_consider):
    """ Output list of loci that satisfy spacing threshold """
    prev_scaffold = ""
    prev_pos = 0
    for row in np.nditer(mapping_results):
        if (len(loci_to_consider) > 0 and 
           row[C_RADLOCUS] not in loci_to_consider):
            continue

        if row[C_REFSCAFFOLD] == prev_scaffold:
            pos_spacing = row[C_REFSTARTPOS] - prev_pos
        else:
            # If first locus in scaffold then take ref_start_pos, to ensure
            # that ref_start_pos > required_spacing away from beginning
            pos_spacing = row[C_REFSTARTPOS]

        if pos_spacing >= required_spacing:
            # Include this locus
            # print("INCLUDE: ", row, pos_spacing)
            prev_scaffold = row[C_REFSCAFFOLD]
            prev_pos = row[C_REFSTARTPOS]
            print('{0}\t{1}\t{2}\t{3}\t#INCLUDE:{4}'.format(row[C_RADLOCUS],
                                                            row[C_REFSCAFFOLD],
                                                            row[C_REFSTARTPOS],
                                                            row[C_FLAG],
                                                            pos_spacing))
        else:
            # Skip this locus (and compare next to previous one)
            print('{0}\t{1}\t{2}\t{3}\t#EXCLUDE:{4}'.format(row[C_RADLOCUS],
                                                            row[C_REFSCAFFOLD],
                                                            row[C_REFSTARTPOS],
                                                            row[C_FLAG],
                                                            pos_spacing))


def main(filename, spacing, max_locus_length, loci_filename):
    mapping_results = load_mapping_results(filename)
    loci_to_consider = load_locus_names(loci_filename)
    required_spacing = spacing + max_locus_length
    output_spaced_loci(mapping_results, required_spacing, loci_to_consider)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename', metavar='filename',
                        help='input file')
    parser.add_argument('spacing', metavar='spacing',
                        help='desired spacing between loci')
    parser.add_argument('max_locus_length', metavar='max_locus_length',
                        help='max. length of loci')
    parser.add_argument('-l', '--loci', dest='loci_filename', 
    					metavar='loci_file',
                        help='file with loci to be considered')
    args = parser.parse_args()
    main(args.filename, int(args.spacing), int(args.max_locus_length), 
         args.loci_filename)
