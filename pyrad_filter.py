#!/usr/bin/env python
"""
Filters PyRAD output file (`.loci`) for those loci (1) present or absent (using
--exclude flag) in supplied list, (2) genotyped for at least X number of
samples, and (3) with at least Y number of informative sites.

Note: can also be used for `.alleles` file but then 2x the number of samples
should be given (assuming diploid individual).
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def meets_matching_criterium(current_locus_id, loci_to_be_included, exclude):
    """ Evaluate if locus meets matching criterium (1) """
    if not exclude:     # Include loci present in list
        return current_locus_id in loci_to_be_included
    elif exclude:       # Include loci absent in list (-i flag)
        return current_locus_id not in loci_to_be_included


def locus_meets_ind_threshold(seqs_in_current_locus, ind_threshold):
    """ Evaluate if locus meets threshold """
    return len(seqs_in_current_locus) >= ind_threshold


def locus_meets_snp_threshold(info_current_locus, snp_threshold):
    """ Evaluate if locus meets SNP threshold """
    return info_current_locus.count('*') >= snp_threshold


def main(pyrad_filename, list_filename, ind_threshold, snp_threshold,
         exclude):

    # Read list of loci to be included and convert to set (unique vals only)
    with open(list_filename) as file:
        lines = [int(line.strip()) for line in file]
    loci_to_be_included = set(lines)

    # Iterate through pyrad file and output loci that meet criteria
    pyrad_file = open(pyrad_filename, 'r')
    seqs_in_current_locus = []
    info_current_locus = ''
    output_count = 0
    for line in pyrad_file:
        if line[0] == '/':
            # Evaluate if locus meets 3 criteria
            info_current_locus = line.strip()
            current_locus_id = int(info_current_locus.split('|')[1])
            if meets_matching_criterium(current_locus_id, loci_to_be_included,
                                        exclude) and \
                locus_meets_ind_threshold(seqs_in_current_locus,
                                          ind_threshold) and \
                locus_meets_snp_threshold(info_current_locus,
                                          snp_threshold):
                # Output to file if meets all criteria
                for output_line in seqs_in_current_locus:
                    print(output_line)
                print(info_current_locus)
                output_count += 1
            # Empty temporary variables
            seqs_in_current_locus = []
            info_current_locus = ''
        else:
            # Store sequences for each locus in list
            seqs_in_current_locus.append(line.strip())
    pyrad_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pyrad_filename', metavar='pyrad_file',
                        help='PyRAD file (`.loci`)')
    parser.add_argument('list_filename', metavar='loci_file',
                        help='text file with PyRAD loci to be included')
    parser.add_argument('ind_threshold', metavar='sample_threshold', type=int,
                        help='min. number of samples genotyped \
                              for a locus to be included')
    parser.add_argument('snp_threshold', metavar='snp_threshold', type=int,
                        help='min. number of SNPs for a locus to be included')
    parser.add_argument('-e', '--exclude', action='store_true',
                        help='use the loci in loci_file as exclusion list')
    args = parser.parse_args()
    main(args.pyrad_filename, args.list_filename,
         int(args.ind_threshold), int(args.snp_threshold), args.exclude)
