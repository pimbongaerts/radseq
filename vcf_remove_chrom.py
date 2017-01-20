#!/usr/bin/env python
"""
Excludes those loci (/CHROMs) in `.vcf` that are listed in exclusion list.
Also outputs a logfile with loci that were listed but not present in `.vcf`.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = '#'


def get_exclusion_list(exclusion_filename):
    """Read exclusion list from text file"""
    exclusion_file = open(exclusion_filename, 'r')
    exclusion_list = []
    for line in exclusion_file:
        locus = line.strip()
        if locus not in exclusion_list:
            exclusion_list.append(locus)
    return exclusion_list


def main(vcf_filename, exclusion_filename):
    # Read exclusion list from text file
    exclusion_list = get_exclusion_list(exclusion_filename)

    # Iterate through VCF and output non-matches
    vcf_file = open(vcf_filename, 'r')
    excluded_loci = []
    for line in vcf_file:
        # Output headers to output file
        if line[0] == HEADER_CHAR:
            print(line, end='')
        else:
            cols = line.split()
            chrom = cols[0].strip()
            # Output SNP if locus (/CHROM) not in exclusion list
            if chrom not in exclusion_list:
                print(line, end='')
            else:
                excluded_loci.append(chrom)

    # Assess whether any items in exclusion list were not found in VCF
    loci_not_found = []
    for locus in exclusion_list:
        if locus not in excluded_loci:
            loci_not_found.append(locus)

    # Output log file if any items in exclusion list were not found
    if len(loci_not_found) > 0:
        output_log = open('{0}_removed.log'.format(vcf_filename), 'w')
        concat_loci_not_found = ','.join(loci_not_found)
        output_log.write(('Following loci ({0}) were '
                          'not found: {1}').format(len(loci_not_found),
                                                   concat_loci_not_found))
        output_log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('exclusion_filename', metavar='exclusion_file',
                        help='text file loci to be excluded')
    args = parser.parse_args()
    main(args.vcf_filename, args.exclusion_filename)
