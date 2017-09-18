#!/usr/bin/env python
"""
Retains only those loci (/CHROMs) in `.vcf` that are given in file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = '#'


def get_inclusion_list(inclusion_filename):
    """Read inclusion list from text file"""
    inclusion_file = open(inclusion_filename, 'r')
    inclusion_list = []
    for line in inclusion_file:
        locus = line.strip()
        if locus not in inclusion_list:
            inclusion_list.append(locus)
    return inclusion_list


def main(vcf_filename, inclusion_filename):
    # Read inclusion list from text file
    inclusion_list = get_inclusion_list(inclusion_filename)

    # Iterate through VCF and output non-matches
    vcf_file = open(vcf_filename, 'r')
    included_loci = []
    for line in vcf_file:
        # Output headers to output file
        if line[0] == HEADER_CHAR:
            print(line, end='')
        else:
            cols = line.split()
            chrom = cols[0].strip()
            # Output SNP if locus (/CHROM) in inclusion list
            if chrom in inclusion_list:
                print(line, end='')
                included_loci.append(chrom)

    # Assess whether any items in inclusion list were not found in VCF
    loci_not_found = []
    for locus in inclusion_list:
        if locus not in included_loci:
            loci_not_found.append(locus)

    # Output log file if any items in inclusion list were not found
    if len(loci_not_found) > 0:
        output_log = open('{0}_notfound.log'.format(vcf_filename), 'w')
        concat_loci_not_found = ','.join(loci_not_found)
        output_log.write(('Following loci ({0}) were '
                          'not found: {1}').format(len(loci_not_found),
                                                   concat_loci_not_found))
        output_log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('inclusion_filename', metavar='inclusion_file',
                        help='text file with loci to be retained')
    args = parser.parse_args()
    main(args.vcf_filename, args.inclusion_filename)
