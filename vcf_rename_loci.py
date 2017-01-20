#!/usr/bin/env python
"""
Renames CHROMS in `.vcf` file according to list with old/new names, and only
outputs those loci that are listed.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = '#'

# TO DO: Add flag that will allow for output of non-listed loci


def get_dict_of_name_changes(locusnames_filename):
    """ Initialise dictionary with old (keys) and new (values)  """
    locusnames_file = open(locusnames_filename, 'r')
    name_changes = {}
    for line in locusnames_file:
        cols = line.replace(',', ' ').split()
        name_changes[cols[0].strip()] = cols[1].strip()
    locusnames_file.close()
    return name_changes


def main(vcf_filename, locusnames_filename):
    # Initialise dictionary with old (keys) and new (values)
    name_changes = get_dict_of_name_changes(locusnames_filename)

    # Open input and output files
    vcf_file = open(vcf_filename, 'r')
    # Iterate through VCF and output to file
    snp_count = total_count = 0
    for line in vcf_file:
        # Output headers to output file
        if line[0] == HEADER_CHAR:
            print(line, end='')
        else:
            cols = line.split()
            old_chrom = cols[0].strip()
            if old_chrom in name_changes:
                new_chrom = name_changes[old_chrom]
                unchanged_cols = '\t'.join(cols[1:])
                new_line = '{0}\t{1}'.format(new_chrom, unchanged_cols)
                print(new_line)
                snp_count += 1
            total_count += 1
    vcf_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('locusnames_filename', metavar='locusnames_file',
                        help='text file (tsv or csv) with old and new name \
                        for each locus (/CHROM)')
    args = parser.parse_args()
    main(args.vcf_filename, args.locusnames_filename)
