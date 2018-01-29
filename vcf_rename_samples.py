#!/usr/bin/env python
"""
Renames sample in `.vcf` file according to list with old/new names; also
outputs samples that are not listed in name_change file.
"""
import argparse
import sys

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_STR = '#CHROM'


def get_dict_of_name_changes(samplenames_filename):
    """ Initialise dictionary with old (keys) and new (values)  """
    rename_file = open(samplenames_filename, 'r')
    name_changes = {}
    for line in rename_file:
        cols = line.replace(',', '\t').split()
        name_changes[cols[0].strip()] = cols[1].strip()
    rename_file.close()
    return name_changes


def main(vcf_filename, samplenames_filename):
    # Initialise dict with new sample names
    name_changes = get_dict_of_name_changes(samplenames_filename)

    # Open logfile for missing samples
    log_file = open('vcf_rename_errorlog.txt', 'w')

    # Open input and output files
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        if line[0:len(HEADER_STR)] == HEADER_STR:
            # Rename each sample in line
            for old_name, new_name in name_changes.items():
                if old_name in line:
                    line = line.replace(old_name, new_name)
                else:
                    log_file.write('Cannot find sample: {0}\n'.format(old_name))
        print(line, end='')

    vcf_file.close()
    log_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('samplenames_filename', metavar='samplenames_file',
                        help='text file (tsv or csv) with old and new name \
                        for each sample (not all samples have to be listed)')
    args = parser.parse_args()
    main(args.vcf_filename, args.samplenames_filename)
