#!/usr/bin/env python
"""
Filter original SNP dataset (`.vcf`) for a particular Fst percentile bin.

Note: order in Fst file needs to correspond with (`.vcf`) file, currently set
(see script CONSTANTS) to work with LOSITAN output file, and output filename
automatically generated from percentile bins.
"""
import sys
import argparse
import numpy as np

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


VCF_FILE_HEADER = '#'
VCF_FILE_EXTENSION = '.vcf'

FST_FILE_HEADER = 'Locus'
FST_FILE_LOCUS_PREFIX = 'SNP_'
FST_FILE_LOCUS_COLUMN_INDEX = 0
FST_FILE_FST_COLUMN_INDEX = 2
FST_FILE_MISSING_VALUE = '-100'


def is_header(line, identifier):
    """ Return True if line is a header """
    return line[:len(identifier)] == identifier


def read_fst_values(fst_filename):
    """ Return dict (key = locus number) and list with FstValues  """
    fst_file = open(fst_filename, 'r')
    fst_values = {}
    fst_values_list = []

    # Iterate through file and store Fst values for each locus
    for line in fst_file:
        if not is_header(line, FST_FILE_HEADER):
            cols = line.replace(',', ' ').split()
            if not cols[FST_FILE_FST_COLUMN_INDEX] == FST_FILE_MISSING_VALUE:
                locus = cols[FST_FILE_LOCUS_COLUMN_INDEX]
                locus_int = int(locus.replace(FST_FILE_LOCUS_PREFIX, ''))
                fst = float(cols[FST_FILE_FST_COLUMN_INDEX])
                fst_values[locus_int] = fst
                fst_values_list.append(fst)
    fst_file.close()
    return fst_values, fst_values_list


def get_output_filename(vcf_filename, min_percentile, max_percentile):
    """ Return output file name """
    output_filename = vcf_filename.lower().replace(VCF_FILE_EXTENSION, '')
    output_filename += '_fst{0}_{1}{2}'.format(min_percentile, max_percentile,
                                               VCF_FILE_EXTENSION)
    return output_filename


def main(vcf_filename, fst_filename, min_percentile, max_percentile):
    # Read Fst values and store in dict/list
    fst_values_dict, fst_values_list = read_fst_values(fst_filename)

    # Calculate Fst percentiles
    min_percentile_fst = np.percentile(fst_values_list,
                                       int(min_percentile))
    max_percentile_fst = np.percentile(fst_values_list,
                                       int(max_percentile))
    print('Fst minimum: {0}'.format(np.amin(fst_values_list)))
    print('Fst maximum: {0}'.format(np.amax(fst_values_list)))
    print('Fst {0}th percentile: {1}'.format(min_percentile,
                                             min_percentile_fst))
    print('Fst {0}th percentile: {1}'.format(max_percentile,
                                             max_percentile_fst))

    # Create output VCF file
    vcf_output_file = open(get_output_filename(vcf_filename, min_percentile,
                                               max_percentile), 'w')

    # Output VCF loci that fall within Fst percentile bin
    locus_index = 0
    retained_snp_count = 0
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        if is_header(line, VCF_FILE_HEADER):
            vcf_output_file.write(line)
        else:
            locus_index += 1
            if locus_index in fst_values_dict:
                fst = fst_values_dict[locus_index]
                if min_percentile_fst <= fst < max_percentile_fst:
                    vcf_output_file.write(line)
                    retained_snp_count += 1

    print('{0} SNPs between {1} and {2} percentile'.format(retained_snp_count,
                                                           min_percentile,
                                                           max_percentile))

    vcf_file.close()
    vcf_output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('fst_filename', metavar='fst_file',
                        help='text file (tsv or csv) with Fst values for \
                        each SNP (same order as in vcf) - currently set \
                        to work with LOSITAN output file')
    parser.add_argument('min_percentile', metavar='min_percentile', type=int,
                        help='min. Fst value for a SNP to be included')
    parser.add_argument('max_percentile', metavar='max_percentile', type=int,
                        help='max. Fst value for a SNP to be included')
    args = parser.parse_args()
    main(args.vcf_filename, args.fst_filename, args.min_percentile,
         args.max_percentile)
