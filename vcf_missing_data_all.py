#!/usr/bin/env python
"""
Wrapper that runs `vcf_missing_data.py` on every `.vcf` file found in a
directory and its subdirectories (defaults to the folder containing this
script). By default only the 10 worst-performing samples are shown per file.

Output is printed to STDOUT, grouped per VCF file.
"""
import os
import sys
import glob
import argparse

import vcf_missing_data

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

DEFAULT_LOWEST_N = 10
VCF_PATTERN = '*.vcf'
FILE_HEADER = '### {0}'


def find_vcf_files(directory):
    """Recursively find all `.vcf` files in directory and subdirectories"""
    pattern = os.path.join(directory, '**', VCF_PATTERN)
    return sorted(glob.glob(pattern, recursive=True))


def main(directory, lowest_n, threshold):
    vcf_files = find_vcf_files(directory)
    if not vcf_files:
        sys.exit('No `.vcf` files found in {0}'.format(directory))

    for vcf_filename in vcf_files:
        print(FILE_HEADER.format(vcf_filename))
        vcf_missing_data.main(vcf_filename, lowest_n, threshold)
        print('')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', metavar='directory', nargs='?',
                        default=os.path.dirname(os.path.abspath(__file__)),
                        help='directory to search for `.vcf` files '
                             '(default: folder containing this script)')
    parser.add_argument('-n', type=int, default=None, metavar='N',
                        help='only output the N lowest records (by %% '
                             'genotyped) per file (default: {0}; ignored when '
                             '-t is used unless set '
                             'explicitly)'.format(DEFAULT_LOWEST_N))
    parser.add_argument('-t', type=float, default=None, metavar='THRESHOLD',
                        help='only output names of samples below this %% '
                             'genotyped threshold, one per line and without '
                             'header (e.g. for a vcftools `--remove` file)')
    args = parser.parse_args()
    # Apply the default N only in the normal table view; when building a
    # threshold-based remove list, show all below-threshold samples unless
    # -n was given explicitly
    lowest_n = args.n
    if lowest_n is None and args.t is None:
        lowest_n = DEFAULT_LOWEST_N
    main(args.directory, lowest_n, args.t)
