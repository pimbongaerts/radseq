#!/usr/bin/env python
"""
Outputs number of reads for each fastq.gz sample to text file, and prints
mean/min/max to STDOUT.

Note: Does not work with all FASTQ formats, and correct depending on OS
`zcat` or `gzcat` needs to be set in COMPRESS_UTIL constant.
"""
import os
import sys
import argparse
import numpy as np
from subprocess import Popen, PIPE

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

COMPRESS_UTIL = 'gzcat'      # use zcat or gzcat


def main(path, output_filename):
    # Initialize list to store read counts and open output file
    read_counts = []
    output_file = open(output_filename, 'w')

    # Iterate over all files and count number of sequences
    for filename in os.listdir(path):
        if filename.endswith('.fastq.gz'):
            # Get linecount through bash
            cmd = '{0} "{1}/{2}" | wc -l'.format(COMPRESS_UTIL, path, filename)
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            stdout, stderr = p.communicate()
            # Calculate
            numseqs = int(int(stdout) / 4)
            read_counts.append(numseqs)
            sample_name = filename.replace('.fastq.gz', '')
            output_file.write('{0}\t{1} reads\n'.format(sample_name, numseqs))
    output_file.close()

    # Output mean
    print('mean: {:,} reads'.format(int(np.mean(read_counts))))
    print('min: {:,} reads'.format(np.min(read_counts)))
    print('max: {:,} reads'.format(np.max(read_counts)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', metavar='path',
                        help='path (for current directory use `.`)')
    parser.add_argument('output_filename', metavar='output_filename',
                        help='name of text output file')
    args = parser.parse_args()
    main(args.path, args.output_filename)
