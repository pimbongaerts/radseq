#!/usr/bin/env python
"""
Calculates the mean for every data-cell in a csv across multiple files.
Note: there should be no other `,` other than as a delimiter.
"""
import os
import glob
import argparse
import random
import numpy

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

DELIM_CHAR = ','

def get_csv_filenames_in_path(path):
    """ Return list with all csv files in path """
    return sorted([fname for fname in glob.glob(os.path.join(path, '*.csv'))])

def main(csv_path, start_column):
    csv_filenames = get_csv_filenames_in_path(csv_path)
    csv_data = {}

    # Store contents of all csv_files in lists/memory
    for csv_filename in csv_filenames:
        csv_data[csv_filename] = open(csv_filename, 'r').read().splitlines()

    # Output calculated mean values to STDOUT (using metadata of first file)
    first_file = csv_data[csv_filenames[0]]
    last_col = len(first_file[0].split(DELIM_CHAR))
    for line_nr, line in enumerate(csv_data[csv_filenames[0]]):
        # Calculate means
        col_means = []
        for column in range((start_column-1), last_col):
            temp_values = []
            for csv_filename in csv_filenames:
                temp_line = csv_data[csv_filename][line_nr]
                temp_values.append(float(temp_line.split(DELIM_CHAR)[column]))
            col_means.append('{:.4f}'.format(numpy.mean(temp_values)))
        # Output new data
        if start_column > 1:
            metadata = line.split(DELIM_CHAR)[0:start_column-1]
            print('{0}{1}{2}'.format(DELIM_CHAR.join(metadata), 
                                     DELIM_CHAR, 
                                     DELIM_CHAR.join(col_means)))
        else:
            print(DELIM_CHAR.join(col_means))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv_path', metavar='csv_path',
                        help='input path with data files (`.csv`), all \
                              containing data in the exact same order.')
    parser.add_argument('start_column', metavar='start_column',
                        help='first column (1-based) that contains data for \
                              which the mean needs to be calculated across \
                              files. The columns before that are included as \
                              metadata and are assumed to be identical across \
                              files')
    args = parser.parse_args()
    main(args.csv_path, int(args.start_column))