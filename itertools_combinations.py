#!/usr/bin/env python
"""
Generate list with all unique pairwise combinations of values in file. Short
script meant to allow use of itertools.combinations in bash.
"""
import argparse
import itertools

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

SEPARATOR = '\t'


def main(filename):
    # Read order of individuals from file
    with open(filename) as file:
        values = [line.replace(',', ' ').split()[0] for line in file]

    # Output all unique combinations
    for value1, value2 in itertools.combinations(values, 2):
        print('{0}{1}{2}'.format(value1, SEPARATOR, value2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename', metavar='filename',
                        help='input file with values')
    args = parser.parse_args()
    main(args.filename)
