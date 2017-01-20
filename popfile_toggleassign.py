#!/usr/bin/env python
"""
Shuffle the assignment of individuals to populations by assigning the indivs
sequentially to the different pops. The assignment is not completely random,
but does generate equal population sizes (which otherwise differ substantially
when using random assignment under originally small population sizes).
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(pop_filename):
    # Read individuals and population assignments
    pop_file = open(pop_filename, 'r')
    indivs = []
    pops = set()
    for line in pop_file:
        cols = line.replace(',', ' ').split()
        indivs.append(cols[0])
        pops.add(cols[1])
    pop_file.close()

    # Create new assignment by multiplying set
    pops_ordered = list(pops) * ((len(indivs) // len(pops)) + 1)

    # Output new individual assignments
    for index, indiv in enumerate(indivs):
        new_assign_pop = pops_ordered[index]
        print('{0}\t{1}'.format(indiv, new_assign_pop))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.pop_filename)
