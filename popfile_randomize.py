#!/usr/bin/env python
"""
Pseudo-randomize the assignment of individuals to populations.

Note: with small population sizes, this can lead to very uneven simulated
population sizes. See also the alternative: `popfile_toggleassign.py`.
Individuals in original popfile should be ordered by population.
"""
import sys
import random
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

    # Output random individual assignment
    for indiv in indivs:
        random_pop = random.sample(pops, 1)[0]
        print('{0}\t{1}'.format(indiv, random_pop))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.pop_filename)
