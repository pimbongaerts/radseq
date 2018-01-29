#!/usr/bin/env python
"""
Outputs a popfile based on cluster assignment file (from e.g. STRUCTURE) and
outputs a popfile based on those assigments (using supplied assignment
threshold).

Note: I use the formatted CLUMPP output (e.g. `clumpp_K4.out.csv`) from the
`structure_mp` wrapper as assignment file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

FIRST_CLUSTER_COL = 2   # zero-based numbering (0 = first column)
POPNAME_PREFIX = 'CLUSTER'
UNASSIGNED = 'UNASSIGNED'
AMBIGUOUS = 'AMBIGUOUS'


def get_pop_assignments(pop_filename):
    """ Read order of samples from file """
    assigned_pops = {}
    if pop_filename:
        pop_file = open(pop_filename, 'r')
        for line in pop_file:
            cols = line.rstrip().replace(',', ' ').split()
            assigned_pops[cols[0]] = cols[1]
    pop_file.close()
    return assigned_pops


def main(assignment_filename, assign_cut_off, pop_filename):
    # Read pop assignments if popfile given
    if pop_filename:
        assigned_pops = get_pop_assignments(pop_filename)

    # Open csv or tsv file with STRUCTURE outcome
    assignment_file = open(assignment_filename, 'r')

    # Iterate through file and output
    for line in assignment_file:
        cols = line.replace(',', '\t').split('\t')
        assigned_cluster = UNASSIGNED
        for cluster_no, cluster in enumerate(cols[FIRST_CLUSTER_COL:], 1):
            if float(cluster) >= assign_cut_off:
                if assigned_cluster != UNASSIGNED:
                    assigned_cluster = AMBIGUOUS
                else:
                    if pop_filename:
                        assigned_pop = assigned_pops[cols[0]]
                        assigned_cluster = '{0}_{1}'.format(assigned_pop,
                                                            cluster_no)
                    else:
                        assigned_cluster = '{0}_{1}'.format(POPNAME_PREFIX,
                                                            cluster_no)
        print('{0}\t{1}'.format(cols[0], assigned_cluster))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assignment_filename', metavar='assignment_file',
                        help='text file (tsv or csv) with assignment \
                        values for each individual (max. 2 clusters); e.g. a \
                        reformatted STRUCTURE output file')
    parser.add_argument('assign_cut_off', metavar='assign_cut_off', type=float,
                        help='min. assignment value for an individual to be \
                        assigned to a cluster')
    parser.add_argument('-p', '--popfile', dest='pop_filename',
                        metavar='pop_filename',
                        help='optional popfile: use original popnames as \
                        assignment prefix')
    args = parser.parse_args()
    main(args.assignment_filename, args.assign_cut_off, args.pop_filename)
