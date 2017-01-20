#!/usr/bin/env python
"""
Concatenates PyRAD/ipyrad sequences from `.loci` file for each individual
and outputs as FASTA (order by popfile).

Note: missing data are filled with gaps (`N`)
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

MISSING_CHAR = 'N'

# TODO: Make provision of popfile optional


def read_pop_file(pop_filename):
    """ Read individuals and population assignments from file """
    pop_file = open(pop_filename, 'r')
    pops_indivs = {}
    for line in pop_file:
        cols = line.replace(',', ' ').split()
        if not cols[1] in pops_indivs:
            pops_indivs[cols[1]] = []
        pops_indivs[cols[1]].append(cols[0])
    pop_file.close()
    return pops_indivs


def main(pyrad_filename, pop_filename):
    # Create dict with pops (keys) and indivs (lists)
    # TODO: Make popfile optional argument
    pops_indivs = read_pop_file(pop_filename)

    # Initialise output dict (sorted by pop)
    output_list = {}
    no_seqs = {}
    for pop in sorted(pops_indivs.keys()):
        for ind in pops_indivs[pop]:
            output_list[ind] = []
            no_seqs[ind] = 0

    # Iterate through pyrad file and output loci that meet criteria
    pyrad_file = open(pyrad_filename, 'r')
    seqs_in_current_locus = {}
    for line in pyrad_file:
        if line[0] == '/':
            # Add all seqs to list when reaching output seq
            for pop in sorted(pops_indivs.keys()):
                for ind in pops_indivs[pop]:
                    if ind in seqs_in_current_locus:
                        output_list[ind].append(seqs_in_current_locus[ind])
                        no_seqs[ind] += 1
                    else:
                        seq = next(iter(seqs_in_current_locus.values()))
                        missing_seq = len(seq) * MISSING_CHAR
                        output_list[ind].append(missing_seq)
            # Empty temporary variable
            seqs_in_current_locus = {}
        else:
            ind, seq = line.split()[:2]
            ind = ind.replace('>', '')  # for old pyrad format
            seqs_in_current_locus[ind] = seq
    pyrad_file.close()

    # Output concatenated FASTA
    for pop in sorted(pops_indivs.keys()):
        for ind in pops_indivs[pop]:
            print('>{0}_{1}_{2}'.format(ind, pop, no_seqs[ind]))
            print(''.join(output_list[ind]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pyrad_filename', metavar='pyrad_file',
                        help='PyRAD file (`.loci`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.pyrad_filename, args.pop_filename)
