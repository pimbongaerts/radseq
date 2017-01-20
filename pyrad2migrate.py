#!/usr/bin/env python
"""
Converts PyRAD `.allele` file to migrate-n input file (population designated
indicated in supplied popfile).

Note: only appropriate for PyRAD `.allele` file (not `.loci`).
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def create_pop_dictionary(pop_filename):
    """ Create dicts of individuals and populations """
    pop_file = open(pop_filename, 'r')
    pops_indivs = {}
    indivs_pops = {}

    for line in pop_file:
        line = line.rstrip()
        cols = line.replace(',', '\t').split()
        pops_indivs.setdefault(str(cols[1]), [])
        pops_indivs[str(cols[1])].append(str(cols[0]))
        indivs_pops[str(cols[0])] = str(cols[1])

    return pops_indivs, indivs_pops


def get_individual_name(line):
    """ Get individual name from sequence line in .alleles file """
    individual = line.split()[0][1:].split('_')[0].strip()
    return individual


def main(input_filename, pop_filename):
    # Create dicts of indivs and pops
    pops_indivs, indivs_pops = create_pop_dictionary(pop_filename)

    # Initialise variables and dictionaries
    number_of_loci = 0
    sequence = ''
    population_output = {}
    locus_lengths = {}
    pop_lines = {}
    pop_locus_hap = {}
    for pop in pops_indivs:
        pop_lines.setdefault(pop, [])
        pop_locus_hap.setdefault(pop, {})

    # Iterate through .allele file
    input_file = open(input_filename, 'r')
    for line in input_file:
        # Parse when reaching locus definition line
        if line[0] == '/':
            # Add information to dict when reaching line with locus definition
            number_of_loci += 1
            cols = line.split('|')
            locusname = cols[1].strip()
            for pop in pops_indivs:
                pop_locus_hap[pop].setdefault(locusname, [])
                pop_locus_hap[pop][locusname] += pop_lines[pop]

            # Store sequence length for locus
            locus_lengths[locusname] = len(sequence)
            sequence = ''

            # Re-initialise locus-population dictionary
            pop_lines = {}
            for population in pops_indivs:
                pop_lines.setdefault(population, [])

        else:
            # Store line item in list for corresponding population
            individual = get_individual_name(line)
            if individual in indivs_pops:
                sequence = line.split()[1]
                pop = indivs_pops[individual]
                pop_lines[pop].append(line[1:].strip())
    input_file.close()

    # Generate header for migrate file
    number_of_populations = len(pops_indivs)
    print('{0} {1} {2}'.format(number_of_populations,
                               number_of_loci, input_filename))
    population = next(iter(pops_indivs.keys()))    # Select 'first' population
    locus_lengths_sorted = []
    for locus in sorted(pop_locus_hap[population]):
        locus_lengths_sorted.append(str(locus_lengths[locus]))
    locus_lengths_concatenated = ' '.join(locus_lengths_sorted)
    print('{0}'.format(locus_lengths_concatenated))

    # Output populations with sequences for migrate file
    for population in sorted(pops_indivs):
        # Output locus counts
        counts = []
        for locus in sorted(pop_locus_hap[population]):
            counts.append(str(len(pop_locus_hap[population][locus])))
        counts_concatenated = ' '.join(counts)
        print('{0} {1}'.format(counts_concatenated, population))
        # Output actual haplotypes for each locus
        for locus in sorted(pop_locus_hap[population]):
            for line in pop_locus_hap[population][locus]:
                print('{0}'.format(line))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_filename', metavar='allele_file',
                        help='PyRAD allele file (.allele)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    args = parser.parse_args()
    main(args.input_filename, args.pop_filename)
