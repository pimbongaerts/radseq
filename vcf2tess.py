#!/usr/bin/env python
"""
Converts `.vcf` file to TESS input files (genotypes and coordinates).
Requires a popfile and a file with coordinates for each population (decimal
]degrees), a then simulates individual coordinates by adding a certain
amount of noise.
Note: outputs individuals in the same order as popfile.
"""
import sys
import argparse
import vcf
import random

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

HOMOZYGOTE_REF = '0'
HETEROZYGOTE = '1'
HOMOZYGOTE_ALT = '2'
MISSING = 'NA'
OUTPUT_DELIM = '\t'


def get_indivs_and_pops(pop_filename):
    """ Read indvs and population assignments from popfile """
    indivs = []         # Stored to maintain popfile order
    indivs_pops = {}    # Stores individual assignments to pops
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        cols = line.replace(',', ' ').split()
        indivs.append(cols[0])
        indivs_pops[cols[0]] = cols[1]
    pop_file.close()
    return indivs, indivs_pops


def get_tess_genotype(vcf_genotype):
    """ Convert VCF to outflank genotype """
    if not vcf_genotype:
        return MISSING
    elif vcf_genotype[0] == '.' and vcf_genotype[2] == '.':
        return MISSING
    elif vcf_genotype[0] == '0' and vcf_genotype[2] == '0':
        return HOMOZYGOTE_REF
    elif vcf_genotype[0] == '0' and vcf_genotype[2] == '1':
        return HETEROZYGOTE
    elif vcf_genotype[0] == '1' and vcf_genotype[2] == '0':
        return HETEROZYGOTE
    elif vcf_genotype[0] == '1' and vcf_genotype[2] == '1':
        return HOMOZYGOTE_ALT
    else:
        sys.exit('Error: unexpected genotype: `{0}`'.format(vcf_genotype))


def get_coords(coord_filename, indivs_pops, noise):
    """ Create coordinates for each indiv from coordfile """
    pop_coords = {}
    coord_file = open(coord_filename, 'r')
    for line in coord_file:
        cols = line.replace(',', ' ').split()
        pop_coords[cols[0]] = [float(cols[1]), float(cols[2])]
    coord_file.close()

    # Create dict with coordinates for each individual
    ind_coords = {}
    for indiv in indivs_pops:
        pop = indivs_pops[indiv]
        if pop in pop_coords:
            noise_gen = random.uniform(noise * 0.1, noise)
            ind_coords[indiv] = [pop_coords[pop][0] + noise_gen,
                                 pop_coords[pop][1] + noise_gen]
        else:
            sys.exit('Error: population `{0}` not in coordfile'.format(pop))
    return ind_coords


def output_coords(output_prefix, indivs, ind_coords):
    """ Output coords for each individual """
    output_coord_file = open('{0}_coords.txt'.format(output_prefix), 'w')
    for indiv in indivs:
        lon, lat = ind_coords[indiv]
        output_coord_file.write('{1}{0}{2}\n'.format(OUTPUT_DELIM, lon, lat))
    output_coord_file.close()


def get_genotypes(vcf_filename, indivs):
    """ Get genotypes for all individuals """
    # Initialize dict for genotypes
    indiv_genotypes = {}
    for indiv in indivs:
        indiv_genotypes[indiv] = []
    # Retrieve genotypes for all indivs
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for record in vcf_reader:
        for indiv in indivs:
            vcf_genotype = record.genotype(indiv)['GT']
            tess_genotype = get_tess_genotype(vcf_genotype)
            indiv_genotypes[indiv].append(tess_genotype)
    return indiv_genotypes


def output_genotypes(output_prefix, indivs, indiv_genotypes):
    """ Output individuals (rows) and genotypes (cols) """
    output_genos_file = open('{0}_tessgenos.txt'.format(output_prefix), 'w')
    for indiv in indivs:
        concat_genotypes = OUTPUT_DELIM.join(indiv_genotypes[indiv])
        output_genos_file.write('{0}\n'.format(concat_genotypes))
    output_genos_file.close()


def main(vcf_filename, pop_filename, coord_filename, noise, output_prefix):
    # Read order of indivs and population assignments
    indivs, indivs_pops = get_indivs_and_pops(pop_filename)

    # Get coords for each individual
    ind_coords = get_coords(coord_filename, indivs_pops, noise)

    # Output simulated coords
    output_coords(output_prefix, indivs, ind_coords)

    # Get genotypes for all individuals in popfile
    indiv_genotypes = get_genotypes(vcf_filename, indivs)

    # Output genotypes
    output_genotypes(output_prefix, indivs, indiv_genotypes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    parser.add_argument('coord_filename', metavar='coord_file',
                        help='text file (tsv or csv) with populations and \
                        their lats and longs (in decimal degrees)')
    parser.add_argument('--noise', '-n', metavar='noise', type=float,
                        default=1e-10,
                        help='max. amount of noise to be added (default = \
                        1e-10)')
    parser.add_argument('output_prefix', metavar='output_prefix',
                        help='name prefix for output files')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.coord_filename,
         args.noise, args.output_prefix)
