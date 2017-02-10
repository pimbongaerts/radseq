#!/usr/bin/env python
"""
Generates artificial crosses between individuals from two indicated (in a
popfile) parentalgroups, and appends crossed individuals to `.vcf` file.
Note: individual SNPs on a single CHROM are independently crossed as if they
are not physically linked - therefore only use when subsampling a single
SNP / CHROM.
"""
import sys
import vcf
import random
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


ALLELES = {0: '00', 1: '01', 2: '11'}
VCF_DIVIDER = '/'
VCF_DELIM = '\t'
HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
CHROM_COL = 0
POS_COL = 1


class CrossedIndividual(object):

    def __init__(self, par1_name, par2_name):
        self.par1_name = par1_name
        self.par2_name = par2_name
        self.gts = {}

    def set_snp_gt(self, chrom, pos, gt):
        if (len(gt) > 3) or (gt[1] not in ('|', '/')):  # Basic verfication
            sys.exit('Unexpected error: incorrect genotype format')
        snp_id = '{0}_{1}'.format(chrom, pos)
        self.gts[snp_id] = gt

    def get_snp_gt(self, chrom, pos):
        snp_id = '{0}_{1}'.format(chrom, pos)
        return self.gts[snp_id]

    def get_parental_name(self):
        parental_name = '{0}x{1}'.format(self.par1_name, self.par2_name)
        return parental_name


def select_parent_combinations(pop_filename, n_crosses):
    """ Determine cross-combinations of parental individuals  """
    # Read individual assignments from file and split in two groups
    pop_file = open(pop_filename, 'r')
    par1_name = par2_name = ''
    par1 = []
    par2 = []
    for line in pop_file:
        cols = line.replace(',', ' ').split()
        indiv = cols[0]
        pop = cols[1]
        if par1_name == '':
            par1_name = pop
            par1.append(indiv)
        elif pop == par1_name:
            par1.append(indiv)
        elif pop == par2_name:
            par2.append(indiv)
        elif par2_name == '':
            par2_name = pop
            par2.append(indiv)
        else:
            sys.exit('Error: more than 2 groups/pops in popfile')
    pop_file.close()

    # Check whether parental groups have min. of required indivs (=n_crosses)
    if n_crosses:
        if min(len(par1), len(par2)) < n_crosses:
            sys.exit('Error: parental pop with insufficient indivs ')
    else:
        n_crosses = min(len(par1), len(par2))

    # Check whether same individual is listed in both parental groups
    if len(set(par1 + par2)) < len(par1 + par2):
        sys.exit('Error: individual(s) listed in both parental pops ')

    # Draw random subsets from each parental group (same amount as n_crosses)
    par1_red = [par1[i] for i in random.sample(range(len(par1)), n_crosses)]
    par2_red = [par1[i] for i in random.sample(range(len(par1)), n_crosses)]

    return par1_red, par2_red


def get_crossed_gt(gt1, gt2):
    """ Create cross between two genotypes """
    if gt1 is not None and gt2 is not None:
        # Only generate genotype when there is no missing parental genotypes
        allele1 = ALLELES[gt1][random.sample(range(2), 1)[0]]
        allele2 = ALLELES[gt2][random.sample(range(2), 1)[0]]
        alleles = sorted([allele1, allele2])
        return '{0}{1}{2}'.format(alleles[0], VCF_DIVIDER, alleles[1])
    else:
        # One of the two genotypes is missing
        return '.{0}.'.format(VCF_DIVIDER)


def generate_crosses(vcf_filename, par1, par2, prefix):
    """ Generate crosses between parental individuals"""
    # Initialize crossed individuals
    crosses = {}
    for i in range(len(par1)):
        cross_id = '{0}{1}'.format(prefix, i)
        crosses[cross_id] = CrossedIndividual(par1[i], par2[i])

    # Generate genotypes for each SNP / individual
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for record in vcf_reader:
        for i in range(len(par1)):
            cross_id = '{0}{1}'.format(prefix, i)
            crossed_gt = get_crossed_gt(record.genotype(par1[i]).gt_type,
                                        record.genotype(par2[i]).gt_type)
            crosses[cross_id].set_snp_gt(record.CHROM, record.POS, crossed_gt)
    return crosses


def append_crosses_to_vcf(vcf_filename, crosses, parentalnames):
    """ Append crosses to vcf """
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        if line[:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
            # Create new header with sample names appended
            if parentalnames:
                cross_names = []
                for cross_id in sorted(crosses.keys()):
                    cross_names.append(crosses[cross_id].get_parental_name())
            else:
                cross_names = sorted(crosses.keys())

            print('{0}{1}{2}'.format(line.strip(), VCF_DELIM,
                                     VCF_DELIM.join(cross_names)))
        elif line[0] != HEADER_CHAR:
            # Output genotypes
            genotypes = []
            for cross_id in sorted(crosses.keys()):
                cols = line.split()
                genotypes.append(crosses[cross_id].get_snp_gt(cols[CHROM_COL],
                                                              cols[POS_COL]))
            print('{0}{1}{2}'.format(line.strip(), VCF_DELIM,
                                     VCF_DELIM.join(genotypes)))
        else:
            print(line, end='')
    vcf_file.close()


def main(vcf_filename, pop_filename, n_crosses, prefix, parentalnames,
         allsnps):
    # Randomly select parental individuals for cross-combinations
    par1, par2 = select_parent_combinations(pop_filename, n_crosses)

    # Create crosses
    crosses = generate_crosses(vcf_filename, par1, par2, prefix)

    # Iterate through all loci in VCF and generate new loci to append
    append_crosses_to_vcf(vcf_filename, crosses, parentalnames)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with the names of the \
                        individuals used for the simulated crosses, and in \
                        the second column which parental population they \
                        belong to (any name can be chosen - as long as there \
                        are exactly two distinct values)')
    parser.add_argument('--n_crosses', '-n', metavar='n_crosses', type=int,
                        help='number of crosses to simulate (should be no \
                        higher than the number of individuals in each of \
                        the two parental populations)')
    parser.add_argument('--prefix', metavar='prefix', default='cross',
                        help='prefix for crosses (used only if \
                        --parentalnames is not set)')
    parser.add_argument('--parentalnames', action='store_true',
                        help='set flag to use names of both parents for cross')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.n_crosses, args.prefix,
         args.parentalnames, args.allsnps)
