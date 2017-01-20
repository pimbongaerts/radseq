#!/usr/bin/env python
"""
Outputs genotype frequencies for specific SNPs in each population, organised
by group.
"""
import sys
import vcf
import numpy
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

COL_INDIV = 0
COL_POP = 1
COL_GROUP = 2
MISSING_GENOTYPE = '.|.'
CHROM_POS_SEPARATOR = '_'
NAME_HOM_REF = 'hom_ref'
NAME_HET = 'het'
NAME_ALT_REF = 'alt_ref'


def create_dicts_from_factor_file(factor_filename):
    """
    Create dictionaries from factor_file: groups (keys) with populations
    (values), and populations (keys) with individuals (values)
    """
    factor_file = open(factor_filename, "r")
    populations_individuals = {}
    groups_populations = {}

    for line in factor_file:
        cols = line.replace(',', ' ').rstrip().split()
        # Add population to group dictionary (if not already added)
        groups_populations.setdefault(cols[COL_GROUP], [])
        if cols[COL_POP] not in groups_populations[cols[COL_GROUP]]:
            groups_populations[cols[COL_GROUP]].append(cols[COL_POP])
        # Add individuals to populations
        populations_individuals.setdefault(cols[COL_POP], [])
        populations_individuals[cols[COL_POP]].append(cols[COL_INDIV])

    return groups_populations, populations_individuals


def create_SNP_inclusion_list(SNP_inclusion_filename):
    """
    Create list with all SNPs to be included
    """
    inclusion_file = open(SNP_inclusion_filename, "r")
    SNPs_to_include = []
    for line in inclusion_file:
        cols = line.replace(',', ' ').rstrip().split()
        SNP_name = '{0}{1}{2}'.format(cols[0], CHROM_POS_SEPARATOR, cols[1])
        if SNP_name not in SNPs_to_include:
            SNPs_to_include.append(SNP_name)
    return SNPs_to_include


def get_genotypes(record, groups_populations, populations_individuals):
    """
    Get genotypes of all individuals for a particular SNP locus (from VCF
    record) organised by population.
    """
    gts_in_pop = {}
    # Iterate through the different groups / pops / indvs and store genotypes
    for group in sorted(groups_populations.keys()):
        # Iterate through the different populations of a group
        for population in groups_populations[group]:
            # Iterate through the different individuals of a population
            gts_in_pop[population] = []
            for individual in populations_individuals[population]:
                # Store SNP genotype in list
                genotype = record.genotype(individual)["GT"]
                gts_in_pop[population].append(genotype if genotype else
                                              MISSING_GENOTYPE)
    return gts_in_pop


def get_gt_counts(genotypes):
    """
    Get genotype counts from string of concatenated genotypes in the following
    order: [major allele homozygote, heterozygote, minor allele heterozygote].
    Only worls with bi-allelic data.
    """
    delimiter = genotypes[0][1]
    gt_counts = []
    concat_genos = ':'.join(genotypes)
    # Major allele homozygote
    gt_counts.append(concat_genos.count('0{0}0'.format(delimiter)))
    # Heterozygote count
    het_count1 = concat_genos.count('1{0}0'.format(delimiter))
    het_count2 = concat_genos.count('0{0}1'.format(delimiter))
    gt_counts.append(het_count1 + het_count2)
    # Minor allele homozygote
    gt_counts.append(concat_genos.count('1{0}1'.format(delimiter)))
    return gt_counts


def main(vcf_filename, factor_filename, SNP_inclusion_filename):
    # Create dictionaries for groups-populations and populations-individuals
    groups_populations, populations_individuals = \
        create_dicts_from_factor_file(factor_filename)
    # Create population list sorted by group
    pop_order = []
    for group in groups_populations:
        for pop in groups_populations[group]:
            pop_order.append(pop)

    # Create list with SNPs to be included
    SNPs_to_include = create_SNP_inclusion_list(SNP_inclusion_filename)

    # Output header
    print('chrom_pos\tgroup\tpop\tgt\tfreq')

    # Read through SNPs in VCF and calculate/output genotype freqs
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for record in vcf_reader:
        SNP_name = '{0}{1}{2}'.format(record.CHROM, CHROM_POS_SEPARATOR,
                                      record.POS)
        # Only consider SNP if listed in SNP_inclusion_file
        if SNP_name in SNPs_to_include:
            # Get all SNP genotypes for each pop
            gts_in_pop = get_genotypes(record, groups_populations,
                                       populations_individuals)
            gt_count = {}
            for group in groups_populations:
                for pop in groups_populations[group]:
                    gt_count[pop] = get_gt_counts(gts_in_pop[pop])
                    counts_ref = gt_count[pop][0]
                    counts_het = gt_count[pop][1]
                    counts_alt = gt_count[pop][2]

                    # Output hom_ref genotype freq
                    print('{0}\t{1}\t{2}\t{3}\t{4}'.format(SNP_name, group,
                                                           pop, NAME_HOM_REF,
                                                           counts_ref))
                    # Output het genotype freq
                    print('{0}\t{1}\t{2}\t{3}\t{4}'.format(SNP_name, group,
                                                           pop, NAME_HET,
                                                           counts_het))
                    # Output hom_alt genotype freq
                    print('{0}\t{1}\t{2}\t{3}\t{4}'.format(SNP_name, group,
                                                           pop, NAME_ALT_REF,
                                                           counts_alt))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('factor_filename', metavar='factor_file',
                        help='text file (tsv or csv) with individuals, their \
                        population assignment and group assignment')
    parser.add_argument('SNP_filename', metavar='SNP_file',
                        help='text file (tsv or csv) with CHROM/POS of each \
                        SNP to be outputted')
    args = parser.parse_args()
    main(args.vcf_filename, args.factor_filename, args.SNP_filename)
