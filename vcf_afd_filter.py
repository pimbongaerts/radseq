#!/usr/bin/env python
"""
Calculate allele frequency differentials between groups, and flag
those loci that have AFDs exceeding threshold between all subgroups
of those groups. Note: so far only used with 3 groups - yet to be tested 
for more.
"""
import sys
import vcf
import itertools
import numpy as np
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

AFD_OUTLIER = 'AFD_OUTLIER'
NON_AFD_OUTLIER = ''

class SnpEvaluator(object):
    """ SNP evaluator object - initialized with group structure,
        and then reused to evaluate individual SNP data
    """

    def __init__(self):
        self.groups = {}
        self.chrom = ''
        self.pos = 0
        self.classification = ''
        self.avg_afds = {}

    def read_sample_groups(self, group_filename):
        """ Initialize sample groups from group_file (or pop_file) """
        group_file = open(group_filename, 'r')
        for line in group_file:
            cols = line.replace(',', ' ').split()
            sample, group, subgroup = cols[0], cols[1], cols[2]
            if group not in self.groups:
                self.groups[group] = Group(group)
            if subgroup not in self.groups[group].subgroups:
                self.groups[group].add_subgroup(subgroup)
            self.groups[group].subgroups[subgroup].add_sample(sample)
        group_file.close()

    def classify(self, vcf_record, afd_threshold):
        """ Classify SNP based on allele frequency differential """
        self.chrom = vcf_record.CHROM
        self.pos = vcf_record.POS
        
        # Obtain counts and  allele frequencies for each subgroup
        for group in sorted(self.groups):
            for subgroup in sorted(self.groups[group].subgroups):
                self.groups[group].subgroups[subgroup].read_gts(vcf_record)
        
        # Conduct all pairwise comparisons between subgroups of groups
        # and flag when all of them have an alllele frequency differential
        # (AFD) that >= afd_threshold
        group_list = []
        for group1, group2 in itertools.combinations(self.groups, 2):
            # Calculate allele frequency differential between all subgroups
            # of those two groups
            group1_subgroups = self.groups[group1].subgroups
            group2_subgroups = self.groups[group2].subgroups
            cross_group_afds = []
            for subgroup1, subgroup2 in itertools.product(group1_subgroups, 
                                                          group2_subgroups):
                afd = abs(self.groups[group1].subgroups[subgroup1].aaf -
                          self.groups[group2].subgroups[subgroup2].aaf)
                cross_group_afds.append(afd)
            # Add average of cross_group_afds to avg_afds
            groups = self.__get_comp_name(group1, group2)
            self.avg_afds[groups] = np.mean(cross_group_afds)
            # Count if all cross-comparisons between two groups are >=threshold
            if all(afd >= afd_threshold for afd in cross_group_afds):
                group_list.extend([group1, group2])
        
        # Flag as outlier when a group is different from all other groups
        # ( = when the group occurs in group_list == group_count - 1)
        outlier_groups = []
        for group in self.groups:
            if group_list.count(group) == len(self.groups) - 1:
                outlier_groups.append(group)

        # Store classification
        if len(outlier_groups):
            outlier_groups_concat = '/'.join(outlier_groups)
            self.classification = '{0}:{1}'.format(AFD_OUTLIER,
                                                   outlier_groups_concat)
        else:
            self.classification = NON_AFD_OUTLIER

    def print_output_header(self):
        """ Output header information """
        output_part1 = ['CHROM', 'POS', 'CLASSIFIER']
        output_part2 = []
        output_part3 = []
        output_part4 = []
        # Header cols for allele freq and count data
        for group in sorted(self.groups):
            for subgroup in sorted(self.groups[group].subgroups):
                output_part2.append('{0}[{1}]_AAF'.format(group, subgroup))
                output_part3.append('{0}[{1}]_N'.format(group, subgroup))
        # Header cols for allele freq differentials
        for group1, group2 in itertools.combinations(self.groups, 2):
            output_part4.append('{0}_AFD'.format(self.__get_comp_name(group1, 
                                                                      group2)))
        # Output header
        print('{0},{1},{2},{3}'.format(','.join(output_part1), 
                                       ','.join(output_part2),
                                       ','.join(output_part3),
                                       ','.join(output_part4)))

    def print_output_summary(self):
        """ Output information on current SNP """
        output_part1 = [self.chrom, str(self.pos), self.classification]
        output_part2 = []
        output_part3 = []
        # Allele freq and count data
        for group in sorted(self.groups):
            for subgroup in sorted(self.groups[group].subgroups):
                aaf = self.groups[group].subgroups[subgroup].aaf
                aaf_formatted = '{0:.4f}'.format(aaf)
                output_part1.append(aaf_formatted)
                output_part2.append(self.groups[group].subgroups[subgroup].aac)
        # Allele freq differentials
        for group1, group2 in itertools.combinations(self.groups, 2):
            afd = self.avg_afds[self.__get_comp_name(group1, group2)]
            output_part3.append('{0:.4f}'.format(afd))

        print('{0},{1},{2}'.format(','.join(output_part1), 
                                   ','.join(output_part2),
                                   ','.join(output_part3)))

    @staticmethod
    def __get_comp_name(group1, group2):
        return '_vs_'.join(sorted([group1, group2]))

class Group(object):
    """ Group """

    def __init__(self, name):
        self.name = name
        self.subgroups = {}
        self.classification = ''
        self.ref_present = False
        self.alt_present = False
        self.always_present_alleles = set()

    def add_subgroup(self, name):
        self.subgroups[name] = Subgroup(name)

class Subgroup(object):
    """ Subgroups """

    def __init__(self, name):
        self.name = name
        self.samples = []
        self.aaf = -1
        self.aac = ''

    def add_sample(self, name):
        self.samples.append(name)

    def read_gts(self, vcf_record):
        """ Get frequency and count of alternative allele for a subgroup """
        alleles = []
        # Extract all alleles from genotypes
        for sample in self.samples:
            genotype = vcf_record.genotype(sample)['GT']
            if genotype:
                alleles.append(genotype[0])
                alleles.append(genotype[2])
        # Extract data from list of alleles
        cat_alleles = ''.join(alleles)
        unique_alleles = ''.join(set(cat_alleles))
        if len(cat_alleles) == 0:
            # Error - only missing data
            sys.exit('Only missing data for subgroup {0}'
                     ' and SNP {1}:{2}'.format(self.name, 
                                               vcf_record.CHROM,
                                               vcf_record.POS))
        elif len(unique_alleles) > 2:
            # Error - SNP has more than 2 alleles'
            sys.exit('More than 2 alleles for subgroup {0}'
                     ' and SNP {1}:{2}'.format(self.name, 
                                               vcf_record.CHROM,
                                               vcf_record.POS))
        else:
            self.aaf = float(cat_alleles.count('1') / len(cat_alleles))
            self.aac = '{0} of {1}'.format(cat_alleles.count('1'), 
                                           len(cat_alleles))


def main(vcf_filename, group_filename, afd_threshold):
    # Initialize evaluator and read in group_file
    snp_evaluator = SnpEvaluator()
    snp_evaluator.read_sample_groups(group_filename)

    # Iterate over vcf and evaluate snps
    snp_evaluator.print_output_header()
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for vcf_record in vcf_reader:
        snp_evaluator.classify(vcf_record, afd_threshold)
        snp_evaluator.print_output_summary()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('group_filename', metavar='group_file',
                        help='text file (tsv or csv) separating individuals \
                              (first column) into groups (second column)')
    parser.add_argument('afd_threshold', metavar='afd_threshold',
                        help='allele frequency differential threshold')
    args = parser.parse_args()
    main(args.vcf_filename, args.group_filename, float(args.afd_threshold))
