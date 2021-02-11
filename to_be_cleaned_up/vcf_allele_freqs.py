#!/usr/bin/env python
"""

"""
import sys
import vcf
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

REF_FIXED = "REF_FIXED"
REF_NEARFIXED = "REF_NEARFIXED"
ALT_FIXED = "ALT_FIXED"
ALT_NEARFIXED = "ALT_NEARFIXED"
ALT_ALWAYS_PRESENT = "ALT_ALWAYS_PRESENT"
OTHER = "OTHER"

class SnpEvaluator(object):
    """ SNP evaluator object - initialized with group structure,
        and then reused to evaluate individual SNP data
    """

    def __init__(self):
        self.groups = {}
        self.chrom = ''
        self.pos = 0
        self.classification = ''

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

    def classify(self, vcf_record, lower_threshold, upper_threshold):
        """ Classify SNP based on alternative allele frequencies """
        self.chrom = vcf_record.CHROM
        self.pos = vcf_record.POS
        classifications = []
        for group in sorted(self.groups):
            # Get allele frequencies and counts for each subgroup
            for subgroup in sorted(self.groups[group].subgroups):
                self.groups[group].subgroups[subgroup].read_gts(vcf_record)
            # Classify group
            self.groups[group].classify(lower_threshold, upper_threshold)
            classifications.append(self.groups[group].classification)
        
        # Classify SNP
        group_count = len(classifications)
        ref_fixed_count = (classifications.count(REF_FIXED) +
                           classifications.count(REF_NEARFIXED))
        alt_fixed_count = (classifications.count(ALT_FIXED) +
                           classifications.count(ALT_NEARFIXED))
        ref_fixed_count_stringent = classifications.count(REF_FIXED)
        alt_fixed_count_stringent = classifications.count(ALT_FIXED)

        # If one group is REF_FIXED or REF_NEARFIXED and all the other groups 
        # are ALT_FIXED or ALT_NEAR_FIXED -> REF FIXED PRIVATE ALLELE
        if (ref_fixed_count == 1 and alt_fixed_count == (group_count - 1)):
            for group in self.groups:
                if self.groups[group].classification in (REF_FIXED, 
                                                         REF_NEARFIXED):
                    break
            self.classification = '{0}: {1}'.format(group, REF_FIXED)
        # If one group is ALT_FIXED or ALT_NEARFIXED and all the other groups 
        # are REF_FIXED or REF_NEAR_FIXED -> ALT FIXED PRIVATE ALLELE
        elif (alt_fixed_count == 1 and ref_fixed_count == (group_count - 1)):
            for group in self.groups:
                if self.groups[group].classification in (ALT_FIXED, 
                                                         ALT_NEARFIXED):
                    break
            self.classification = '{0}: {1}'.format(group, ALT_FIXED)
        else:
            self.classification = '.'

    def print_output_header(self):
        """ Output header information """
        output_part1 = ['CHROM', 'POS', 'CLASSIFIER']
        output_part2 = []
        output_part3 = []
        for group in sorted(self.groups):
            for subgroup in sorted(self.groups[group].subgroups):
                output_part2.append('{0}[{1}]_AAF'.format(group, subgroup))
                output_part3.append('{0}[{1}]_N'.format(group, subgroup))
        print('{0},{1},{2}'.format(','.join(output_part1), 
                                   ','.join(output_part2),
                                   ','.join(output_part3)))

    def print_output_summary(self):
        """ Output information on current SNP """
        output_part1 = [self.chrom, str(self.pos), self.classification]
        output_part2 = []
        for group in sorted(self.groups):
            for subgroup in sorted(self.groups[group].subgroups):
                aaf = self.groups[group].subgroups[subgroup].aaf
                aaf_formatted = '{0:.2f}'.format(aaf)
                output_part1.append(aaf_formatted)
                output_part2.append(self.groups[group].subgroups[subgroup].aac)
        print('{0},{1}'.format(','.join(output_part1), ','.join(output_part2)))

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

    def classify(self, lower_threshold, upper_threshold):
        """ Classify group based on allele frequencies of each subgroup """
        allele_freqs = []

        for subgroup in sorted(self.subgroups):
            allele_freqs.append(self.subgroups[subgroup].aaf)
        
        if all(freq == 0 for freq in allele_freqs):
            self.classification = REF_FIXED
        elif all(freq <= lower_threshold for freq in allele_freqs):
            self.classification = REF_NEARFIXED
        elif all(freq == 1 for freq in allele_freqs):
            self.classification = ALT_FIXED
        elif all(freq >= upper_threshold for freq in allele_freqs):
            self.classification = ALT_NEARFIXED
        elif all(freq > 0 for freq in allele_freqs):
            self.classification = ALT_ALWAYS_PRESENT
        else:
            self.classification = OTHER

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


def main(vcf_filename, group_filename, lower_threshold, upper_threshold):
    # Initialize evaluator and read in group_file
    snp_evaluator = SnpEvaluator()
    snp_evaluator.read_sample_groups(group_filename)

    # Iterate over vcf and evaluate snps
    snp_evaluator.print_output_header()
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for vcf_record in vcf_reader:
        snp_evaluator.classify(vcf_record, lower_threshold, upper_threshold)
        snp_evaluator.print_output_summary()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('group_filename', metavar='group_file',
                        help='text file (tsv or csv) separating individuals \
                              (first column) into groups (second column)')
    parser.add_argument('lower_threshold', metavar='lower_threshold',
                        help='lower_threshold')
    parser.add_argument('upper_threshold', metavar='upper_threshold',
                        help='upper_threshold')
    args = parser.parse_args()
    main(args.vcf_filename, args.group_filename, 
         float(args.lower_threshold), float(args.upper_threshold))
