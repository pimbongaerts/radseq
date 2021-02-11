#!/usr/bin/env python
"""
Creates a genotype matrix
"""
import sys
import vcf
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

COL_INDIV = 0
COL_POP = 1
COL_ASSIGN = 2
MISSING_DATA = -9
ERROR_FLAG = "ERROR"


def read_pop_file(pop_filename):
    """ Read individuals and assignments from pop_file """
    indiv_order = []
    indiv_pops = {}
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        cols = line.replace(',', ' ').split()
        indiv_order.append(cols[0])
        indiv_pops[cols[0]] = cols[1]
    return indiv_order, indiv_pops

def read_locus_file(locus_filename):
    """ Read loci and assignments from locus_file"""
    locus_order = []
    locus_groups = {}
    locus_file = open(locus_filename, 'r')
    for line in locus_file:
        cols = line.replace(',', ' ').split()
        locus_name = '\t'.join([cols[0], str(cols[1])])
        locus_order.append(locus_name)
        locus_groups[locus_name] = cols[2]
    return locus_order, locus_groups

def get_all_gts(record, list_of_indivs, switch_ref_alt):
    """ Outputs gts of all indivs for current SNP to output_file """
    gts_for_snp = []
    for indiv in list_of_indivs:
        gt = record.genotype(indiv).gt_type
        if gt is None:
            # Missing data char if gt is missing
            gt = MISSING_DATA
        elif switch_ref_alt:
            # Switch REF/ALT if necessary
            if gt == 0:
                gt = 2
            elif gt == 2:
                gt = 0
        gts_for_snp.append(str(gt))
    return gts_for_snp

def ref_alt_switch_necessary(record, par2_inds):
    """ Assess whether REF/ALT need to be switched (so that par1 ~ ref) """
    gts = []
    for indiv in par2_inds:
        gts.append(str(record.genotype(indiv).gt_type))
    concat_gts = ''.join(gts)
    return concat_gts.count('0') > concat_gts.count('2')


def obtain_snp_genotypes(vcf_filename, indiv_order, locus_order):
    """ Obtain genotypes for SNPs that meet criteria & listed indivs """
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    snp_gts = {}
    for record in vcf_reader:
        marker = '\t'.join([record.CHROM, str(record.POS)])
        if marker in locus_order:
            snp_gts[marker] = get_all_gts(record, indiv_order, False)
    return snp_gts


def main(vcf_filename, pop_filename, locus_filename):
    # Obtain order and grouping of individuals
    indiv_order, indiv_pops = read_pop_file(pop_filename)

    # Obtain order and grouping of loci
    locus_order, locus_groups = read_locus_file(locus_filename)

    # Obtain genotypes for SNPs that meet criteria & listed indivs
    snp_gts = obtain_snp_genotypes(vcf_filename, indiv_order, locus_order)

    # Output header
    print('chrom\tpos\ttype\t{0}'.format('\t'.join(indiv_order)))

    # Output SNPs
    for snp in locus_order:
        snp_formatted = '{0}\t{1}'.format(snp.split()[0], snp.split()[1])
        print('{0}\t{1}\t{2}'.format(snp_formatted,
                                     locus_groups[snp],
                                        '\t'.join(snp_gts[snp])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with group assignments \
                              of individuals, and in the order that is \
                              required for the genotype matrix')
    parser.add_argument('locus_filename', metavar='locus_file',
                        help='text file (tsv or csv) indicating which loci \
                        should be included (first column), what group \
                        they belong to (second column), and their order.')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.locus_filename)
