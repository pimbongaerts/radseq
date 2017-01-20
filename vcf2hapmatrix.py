#!/usr/bin/env python
"""
Converts `.vcf` file to Tag Haplotype Matrix (with Chrom), with order of
individuals as indicated in optional file.

Note: not yet properly tested. SNPs of same CHROM (first column) in `.vcf`
should be grouped together/sequentially, and all individuals need to be listed
in order_file.
"""
import sys
import vcf
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


VCF_HEADER = '#CHROM'


def list_from_indivorder_file(indivorder_filename):
    """ Initialize list of indvs from textfile or popfile """
    indv_order = []
    indivorder_file = open(indivorder_filename, 'r')
    for line in indivorder_file:
        indv_order.append(line.rstrip().split()[0])  # ignore second col
    return indv_order


def list_from_vcf_file(vcf_filename):
    """ Initialize list of indvs from VCF file """
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        # Find VCF header and extract names and order of individuals
        if line[0:len(VCF_HEADER)] == VCF_HEADER:
            cols = line.split()
            indv_order = cols[9:]
            break
    vcf_file.close()
    return indv_order


def get_haplotypes_output(indv_order, haplotypes):
    """ Generate output line with haplotypes (order as in indv_order) """
    haplotypes_list = []
    for indiv in indv_order:
        haplotype_a = ''.join(haplotypes[indiv][0])
        haplotype_b = ''.join(haplotypes[indiv][1])
        # Only include genotype when no missing data present
        if '.' not in haplotype_a:
            haplotypes_list.append('{0}/{1}'.format(haplotype_a, haplotype_b))
        else:
            haplotypes_list.append('')
    # Join all genotypes into single tab-separated string
    haplotypes_output = '{0}\n'.format('\t'.join(haplotypes_list))
    return haplotypes_output


def main(vcf_filename, indivorder_filename):
    # Initialize list with order of individuals
    if not indivorder_filename:
        # Use order in VCF file if separate file not given
        indv_order = list_from_vcf_file(vcf_filename)
    else:
        # Use order as indicated in separate file
        indv_order = list_from_indivorder_file(indivorder_filename)

    # Output header
    print('Chr\t{0}'.format('\t'.join(indv_order)))

    # Iterate over loci/snps for each indv in VCF file
    previous_chrom = ""
    haplotypes = {}
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for record in vcf_reader:
        # Output haplotype when reaching new CHROM
        if previous_chrom != '' and record.CHROM != previous_chrom:
            haplotypes_output = get_haplotypes_output(indv_order, haplotypes)
            print('{0}\t{1}'.format(previous_chrom, haplotypes_output), end='')
            haplotypes.clear()

        # Store all SNP genotypes in dict (appending when same CHROM)
        for indiv in indv_order:
            if record.genotype(indiv).gt_bases:
                snp_genotypes = record.genotype(indiv).gt_bases \
                                      .replace('|', '/').split("/")
            else:
                snp_genotypes = ['.', '.']
            if indiv not in haplotypes:
                # Create dict item (indiv = key) with a list for each hap
                haplotypes[indiv] = {}
                haplotypes[indiv][0] = []
                haplotypes[indiv][1] = []
            haplotypes[indiv][0].append(snp_genotypes[0])
            haplotypes[indiv][1].append(snp_genotypes[1])

        previous_chrom = record.CHROM

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-o', '--order_file', metavar='order_file',
                        help='text file with preferred output order of \
                        individuals')
    args = parser.parse_args()
    main(args.vcf_filename, args.order_file)
