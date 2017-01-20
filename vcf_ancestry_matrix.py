#!/usr/bin/env python
"""
Creates a genotype matrix for loci that have a large allele frequency
difference between two genetic clusters (as identified with e.g. STRUCTURE).
The script takes both a `.vcf` file and a text file with the assignment
probabilities as input. An assignment threshold (e.g. 0.98) needs to be
supplied to identify the reference individuals in the two clusters, and an
allele frequency cut-off needs to be supplied to identify divergent loci. An
optional file can be supplied with a list of loci that need to be included
regardless (e.g. previously identified outliers).

Note: I use the formatted CLUMPP output (`clumpp_K2.out.csv`) from the
`structure_mp` wrapper as assignment file (max. of 2 clusters).
"""
import sys
import vcf
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

COL_INDIV = 0
COL_POP = 1
COL_ASSIGN = 2
MISSING_DATA = -9
ERROR_FLAG = "ERROR"


def read_assignment_file(assignment_filename, assign_cut_off):
    """ Assess assignment of each individual and store in dict/list """
    par1_inds = []
    par2_inds = []
    all_inds = []
    assignment_file = open(assignment_filename, 'r')
    for line in assignment_file:
        cols = line.replace('\t', ',').split(',')
        if float(cols[COL_ASSIGN]) > float(assign_cut_off):
            # Parental individuals 1
            par1_inds.append(cols[COL_INDIV])
        elif float(cols[COL_ASSIGN]) < (1 - float(assign_cut_off)):
            # Parental individuals 2
            par2_inds.append(cols[COL_INDIV])
        # All individuals
        all_inds.append(cols[COL_INDIV])

    # Include parentals to admixed individuals if include_parentals = true
    return par1_inds, par2_inds, all_inds


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


def calculate_allele_frequency(record, indivs):
    """ Calculate allele frequency """
    alleles = []
    for indiv in indivs:
        gt = record.genotype(indiv)['GT']
        if gt:
            alleles.append(gt[0])
            alleles.append(gt[2])
    concat_alleles = ''.join(alleles)
    unique_alleles = ''.join(set(concat_alleles))
    if len(concat_alleles) == 0:
        # Error - only missing data
        allele_freq = ERROR_FLAG
    elif len(unique_alleles) > 2:
        # Error - SNP has more than 2 alleles'
        allele_freq = ERROR_FLAG
    else:
        allele_freq = float(concat_alleles.count('1') / len(concat_alleles))
    return allele_freq


def meets_freq_cut_off(record, par1_inds, par2_inds, freq_cut_off):
    """ Assess whether abs difference in allele freqs exceeds threshold """
    p1_allele_freq = calculate_allele_frequency(record, par1_inds)
    p2_allele_freq = calculate_allele_frequency(record, par2_inds)
    if p1_allele_freq == ERROR_FLAG or p2_allele_freq == ERROR_FLAG:
        return False
    else:
        return abs(p1_allele_freq - p2_allele_freq) >= float(freq_cut_off)


def ref_alt_switch_necessary(record, par2_inds):
    """ Assess whether REF/ALT need to be switched (so that par1 ~ ref) """
    gts = []
    for indiv in par2_inds:
        gts.append(str(record.genotype(indiv).gt_type))
    concat_gts = ''.join(gts)
    return concat_gts.count('0') > concat_gts.count('2')


def obtain_snp_genotypes(vcf_filename, outlier_snps, assignments,
                         freq_cut_off):
    """ Obtain genotypes for SNPs that meet criteria & listed indivs """
    par1_inds, par2_inds, all_inds = assignments
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    snp_gts = {}
    for record in vcf_reader:
        marker = '{0}\t{1}'.format(record.CHROM, record.POS)

        # Assess if SNP is either listed as outlier or exceeds threshold
        if marker in outlier_snps or \
           meets_freq_cut_off(record, par1_inds, par2_inds, freq_cut_off):

            # Assess whether REF/ALT alleles need to be switched
            switch_ref_alt = ref_alt_switch_necessary(record, par2_inds)

            # If so, store in dictionary
            snp_gts[marker] = get_all_gts(record, all_inds, switch_ref_alt)

    return snp_gts


def main(vcf_filename, assignment_filename, inclusion_file, assign_cut_off,
         freq_cut_off):
    # Assess assignment of individuals in dataset and store in lists
    assignments = read_assignment_file(assignment_filename, assign_cut_off)
    par1_inds, par2_inds, all_inds = assignments

    # Read list out outlier snps that should be included regardless of freq
    if inclusion_file:
        with open(inclusion_file) as file:
            outlier_snps = [line.strip() for line in file]
    else:
        outlier_snps = []

    # Obtain genotypes for SNPs that meet criteria & listed indivs
    snp_gts = obtain_snp_genotypes(vcf_filename, outlier_snps, assignments,
                                   freq_cut_off)

    # Output header
    print('chrom_pos\ttype\t{0}'.format('\t'.join(all_inds)))

    # Output outlier SNPs first
    for outlier_snp in outlier_snps:
        outlier_snp_formatted = '{0}_{1}'.format(outlier_snp.split()[0],
                                                 outlier_snp.split()[1])
        print('{0}\toutlier\t{1}'.format(outlier_snp_formatted,
                                         '\t'.join(snp_gts[outlier_snp])))

    # Output other SNPs
    for snp in sorted(snp_gts.keys()):
        if snp not in outlier_snps:
            snp_formatted = '{0}_{1}'.format(snp.split()[0], snp.split()[1])
            print(
                '{0}\tdeltap\t{1}'.format(snp_formatted,
                                          '\t'.join(snp_gts[snp])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('assignment_filename', metavar='assignment_file',
                        help='text file (tsv or csv) with assignment \
                        values for each individual (max. 2 clusters); e.g. a \
                        reformatted STRUCTURE output file')
    parser.add_argument('assign_cut_off', metavar='assign_cut_off', type=float,
                        help='min. assignment value for an individual to be \
                        included in the allele frequency calculation (i.e. \
                        putative purebred')
    parser.add_argument('freq_cut_off', metavar='freq_cut_off', type=float,
                        help='min. allele frequency difference between the 2 \
                        clusters for a locus to be included in the output')
    parser.add_argument('--include', '-i', metavar='inclusion_file',
                        help='text file with loci to be included \
                        in output regardless of allele frequency differences')
    args = parser.parse_args()
    main(args.vcf_filename, args.assignment_filename,
         args.include, args.assign_cut_off, args.freq_cut_off)
