#!/usr/bin/env python
"""
Converts `.vcf` file to INTROGRESS input files (4 files). Splits data into
three categories: parental1, parental2 and admixed based on cluster assignment
(provided in separate file; e.g. STRUCTURE output) and given threshold, and
outputs data for loc that exceed a certain frequency difference between the
two 'parental' categories.

Note: not yet properly tested. also see similar `vcf_ancestry_matrix.py`
script. I use the formatted CLUMPP output (`clumpp_K2.out.csv`) from the
`structure_mp` wrapper as assignment file (max. of 2 clusters).
"""
import sys
import argparse
import vcf

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

ADMIX_LABEL = "ADMIX"
COL_INDIV = 0
COL_POP = 1
COL_ASSIGN = 2
ERROR_FLAG = "ERROR"


def read_assignment_file(
        assignment_filename, assign_cut_off, include_parentals):
    """ Assess assignment of each individual and store in dict/list """
    assignment_file = open(assignment_filename, 'r')
    par1_indivs = []
    par2_indivs = []
    admixed_indivs = []
    indiv_pops = {}
    for line in assignment_file:
        cols = line.split(',')
        indiv_pops[cols[COL_INDIV]] = cols[COL_POP]
        if float(cols[COL_ASSIGN]) > float(assign_cut_off):
            # Parental individuals 1
            par1_indivs.append(cols[COL_INDIV])
        elif float(cols[COL_ASSIGN]) < (1 - float(assign_cut_off)):
            # Parental individuals 2
            par2_indivs.append(cols[COL_INDIV])
        else:
            # Admixed individuals
            admixed_indivs.append(cols[COL_INDIV])

    # Include parentals to admixed individuals if include flag is set
    if include_parentals:
        final_admixed_indivs = par1_indivs + admixed_indivs + par2_indivs
    else:
        final_admixed_indivs = admixed_indivs

    return par1_indivs, par2_indivs, final_admixed_indivs, indiv_pops


def output_individuals(record, list_of_indivs, output_file):
    """ Outputs genotypes of all indivs for current SNP to output_file """
    genotypes_for_snp = []
    for indiv in list_of_indivs:
        genotype = record.genotype(indiv).gt_bases
        if genotype:
            introgress_genotype = '{0}/{1}'.format(genotype[0], genotype[2])
        else:
            introgress_genotype = 'NA/NA'
        genotypes_for_snp.append(introgress_genotype)
    concat_genotypes_for_snp = ','.join(genotypes_for_snp)
    output_file.write('{}\n'.format(concat_genotypes_for_snp))


def calculate_allele_frequency(record, indivs):
    """ Calculate allele frequency """
    # Create list of all alleles in population
    alleles = []
    for indiv in indivs:
        genotype = record.genotype(indiv)['GT']
        if genotype:
            alleles.append(genotype[0])
            alleles.append(genotype[2])
    concat_alleles = ''.join(alleles)

    # Assess whether there are more than two alleles
    unique_alleles = ''.join(set(concat_alleles))
    if len(concat_alleles) == 0:
        error_message = 'Error - only missing data for'
        print('{0}: {1}_{2}'.format(error_message, record.CHROM, record.POS))
        allele_freq = ERROR_FLAG
    elif len(unique_alleles) > 2:
        error_message = 'Error - locus has more than 2 alleles'
        print('{0}: {1}_{2}'.format(error_message, record.CHROM, record.POS))
        allele_freq = ERROR_FLAG
    else:
        allele_freq = float(concat_alleles.count('1') / len(concat_alleles))
    return allele_freq


def allele_frequency_difference(record, par1_indivs,
                                par2_indivs, freq_diff_cut_off):
    """ Calculates absolute difference in allele frequencies """
    p1_allele_freq = calculate_allele_frequency(record, par1_indivs)
    p2_allele_freq = calculate_allele_frequency(record, par2_indivs)
    if p1_allele_freq == ERROR_FLAG or p2_allele_freq == ERROR_FLAG:
        return ERROR_FLAG
    else:
        return abs(p1_allele_freq - p2_allele_freq)


def main(vcf_filename, assignment_filename, assign_cut_off, freq_diff_cut_off,
         output_prefix, include_parentals):
    # Assess assignment of individuals in dataset and store in lists
    (par1_indivs, par2_indivs, admixed_indivs,
     indiv_pops) = read_assignment_file(assignment_filename, assign_cut_off,
                                        include_parentals)

    # Open output files
    parental1_file = open('{}_par1'.format(output_prefix), 'w')
    parental2_file = open('{}_par2'.format(output_prefix), 'w')
    admixed_file = open('{}_admix'.format(output_prefix), 'w')
    loci_file = open('{}_loci'.format(output_prefix), 'w')

    # Write output file headers
    loci_file.write('locus,type,lg,marker pos.\n')
    admixed_file.write(
        '{}\n'.format(
            ','.join(
                ['ADMIX'] *
                len(admixed_indivs))))
    admixed_file.write('{}\n'.format(','.join(admixed_indivs)))

    # Iterate over loci in VCF file and output to files
    output_count = 0
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    for record in vcf_reader:
        # Calculate absolute difference in allele frequencies
        freq_diff = allele_frequency_difference(record, par1_indivs,
                                                par2_indivs,
                                                freq_diff_cut_off)
        # Only include locus if above threshold frequency_diff_cut_off
        if freq_diff is not ERROR_FLAG and freq_diff > float(
                freq_diff_cut_off):
            # Output locus information (CHROM as lg and POS as marker pos.)
            marker = '{0}_{1}'.format(record.CHROM, record.POS)
            loci_file.write('{0},C,{1},{1}.{2}\n'.format(marker, record.CHROM,
                                                         record.POS))
            # Output genotypes of individuals to different files
            output_individuals(record, par1_indivs, parental1_file)
            output_individuals(record, par2_indivs, parental2_file)
            output_individuals(record, admixed_indivs, admixed_file)
            output_count += 1

    # Close output files
    parental1_file.close()
    parental2_file.close()
    admixed_file.close()
    loci_file.close()

    # Summary message
    print('Outputted {0} loci'.format(output_count))

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
                        putative purebred)')
    parser.add_argument('freq_cut_off', metavar='freq_cut_off', type=float,
                        help='min. allele frequency difference between the 2 \
                        clusters for a locus to be included in the output')
    parser.add_argument('output_prefix', metavar='output_prefix',
                        help='prefix for output files')
    parser.add_argument('--include', '-i', action='store_true',
                        help='set this flag if parental pops need to be \
                        included in output')
    args = parser.parse_args()
    main(args.vcf_filename, args.assignment_filename, args.assign_cut_off,
         args.freq_cut_off, args.output_prefix, args.include)
