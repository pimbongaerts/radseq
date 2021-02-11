#!/usr/bin/env python
"""
Filters `.vcf` file for SNPs that are genotyped for a minimum number of
individuals in each of the populations (rather than overall proportion of
individuals). This can help to guarantee a minimum number of individuals
to calculate population-based statistics, and eliminate loci that might
be suffering from locus drop-out in particular populations.

Note: only individuals that are listed in popfile are taken into account to
determine number of individuals genotyped (but all indivs are outputted).
"""
import sys
import vcf
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


MISSING_DATA_CHAR = '.'
HEADER_CHAR = '#'


def dict_from_popfile(pop_filename):
    """ Initialise dict of pops with lists of indvs from popfile """
    pops_indivs = {}
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        cols = line.rstrip().split()
        indiv = cols[0]
        pop = cols[1]
        if pop not in pops_indivs:
            pops_indivs[pop] = []
        pops_indivs[pop].append(indiv)
    return pops_indivs


def num_in_pop_genotyped(record, indvs):
    """ Count number of indivs genotyped in pop for current SNP """
    count_genotyped = count_missing = 0
    for indiv in indvs:
        vcf_genotype = record.genotype(indiv)['GT']
        if not vcf_genotype or MISSING_DATA_CHAR in vcf_genotype:
            count_missing += 1
        elif any(char.isdigit() for char in vcf_genotype):
            count_genotyped += 1
        else:
            sys.exit('Unexpected value: {0}:{1}'.format(indiv, record.CHROM))
    return int(count_genotyped)


def prop_of_overall_genotyped(record, pops_indivs):
    """ Calculate proportion of indivs genotyped overall for current SNP """
    count_genotyped = count_missing = 0
    for pop in pops_indivs:
        for indiv in pops_indivs[pop]:
            vcf_genotype = record.genotype(indiv)['GT']
            if not vcf_genotype or MISSING_DATA_CHAR in vcf_genotype:
                count_missing += 1
            elif any(char.isdigit() for char in vcf_genotype):
                count_genotyped += 1
            else:
                sys.exit('Unexpected value: {0}:{1}'.format(indiv,
                                                            record.CHROM))
    return float(count_genotyped / (count_genotyped + count_missing))


def main(vcf_filename, pop_filename, threshold, output_filename):
    # Read list of indvs and population assignments from popfile
    print('Reading popfile {0}...'.format(pop_filename))
    pops_indivs = dict_from_popfile(pop_filename)

    # Iterate through loci and output only those that meet the minimum number
    # of samples genotyped for each population in pop_filename
    print('Opening vcf {0}...'.format(vcf_filename))
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    snps_failing_threshold = []
    min_overall_genotyped = 1.0
    overall_snp_count = 0
    previous_CHROM = ''
    print('Evaluating SNPs in vcf...')
    for record in vcf_reader:
        failing_snp = False
        overall_snp_count += 1
        if record.CHROM != previous_CHROM:
            previous_CHROM = record.CHROM
        # Assess all pops until pop is reached that fails threshold
        for pop in pops_indivs:
            if num_in_pop_genotyped(record, pops_indivs[pop]) < threshold:
                failing_snp = True
                break
        # If failing treshold add to list of SNPs to be removed
        if failing_snp:
            snps_failing_threshold.append('{0}\t{1}'.format(record.CHROM,
                                                            record.POS))
        # Otherwise assess overall percentage genotyped to establish minimum
        else:
            overall_genotyped = prop_of_overall_genotyped(record, pops_indivs)
            min_overall_genotyped = min(min_overall_genotyped,
                                        overall_genotyped)

    # Remove failing SNPs from VCF
    vcf_file = open(vcf_filename, 'r')
    new_vcf_file = open(output_filename, 'w')
    for line in vcf_file:
        # Output header
        if line[0] == HEADER_CHAR:
            new_vcf_file.write(line)
        else:
            # Output SNP only if not in list of SNPs to be removed
            cols = line.rstrip().split()
            current_snp = '{0}\t{1}'.format(cols[0], cols[1])
            if current_snp not in snps_failing_threshold:
                new_vcf_file.write(line)
    vcf_file.close()
    new_vcf_file.close()

    # Output summary
    print(('{0} out of {1} SNPs '
           'failing {2} threshold').format(len(snps_failing_threshold),
                                           overall_snp_count, threshold))
    print(('Lowest overall proportion genotyped '
           'of remaining SNPs: {0}').format(min_overall_genotyped))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    parser.add_argument('threshold', metavar='min_proportion', type=float,
                        help='proportion of individuals required to be \
                        genotyped in each population for a SNP to be \
                        included (e.g `0.8` for 80 percent of individuals)')
    parser.add_argument('output_filename', metavar='output_filename',
                        help='name of output file (`.vcf`)')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.threshold,
         args.output_filename)
