#!/usr/bin/env python
"""
Converts `.vcf` file to OutFLANK input files (3 different files).

Note: not yet properly tested.
"""
import sys
import argparse
import vcf

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def dict_from_popfile(pop_filename):
    # Read list of indvs and population assignments from popfile
    indvs_pops = {}
    pop_file = open(pop_filename, 'r')
    for line in pop_file:
        line = line.rstrip()
        cols = line.split('\t')
        indvs_pops[cols[0]] = cols[1]

    return indvs_pops


def get_outflank_genotype(vcf_genotype):
    # Convert VCF to outflank genotype
    if not vcf_genotype:
        return '9'		        # Missing data (=9)
    elif vcf_genotype[0] == '0' and vcf_genotype[2] == '0':
        return '0'              # Homozygote for 0 (=0)
    elif vcf_genotype[0] == '0' and vcf_genotype[2] == '1':
        return '1'              # Heterozygote 0/1 (=1)
    elif vcf_genotype[0] == '1' and vcf_genotype[2] == '0':
        return '1'	            # Heterozygote 0/1 (=1)
    elif vcf_genotype[0] == '1' and vcf_genotype[2] == '1':
        return '2'              # Homozygote for 1 (=2)
    else:
        sys.exit('Unexpected genotype in vcf: {0}\n'.format(vcf_genotype))


def main(vcf_filename, pop_filename, output_prefix):
    # Read popfile with list of indvs and population assignments from
    # popfile (only those indvs in popfile are included)
    indvs_pops = dict_from_popfile(pop_filename)

    # Initialise dictionary with indvs (keys) and list of genotypes (values)
    outflank_genotypes = {}
    for individual in indvs_pops:
        outflank_genotypes[individual] = []

    # Iterate over loci/indvs in VCF file
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    outflank_loci_file = open('{0}_outflank_loci.txt'.format(output_prefix),
                              'w')
    for record in vcf_reader:
        outflank_loci_file.write('{0}:{1}\n'.format(record.CHROM, record.POS))
        for individual in indvs_pops:
            vcf_genotype = record.genotype(individual)['GT']
            outflank_genotype = get_outflank_genotype(vcf_genotype)
            outflank_genotypes[individual].append(outflank_genotype)
    outflank_loci_file.close()

    # Output genotypes for each individual
    outflank_file = open('{0}_outflank_data.txt'.format(output_prefix), 'w')
    outflank_pops_file = open('{0}_outflank_pops.txt'.format(output_prefix),
                              'w')
    for individual in sorted(indvs_pops.keys()):
        outflank_pops_file.write('{0}\n'.format(indvs_pops[individual]))
        concat_outflank_genotypes = ' '.join(outflank_genotypes[individual])
        outflank_file.write('{0}\n'.format(concat_outflank_genotypes))
    outflank_file.close()
    outflank_pops_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        help='text file (tsv or csv) with individuals and \
                              populations')
    parser.add_argument('output_prefix', metavar='output_prefix',
                        help='prefix for output files')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.output_prefix)
