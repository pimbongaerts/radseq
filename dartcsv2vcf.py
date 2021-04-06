#!/usr/bin/env python
"""
dartcsv2vcf.py: Converts DArT-seq CSV file to basic vcf
"""
import sys
import argparse

VCF_HEADER = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

CSV_COMMENT = '*'
CSV_HEADER = 'AlleleID'

CSV_COL_ALLELEID= 0
CSV_COL_CLONEID = 1
CSV_COL_ALLELE_SEQUENCE = 2
CSV_COL_TRIMMED_SEQUENCE = 3
CSV_COL_SNP = 4
CSV_COL_SNP_POSITION = 5
CSV_COL_FIRST_SAMPLE = 18

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2021 Pim Bongaerts'
__license__ = 'GPL'


def get_output_filename(csv_filename, extension):
    # Return output filename based on original csv filename
    return '{0}.{1}'.format(csv_filename.replace('.csv', ''), extension)

def get_dict_of_name_changes(samplenames_filename):
    """ Initialise dictionary with old (keys) and new (values)  """
    rename_file = open(samplenames_filename, 'r')
    name_changes = {}
    for line in rename_file:
        cols = line.replace(',', '\t').split()
        name_changes[cols[0].strip()] = cols[1].strip()
    rename_file.close()
    return name_changes

def get_new_sample_names(old_samples, name_changes):
    """ Change samples to new sample names """
    sample_names = []
    for sample in old_samples:
        indiv_key = sample.strip()
        if sample.strip() in name_changes:
            sample_names.append(name_changes[sample.strip()])
        else:
            print('Not renamed: {0}'.format(sample))
            sample_names.append(sample.strip())
    return sample_names

def convert_to_vcf_genotype(genotype):
    """ Convert DaRT to VCF genotype """
    dart_gt = genotype.strip()
    if dart_gt == '-':
        vcf_gt = './.'
    elif dart_gt == '0':
        vcf_gt = '0/0'
    elif dart_gt == '1':
        vcf_gt = '1/1'
    elif dart_gt == '2':
        vcf_gt = '0/1'
    else:
        print(dart_gt)
        vcf_gt = 'err'
    return vcf_gt

def get_vcf_genotypes(dart_genotypes):
    """ Get VCF genotypes"""
    return '\t'.join([convert_to_vcf_genotype(dart_gt) for dart_gt in dart_genotypes])

def output_vcf_line(vcf_file, chrom, pos, ref, alt, vcf_genotypes_concat):
    """ Output VCF line """
    vcf_file.write('{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.\tGT\t{4}\n'.format(chrom, pos, ref, alt, vcf_genotypes_concat))

def main(dartcsv_filename, samplenames_filename):
    # Generate look-up dict with sample names
    name_changes = get_dict_of_name_changes(samplenames_filename)
    
    # Open VCF output file
    vcf_file = open(get_output_filename(dartcsv_filename, 'vcf'), 'w')
    fasta_file = open(get_output_filename(dartcsv_filename, 'fa'), 'w')

    # Iterate over CSV input file and output to VCF
    csv_file = open(dartcsv_filename, 'r')
    for line in csv_file:
        cols = line.split(',')
        if line[0:len(CSV_HEADER)] == CSV_HEADER:
            sample_names = get_new_sample_names(cols[CSV_COL_FIRST_SAMPLE:], name_changes)
            # Output VCF header
            vcf_file.write('{0}\t{1}\n'.format(VCF_HEADER, '\t'.join(sample_names)))
        elif line[0] != CSV_COMMENT:
            # Obtain values for VCF output
            ref_allele, divider, alt_allele = cols[CSV_COL_SNP].split(':')[1]
            # Convert genotypes to VCF format
            vcf_genotypes_concat = get_vcf_genotypes(cols[CSV_COL_FIRST_SAMPLE:])
            # Output SNP line to VCF
            output_vcf_line(vcf_file, cols[CSV_COL_CLONEID], cols[CSV_COL_SNP_POSITION], 
                            ref_allele, alt_allele, vcf_genotypes_concat)
            # Output allele sequence to FASTA
            fasta_file.write('>{0}\n'.format(cols[CSV_COL_CLONEID]))
            fasta_file.write('>{0}\n'.format(cols[CSV_COL_ALLELE_SEQUENCE]))
   
    csv_file.close()      
    vcf_file.close()
    fasta_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dartcsv_filename', metavar='dartcsv_filename',
                        help='DArT-seq csv filename`)')
    parser.add_argument('samplenames_filename', metavar='samplenames_file',
                        help='text file (tsv or csv) with old and new name \
                        for each sample (not all samples have to be listed)')
    args = parser.parse_args()
    main(args.dartcsv_filename, args.samplenames_filename)