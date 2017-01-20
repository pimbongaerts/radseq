#!/usr/bin/env python
"""
Contrast all samples in `.vcf` file against certain reference sample(s) (e.g.
outgroup samples), to assess for fixed / private alleles.
"""
import sys
import argparse
import vcf

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


# TODO: NEED TO IMPROVE OUTPUT DESCRIPTION

def get_genotypes(record, reference_samples):
    """ Get all alleles for a particular SNP """
    gts_ref = []
    gts_other = []

    for sample in record.samples:
        if sample['GT'] is not None:
            if sample.sample in reference_samples:
                gts_ref.append(sample['GT'][0])
                gts_ref.append(sample['GT'][2])
            else:
                gts_other.append(sample['GT'][0])
                gts_other.append(sample['GT'][2])
    return ''.join(gts_ref), ''.join(gts_other)


def is_monomorphic(genotypes):
    return genotypes == genotypes[0] * len(genotypes)


def main(vcf_filename, reference_samples):
    # Open VCF file
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))

    # Initalize counters
    count_total = count_fixed_both = 0
    count_fixed_others = count_fixed_ref = 0

    # Iterate through SNPs and evaluate reference versus the others
    for record in vcf_reader:
        gts_ref_concat, gts_other_concat = get_genotypes(record,
                                                         reference_samples)
        if gts_ref_concat:
            if (is_monomorphic(gts_other_concat) and
                    is_monomorphic(gts_ref_concat)):
                if gts_other_concat[0] not in gts_ref_concat:
                    # Alternatively fixed in both reference and others
                    count_fixed_both += 1
                else:
                    print('Error: monomorphic locus found')
            elif (is_monomorphic(gts_other_concat) and
                  len(gts_other_concat) >= 4):
                # Fixed in others
                count_fixed_others += 1
            elif (is_monomorphic(gts_ref_concat) and
                  len(gts_ref_concat) >= 4):
                # Fixed in reference
                count_fixed_ref += 1
        count_total += 1

    # Summary
    print('{0} SNPs found in VCF file.'.format(count_total))
    print(('{0} SNPs were fixed in the non-reference'
           ' samples,').format(count_fixed_both +
                               count_fixed_others))
    print(('of which {0} SNPs were alternatively fixed '
           'in the reference samples').format(count_fixed_both))
    print('and {0} SNPs were variable in the reference samples.'.format(
        count_fixed_others))
    print('{0} SNPs are fixed in reference samples'
          ' (only).'.format(count_fixed_ref))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('reference_samples', metavar='reference_samples',
                        nargs='*', help='sample(s) against which the \
                        remainder of the dataset will be compared')
    args = parser.parse_args()
    main(args.vcf_filename, args.reference_samples)
