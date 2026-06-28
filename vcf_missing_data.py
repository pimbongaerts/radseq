#!/usr/bin/env python
"""
Outputs list of missing data (# and % of SNPs) for each sample in VCF, to
identify poor-performing samples to eliminate prior to SNP filtering.

Takes vcf_filename as argument. Outputs to STDOUT (no output file).

"""
import sys
import os
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
MISSING_CHAR = "."
OUTPUT_HEADER = 'INDIVIDUAL\tMISS\tGENO\tTOTAL\t% GENOTYPED'
OUTPUT_SUFFIX = '_samples_to_remove.txt'


def main(vcf_filename, lowest_n=None, threshold=None, save_to_file=False):
    # Read in genotypes for all individuals
    individuals = {}
    genotypes = {}
    with open(vcf_filename, 'r') as vcf_file:
        for line in vcf_file:
            line = line.strip()
            if line[0:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
                cols = line.split('\t')
                # Store individual names with col_index as key
                for col_index, col in enumerate(cols):
                    if col_index >= 9:
                        individual_name = col
                        individuals[col_index] = individual_name
                        genotypes[individual_name] = []
            elif not line[0:len(HEADER_CHAR)] == HEADER_CHAR:
                cols = line.split('\t')
                # Store genotypes for each individual
                for col_index, col in enumerate(cols):
                    if col_index >= 9:
                        genotype = col
                        individual = individuals[col_index]
                        genotypes[individual].append(genotype[0:3])

    # Assess missing data for each individual
    records = []
    for individual in sorted(genotypes.keys()):
        # A genotype is counted as missing if either allele is missing
        # (e.g. `./.` or partial calls such as `./0`)
        missing_count = sum(1 for genotype in genotypes[individual]
                            if MISSING_CHAR in genotype)
        total_count = len(genotypes[individual])
        genotyped_count = total_count - missing_count
        if total_count > 0:
            perc_count = round((genotyped_count / total_count) * 100, 2)
        else:
            perc_count = 'NA'
        records.append((individual, missing_count, genotyped_count,
                        total_count, perc_count))

    # Numeric sort key for % genotyped; samples with no loci ('NA') rank lowest
    def perc_key(record):
        return float('-inf') if record[4] == 'NA' else record[4]

    # If requested, only keep samples below the threshold (% genotyped)
    if threshold is not None:
        records = [record for record in records if perc_key(record) < threshold]

    # Sort by % genotyped ascending (worst-performing samples first) when
    # limiting to the n lowest or producing a threshold-based remove list
    if lowest_n is not None or threshold is not None:
        records.sort(key=perc_key)
    if lowest_n is not None:
        records = records[:lowest_n]

    # Assemble output lines
    if threshold is not None:
        # Bare sample names (no header) for use as a vcftools --remove file
        output_lines = [record[0] for record in records]
    else:
        output_lines = [OUTPUT_HEADER]
        output_lines += ['{0}\t{1}\t{2}\t{3}\t{4}'.format(*record)
                         for record in records]

    # Output to STDOUT, or to a `<vcf_basename>_samples_to_remove.txt` file
    if save_to_file:
        output_filename = os.path.splitext(vcf_filename)[0] + OUTPUT_SUFFIX
        with open(output_filename, 'w') as output_file:
            output_file.write('\n'.join(output_lines) + '\n')
        sys.stderr.write('Output written to {0}\n'.format(output_filename))
    else:
        print('\n'.join(output_lines))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-n', type=int, default=None, metavar='N',
                        help='only output the N lowest records (by %% '
                             'genotyped), i.e. the worst-performing samples')
    parser.add_argument('-t', type=float, default=None, metavar='THRESHOLD',
                        help='only output names of samples below this %% '
                             'genotyped threshold, one per line and without '
                             'header (e.g. for a vcftools `--remove` file)')
    parser.add_argument('-s', action='store_true',
                        help='save output to a file named after the VCF with '
                             '`{0}` appended (e.g. `STEPHANOCOENIA.vcf` -> '
                             '`STEPHANOCOENIA{0}`) instead of writing to '
                             'STDOUT'.format(OUTPUT_SUFFIX))
    args = parser.parse_args()
    main(args.vcf_filename, args.n, args.t, args.s)
