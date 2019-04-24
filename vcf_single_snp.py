#!/usr/bin/env python
"""
Reduces `.vcf` file to a single 'random' SNP per locus/chrom. Use for analyses
that require SNPs that are not physically linked. (although note that they of
course still may be - particularly so when dealing with short loci)
"""
import sys
import argparse
import random

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016-18 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"
COL_CHROM = 0
COL_POS = 1


def main(vcf_filename, distance_threshold):
    temp_SNPs = []; previous_CHROM = ''; previous_POS = -1

    input_file = open(vcf_filename, 'r')
    for line in input_file:
        if line[0] == HEADER_CHAR:
            print('{0}'.format(line), end='')
        else:
            cols = line.split()
            current_CHROM = cols[COL_CHROM]
            current_POS = int(cols[COL_POS])

            if (previous_CHROM not in ('', current_CHROM) or 
               (previous_POS != -1 and
                  current_POS >= (previous_POS + distance_threshold))):
                # When reaching new CHROM or a position >= than 
                # DISTANCE_THRESHOLD away from previous POS:
                # select one random SNP from list for each replicate
                print(temp_SNPs[random.randint(0, len(temp_SNPs) - 1)])
                temp_SNPs = []
            # Generate list of all SNPs in CHROM/block
            temp_SNPs.append(line.rstrip())
            previous_CHROM = current_CHROM
            previous_POS = current_POS
    
    input_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('-d', '--distance', dest='distance_threshold',
                        metavar='distance_threshold',
                        default = 2500,
                        help='optional custom distance threshold between SNPs \
                              (default is 2500; not relevant for short \
                              de-novo loci not mapped to reference scaffolds)')
    args = parser.parse_args()
    main(args.vcf_filename, int(args.distance_threshold))