#!/usr/bin/env python
"""
Create FASTA file with a representative sequence (using first sample) or
all sequences (when --all_seqs flag is set) for each locus in pyRAD/ipyrad
`.loci` or `.allele` file.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(input_filename, all_seqs):
    input_file = open(input_filename, 'r')
    current_seqlength = 0
    current_sequence = ''
    sequences_in_loci = {}

    for line in input_file:
        if line[0] == '/':
            # Retrieve locus number
            cols = line.split('|')
            locus_number = cols[1].strip()
            # Write sequence(s) to output file
            if all_seqs:
                for sample in sorted(sequences_in_loci):
                    print('>{0}_{1}'.format(locus_number, sample))
                    print('{0}'.format(sequences_in_loci[sample]))                    
                sequences_in_loci.clear()
            else:
                print('>{0}'.format(locus_number))
                print('{0}'.format(current_sequence))
                current_sequence = ''
        else:
            # Retrieve sequence
            cols = line.split()
            sequence = cols[1].strip()
            if all_seqs:
                # Store all sequences in dict
                sample = cols[0].replace('>', '')
                sequences_in_loci[sample] = sequence
            elif current_sequence == '':
                # Store first sequence in locus
                current_sequence = sequence

    input_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_filename', metavar='pyrad_file',
                        help='PyRAD allele file (`.loci` or `.allele`)')
    parser.add_argument('-a', '--all_seqs', action='store_true',
                        help='set flag to output all sequences')
    args = parser.parse_args()
    main(args.input_filename, args.all_seqs)
