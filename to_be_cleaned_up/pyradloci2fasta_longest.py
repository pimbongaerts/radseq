#!/usr/bin/env python
"""
Create FASTA file with a representative sequence (using the first sample with
the longest sequence and the smallest number of gapped sites) for
each locus in pyRAD or ipyrad .loci file.

Takes .loci filename as argument.

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

import sys


def main(input_filename):
    input_file = open(input_filename, 'r')
    output_file = open(input_filename.replace('.loci', '') + '.fa', 'w')

    selected_sequence = ''

    for line in input_file:
        if line[0] == '/':
            # Retrieve locus number and write to output file
            cols = line.split('|')
            locus_number = cols[1].strip()
            output_file.write('>{0}\n'.format(locus_number))
            output_file.write('{0}\n'.format(selected_sequence))
            selected_sequence = ''
        else:
            # Select sequence if longer than previous sequence and with
            # smaller number of gaps (unless first sample)
            cols = line.split()
            current_sequence = cols[1].strip()
            if selected_sequence == '':
                selected_sequence = current_sequence
            elif (len(current_sequence) >= len(selected_sequence) and
                  current_sequence.count('-') < selected_sequence.count('-')):
                selected_sequence = current_sequence

    input_file.close()
    output_file.close()

main(sys.argv[1])
