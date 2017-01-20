#!/usr/bin/env python
"""
Create FASTA file with a single sequence for each locus. Only sequences that 
are from one of the indicated reference samples are used, except for when 
none of these samples are genotyped (then the last encountered one is used).

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

import sys

def main(inputFilename, reference_samples):	
	inputFile = open(inputFilename, 'r')
	outputFile = open(inputFilename.replace('.loci','') + '.fa', 'w')
	sequence = alt_sequence = ''
	locuscount = 0
	
	for line in inputFile:
		if sequence == '' and line[0] == '>':
			# Only store sequence if sample in reference_samples
			cols = line.split()
			sample = cols[0][1:].strip()
			if sample in reference_samples:
			    sequence = cols[1].strip()
			alt_sequence = cols[1].strip()
		elif line[0] == '/':
			# Retrieve locus number and write to output file
			cols = line.split('|')
			locus_number = cols[1].strip()
			outputFile.write('>{0}\n'.format(locus_number))
			if sequence == '':
			    sequence_no_gaps = alt_sequence.replace('-','')
			else:
			    sequence_no_gaps = sequence.replace('-','')
			outputFile.write('{0}\n'.format(sequence_no_gaps))
			locuscount += 1
			sequence = ''
			
	inputFile.close()
	outputFile.close()

main(sys.argv[1],sys.argv[2:])