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
			print('>{0}'.format(locus_number))
			if sequence == '':
			    final_sequence = alt_sequence #.replace('-','') to remove gaps
			else:
			    final_sequence = sequence  #.replace('-','') to remove gaps
			print('{0}'.format(final_sequence))
			locuscount += 1
			sequence = ''
			
	inputFile.close()

main(sys.argv[1],sys.argv[2:])