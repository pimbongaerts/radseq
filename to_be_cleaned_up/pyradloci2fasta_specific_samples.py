#!/usr/bin/env python
"""
Creates a FASTA file with a representative sequence for each locus in a pyRAD 
`.loci` file. It will only include loci that are genotyped for one or more 
indicated reference samples. The sequence of the first encountered sample will 
be used as a reference.

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

import sys

def main(inputFilename, reference_samples):	
	inputFile = open(inputFilename, 'r')
	outputFile = open(inputFilename.replace('.loci','') + '.fa', 'w')
	sequence = ''
	locuscount = 0
	
	for line in inputFile:
		if sequence == '' and line[0] == '>':
			# Only store sequence if sample in reference_samples
			cols = line.split()
			sample = cols[0][1:].strip()
			if sample in reference_samples:
			    sequence = cols[1].strip()
		elif line[0] == '/' and not sequence == '':
			# Retrieve locus number and write to output file
			cols = line.split('|')
			locus_number = cols[1].strip()
			outputFile.write('>{0}\n'.format(locus_number))
			sequence_no_gaps = sequence.replace('-','')
			outputFile.write('{0}\n'.format(sequence_no_gaps))
			locuscount += 1
			sequence = ''
			
	inputFile.close()
	outputFile.close()

main(sys.argv[1],sys.argv[2:])