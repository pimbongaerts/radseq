#!/usr/bin/env python
"""
Adds ordering indexes to a file with CHROM and POS columns, based on RaGOO mapping to
a chromosome-level assembly. These ordering indexes (chromosome, scaffold and SNP order)
can be used for Manhattan-style plots, without having to remap coordinates of the vcf.

CHROM and POS columns should be the first and second column in the input file. In the output, 
three extra columns are inserted (after the CHROM and POS columns) that correspond to the 
chromosome ID, scaffold order on that chromosome, and position order within the scaffold 
(based on the RaGOO mapping orientation).
"""
import sys
import argparse
import glob
import os

class Scaffold(object):
    """ Mapping info for Scaffold """

    def __init__(self, chrom, scaffold_orderid, forward_orientation):
        self.chrom = chrom
        self.scaffold_orderid = scaffold_orderid
        self.forward_orientation = forward_orientation

def read_ragoo_orderings_files(ragoo_orderings_path):
	""" Create dictionary with RaGOO mapping information for each scaffold """
	scaffolds = {}
	for filename in glob.glob(os.path.join(ragoo_orderings_path, '*.txt')):
		orderings_file = open(filename, 'r')
		chrom = int(filename.split('/')[2].split('_')[0][3:])
		scaffold_orderid = 0
		for line in orderings_file:
			scaffold_orderid += 1
			scaffold_name, forward_orientation = line.split()[0:2]
			scaffolds[scaffold_name]= Scaffold(chrom, scaffold_orderid, forward_orientation)
		orderings_file.close()
	return scaffolds

def output_scaffold_block(scaffold_block, scaffold_name, scaffolds):
	""" Output all lines belonging to a scaffold with additional index columns """
	# Process all SNPs in scaffold array
	if scaffolds[scaffold_name].forward_orientation == "+":
		snp_orderid = 0
		step = 1
	else:
		snp_orderid = len(scaffold_block) + 1
		step = -1
	for scaffold_block_line in scaffold_block:
		snp_orderid += step
		cols = scaffold_block_line.split()
		print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(cols[0], 
													cols[1], 
													scaffolds[scaffold_name].chrom,
													scaffolds[scaffold_name].scaffold_orderid,
													snp_orderid,
													'\t'.join(cols[2:])))	

def output_indexed_input_file(filename, scaffolds):
	""" Output file with additional index columns """
	input_file = open(filename, 'r')
	previous_scaffold_name = ''
	scaffold_block = []
	for line in input_file:
		if line[0] == 'C' or line[0] == "#":
			pass # do not output headers
		elif line.split()[0] == previous_scaffold_name or previous_scaffold_name == '':
			# Add all lines relating to one scaffold in an array
			scaffold_block.append(line)
			previous_scaffold_name = line.split()[0]
		else:
			output_scaffold_block(scaffold_block, previous_scaffold_name, scaffolds)
			scaffold_block = []
			scaffold_block.append(line)
			previous_scaffold_name = line.split()[0]
	# Process last scaffold block
	output_scaffold_block(scaffold_block, previous_scaffold_name, scaffolds)

def main(filename, ragoo_orderings_path):
	# Create dictionary with RaGOO mapping information for each scaffold
	scaffolds = read_ragoo_orderings_files(ragoo_orderings_path)

	# Read through input file, and append three columns (chromosome, scaffold and )
	output_indexed_input_file(filename, scaffolds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename', metavar='filename',
                        help='input file')
    parser.add_argument('ragoo_orderings_path', metavar='ragoo_orderings_path',
                        help='path with ragoo orderings files')
    args = parser.parse_args()
    main(args.filename, args.ragoo_orderings_path)