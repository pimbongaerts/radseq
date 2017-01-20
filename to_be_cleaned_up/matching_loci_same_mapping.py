#!/usr/bin/env python
"""
Only output loci (first column) that are identified in both location and are 
mapping to the same reference (second column)

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'
import sys

HEADER_CHAR = '#'

def is_header(line):
	""" Return True if line is a header """
	return line[:len(HEADER_CHAR)] == HEADER_CHAR

def get_mapping_from_file(filename):
	""" Return values if line is a header """
	mapping_file = open(filename, 'r')
	mapping_dict = {}
	for line in mapping_file:
		cols = line.split(',')
		locus = cols[0].strip()
		mapping = cols[1].split()[0].strip()
		if not is_header(line):
			mapping_dict[locus] = mapping
	mapping_file.close()
	return mapping_dict

def main(incl_filename1, incl_filename2, output_filename):
	# Store values of incl file in list
	mapping_dict1 = get_mapping_from_file(incl_filename1)

	# Store values of incl file in list
	mapping_dict2 = get_mapping_from_file(incl_filename2)
	
	# Iterate through values and only keep those present in both lists
	output_file = open(output_filename, 'w')
	match_count = 0
	for locus in sorted(mapping_dict1):
	    # Check if locus occurs also in file2
	    if locus in mapping_dict2:
	        # Check if locus maps to same location in file2
	        print(mapping_dict1[locus], mapping_dict2[locus])
	        if mapping_dict1[locus] == mapping_dict2[locus]:
	            mapping = mapping_dict1[locus]
	            output_file.write('{0}\t{1}_{0}\n'.format(locus, 
	                                                      mapping))
	            match_count += 1
	output_file.close()
	
	# Output stats
	print('File `{0}`: {1} values'.format(incl_filename1,
	                                    len(mapping_dict1)))
	print('File `{0}`: {1} values'.format(incl_filename2,
	                                    len(mapping_dict2)))	
	print('Output file `{0}`: {1} unique values'.format(output_filename, 
	                                    match_count))

main(sys.argv[1],sys.argv[2],sys.argv[3])