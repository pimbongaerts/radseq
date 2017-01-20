#!/usr/bin/env python
"""
Generate list of loci that are identified in both input files 
(first column in comma-separated values file)

"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

import sys

HEADER_CHAR = '#'

def is_header(line):
	""" Return True if line is a header """
	return line[:len(HEADER_CHAR)] == HEADER_CHAR

def get_values_from_file(file):
	""" Return values if line is a header """
	value_list = []
	for line in file:
		cols = line.split(',')
		if not is_header(line):
			value_list.append(int(cols[0]))
	return value_list

def main(incl_filename1, incl_filename2, output_filename):
	# Store values of incl file in list
	incl_file1 = open(incl_filename1, 'r')
	incl_list1 = get_values_from_file(incl_file1)
	incl_file1.close()

	# Store values of incl file in list
	incl_file2 = open(incl_filename2, 'r')
	incl_list2 = get_values_from_file(incl_file2)
	incl_file2.close()
	
	# Iterate through values and only keep those present in both lists
	incl_set1 = set(incl_list1)
	incl_set2 = set(incl_list2)
	match_count = 0
	output_file = open(output_filename, 'w')
	for value in sorted(incl_set1):
	    if value in incl_set2:
	        output_file.write('{0}\n'.format(value))
	        match_count += 1
	output_file.close()
	
	# Output stats
	print('File `{0}`: {1} unique values (total: {2})'.format(incl_filename1,
	                                    len(incl_set1), len(incl_list1)))
	print('File `{0}`: {1} unique values (total: {2})'.format(incl_filename2,
	                                    len(incl_set2), len(incl_list2)))	
	print('Output file `{0}`: {1} unique values'.format(output_filename, 
	                                    match_count))

main(sys.argv[1],sys.argv[2],sys.argv[3])