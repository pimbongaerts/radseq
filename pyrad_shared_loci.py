#!/usr/bin/env python
"""
Calculates the mean/min/max of shared loci for each sample
"""
import sys
import argparse
import numpy

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

def get_loci_for_samples(loci_filename):
  """ Obtain genotyped loci for each sample """
    sample_loci = {}
    temp_sample_list = []
    loci_file = open(loci_filename, 'r')
    for line in loci_file:
      if line[0] == '/':
        locus = line.split('|')[1]
        for sample in temp_sample_list:
          if sample not in sample_loci:
            sample_loci[sample] = []
          sample_loci[sample].append(locus)
        temp_sample_list.clear()
      else:
        sample = line.split()[0]
        temp_sample_list.append(sample)
    loci_file.close()
    return sample_loci

def main(loci_filename):
    # Obtain genotyped loci for each sample
    sample_loci = get_loci_for_samples(loci_filename)
    
    # Output header to STDOUT
    print('sample\tmax\tavg')

    # Compare genotyped loci between samples
    for sample1 in sample_loci:
      number_of_shared_loci = []
      for sample2 in sample_loci:
        intersect = len(set(sample_loci[sample1]) & set(sample_loci[sample2]))
        number_of_shared_loci.append(intersect)
      # Output results for sample
      print('{0}\t{1}\t{2}'.format(sample1, 
                                        numpy.max(number_of_shared_loci), 
                                        numpy.mean(number_of_shared_loci)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('loci_filename', metavar='loci_file',
                        help='(i)pyrad .loci file')
    args = parser.parse_args()
    main(args.loci_filename)