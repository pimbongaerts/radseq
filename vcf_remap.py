#!/usr/bin/env python
"""
Remaps variants in VCF format to new CHROM and POS as obtained through the
`mapping_get_bwa_matches.py` scripts. Positions are rough estimates because:
(1) new position is simply an offset of the mapping position + 0-based
position in locus (and e.g. do not take into account reference insertions),
 (2) one standard contig length is used to determine pos in reverse mapping
 reads (flag 16).
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

VCF_HEADER_CHAR = '#'
STR_LENGTH = 255
C_RAD_LOCUS = 0
C_REF_CONTIG = 1
C_REF_POS = 2
C_FLAG = 3


class MappingInfo(object):
    """ Mapping info for particular locus """

    def __init__(self, ref_contig, ref_pos, flag):
        self.ref_contig = ref_contig
        self.ref_pos = int(ref_pos)
        self.flag = int(flag)

    def remap(self, pos, locus_len):
        """ Obtain remapped CHROM and POS values for locus """
        new_chrom = self.ref_contig
        if self.flag == 0:       # Forward
            return new_chrom, (self.ref_pos + pos - 1)
        elif self.flag == 16:     # Reverse
            return new_chrom, (self.ref_pos + locus_len - pos - 1)
        else:
            sys.exit('Unsupported flag: {0}'.format(flag))


def load_mapping_results(mapping_filename):
    """ Load mapping results into dictionary of RadLocus objects """
    mapping_info = {}
    mapping_file = open(mapping_filename, 'r')
    for line in mapping_file:
        cols = line.replace(',', ' ').split()
        mapping_info[cols[C_RAD_LOCUS]] = MappingInfo(cols[C_REF_CONTIG],
                                                      cols[C_REF_POS],
                                                      cols[C_FLAG])
    mapping_file.close()
    return mapping_info


def remap_vcf_file(vcf_filename, mapping_info, locus_len):
    """ Iterate over vcf file and replace with new CHROM/POS values """
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        # Output headers as they are
        if line[0] == VCF_HEADER_CHAR:
            print(line.strip())
            continue
        # Replace CHROM and POS in datarows
        cols = line.split()
        locus_name, pos = cols[0:2]
        rest_of_line = '\t'.join(cols[2:])
        new_chrom, new_pos = mapping_info[locus_name].remap(int(pos),
                                                            locus_len)
        print('{0}\t{1}\t{2}'.format(new_chrom, new_pos, rest_of_line))
    vcf_file.close()

def main(vcf_filename, mapping_filename, locus_len):
    mapping_info = load_mapping_results(mapping_filename)
    remap_vcf_file(vcf_filename, mapping_info, locus_len)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='vcf input file')
    parser.add_argument('mapping_filename', metavar='mapping_file',
                        help='file with mapping results')
    parser.add_argument('locus_len', metavar='locus_length',
                        help='length of query loci')
    args = parser.parse_args()
    main(args.vcf_filename, args.mapping_filename, int(args.locus_len))
