#!/usr/bin/env python
"""
Wrapper for PGDspider on Mac OS to convert `.vcf` files to various formats.

Note : set PGDSPIDER_PATH constant before use, and make script executable
in terminal with `$ chmod +x vcf_spider.py`.
"""
import sys
import os
import argparse
from subprocess import call
from time import sleep

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


PGDSPIDER_PATH = '/Users/pbongaerts/tools/PGDSpider_2.1.1.5/PGDSpider2-cli.jar'
PLOIDY = 'DIPLOID'

EXT_STRUCTURE = '.structure'
EXT_STRUCTURE_SHORT = '.str'
EXT_BAYESCAN = '.bayescan'
EXT_GENEPOP = '.genepop'
EXT_LOSITAN = '.lositan'
EXT_ARLEQUIN = '.arlequin'
EXT_ARLEQUIN_SHORT = '.arp'
EXT_PLINK = '.plink'
EXT_PLINK_SHORT = '.ped'


def write_spid(pop_filename, output_filename):
    """Write .spid file with conversion settings"""
    spid_file_name = '{0}/conversion.spid'.format(os.getcwd())
    spid_file = open(spid_file_name, 'w')
    # Input format
    spid_file.write('PARSER_FORMAT=VCF\n\n')
    # Popfile
    if '\\' not in pop_filename:
        pop_file_path = '{0}/{1}'.format(os.getcwd(), pop_filename)
    else:
        pop_file_path = pop_filename
    spid_file.write('VCF_PARSER_POP_QUESTION=true\n')
    spid_file.write('VCF_PARSER_POP_FILE_QUESTION=%s\n' % pop_file_path)
    # Etc
    spid_file.write('VCF_PARSER_PLOIDY_QUESTION=%s\n' % PLOIDY)
    spid_file.write('VCF_PARSER_MONOMORPHIC_QUESTION=true\n')
    spid_file.write('VCF_PARSER_PL_QUESTION=false\n')
    spid_file.write('VCF_PARSER_QUAL_QUESTION=\n')
    spid_file.write('VCF_PARSER_GTQUAL_QUESTION=\n')
    spid_file.write('VCF_PARSER_IND_QUESTION=\n')
    spid_file.write('VCF_PARSER_REGION_QUESTION=\n')
    spid_file.write('VCF_PARSER_READ_QUESTION=\n\n')

    # Output format
    filename, file_extension = os.path.splitext(output_filename)
    if file_extension in (EXT_STRUCTURE, EXT_STRUCTURE_SHORT):
        spid_file.write('WRITER_FORMAT=STRUCTURE\n\n')
        spid_file.write('STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP\n')
        spid_file.write('STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false\n')
    elif file_extension == EXT_BAYESCAN:
        spid_file.write('WRITER_FORMAT=GESTE_BAYE_SCAN\n\n')
        spid_file.write('GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP\n')
    elif file_extension in (EXT_GENEPOP, EXT_LOSITAN):
        spid_file.write('WRITER_FORMAT=GENEPOP\n\n')
        spid_file.write('GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP\n')
    elif file_extension in (EXT_ARLEQUIN, EXT_ARLEQUIN_SHORT):
        spid_file.write('WRITER_FORMAT=ARLEQUIN\n\n')
        spid_file.write('ARLEQUIN_WRITER_DATA_TYPE_QUESTION=SNP\n')
    elif file_extension in (EXT_PLINK, EXT_PLINK_SHORT):
        spid_file.write('WRITER_FORMAT=PED\n\n')
        map_file = output_filename.replace(file_extension, '.map')
        spid_file.write('PED_WRITER_MAP_FILE_QUESTION=%s\n' % map_file)
        spid_file.write('PED_WRITER_ZERO_CHAR_QUESTION=\n')
        spid_file.write('PED_WRITER_LOCUS_COMBINATION_QUESTION=\n')
        spid_file.write('PED_WRITER_MAP_QUESTION=true\n')
    else:
        sys.exit('Unrecognised file extension: {0}'.format(file_extension))
    spid_file.close()
    return spid_file_name


def main(vcf_filename, pop_filename, output_filename):
    # Write .spid file with conversion settings
    print('create .spid conversion file...')
    sys.stdout.flush()
    spid_file_name = write_spid(pop_filename, output_filename)
    # Execute PGDspider
    print('execute PGDspider...')
    sys.stdout.flush()
    vcf_path = '{0}/{1}'.format(os.getcwd(), vcf_filename)
    java_call = ('java -Xmx1024m -Xms512m -jar {0} -inputfile \"{1}\" '
                 '-outputfile \"{2}\" -spid \"{3}\"').format(PGDSPIDER_PATH,
                                                             vcf_path,
                                                             output_filename,
                                                             spid_file_name)
    call([java_call], shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_filename',
                        help='original vcf file')
    parser.add_argument('pop_filename', metavar='pop_filename',
                        help='pop filename (.txt)')
    parser.add_argument('output_filename', metavar='output_filename',
                        help='output filename (extension used to determine \
                        file format (.genepop, .bayescan, .structure or \
                        .arlequin)')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename, args.output_filename)
