#!/usr/bin/env python
"""
Extracts a list of loci that have a blastn e-value below a certain threshold,
and outputs the (first) matching reference locus, as well as the alignment
length, nident, e-value and bitscore. It also compiles a set of all tax_ids,
which it uses to connect with the NCBI taxonomy database to get phylum ids
for each match using Entrez. Results are outputted to file with the chosen
e-value as post-fix, and STDOUT gives minimum alignment stats for filtered
loci.

Note: fields in input file should be (in this order): query id, subject id,
alignment length, identity, perc. identity, evalue, bitscore, staxids, stitle.
This can be achieved by using `blastn` with the following argument: `-outfmt 7
qseqid sseqid length nident pident evalue bitscore staxids stitle`.
"""
import sys
import argparse
from Bio import Entrez

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

COMMENT_CHARACTER = '#'
COL_QSEQID = 0
COL_SSEQID = 1
COL_LENGTH = 2
COL_NIDENT = 3
COL_PIDENT = 4
COL_EVALUE = 5
COL_BITSCORE = 6
COL_TAXIDS = 7
COL_TITLE = 8
TAX_ID_SEP = ';'
NCBI_SCINAME = 'ScientificName'


def main(input_filename, max_evalue, email):

    # Iterate through blastn file and output matches to list
    print('Parsing matches from {0}...'.format(input_filename))
    input_file = open(input_filename, 'r')
    output_list = []
    taxids = set()
    # Blastn stats
    min_length = min_nident = min_pident = match_count = 0
    # Iterate through blastn file
    for line in input_file:
        if line[:len(COMMENT_CHARACTER)] == COMMENT_CHARACTER:
            match_flag = False
        elif not match_flag:
            cols = line.split('\t')
            if float(cols[COL_EVALUE]) <= float(max_evalue):
                match_flag = True
                # Only retrieve first taxid when multiple given
                if TAX_ID_SEP in cols[COL_TAXIDS]:
                    tax_id = int(cols[COL_TAXIDS].split(TAX_ID_SEP)[0])
                else:
                    tax_id = int(cols[COL_TAXIDS])
                taxids.add(tax_id)
                # Write to output list if below e-value threshold
                output_list.append(('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'
                                    '\t{6}').format(cols[COL_QSEQID],
                                                    cols[COL_SSEQID],
                                                    cols[COL_LENGTH],
                                                    cols[COL_NIDENT],
                                                    cols[COL_EVALUE],
                                                    cols[COL_BITSCORE],
                                                    tax_id))

                # Keep track of min. length, nident and pident values
                if min_length == 0 or float(cols[COL_LENGTH]) < min_length:
                    min_length = float(cols[COL_LENGTH])
                if min_nident == 0 or float(cols[COL_NIDENT]) < min_nident:
                    min_nident = float(cols[COL_NIDENT])
                if min_pident == 0 or float(cols[COL_PIDENT]) < min_pident:
                    min_pident = float(cols[COL_PIDENT])
                match_count += 1
    input_file.close()

    # Look up phylum for each taxid
    print('Looking up phyla for {0} unique taxa...'.format(len(taxids)))
    Entrez.email = email
    tax_phyla = {}
    not_found = []
    for taxid in taxids:
        record = Entrez.read(Entrez.efetch(id=int(taxid),
                                           db="taxonomy"))
        phylum = [d for d in record[0]['LineageEx'] if d['Rank'] in ['phylum']]
        if phylum:
            tax_phyla[taxid] = phylum[0][NCBI_SCINAME]
        else:
            tax_phyla[taxid] = 'NOT_FOUND'
        print(taxid, tax_phyla[taxid])

    # Output to file
    print('Outputting matches to file...'.format(len(taxids)))
    output_filename = '{0}_match{1}.txt'.format(
        input_filename.replace('.txt', ''),
        max_evalue)
    output_file = open(output_filename, 'w')
    for line in output_list:
        taxid = int(line.split()[6])
        output_file.write('{0}\t{1}\n'.format(line, tax_phyla[taxid]))
    output_file.close()

    # Output match summary
    print(('\nMatches: {0} | Min.length: {1} bp | Min. nident: {2} bp | '
           'Min. pident: {3} %').format(match_count, min_length,
                                        min_nident, min_pident))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('blastn_filename', metavar='blastn_file',
                        help='blastn output file with the following fields \
                        (in that order): query id, subject id, alignment \
                        length, identity, perc. identity, evalue, bitscore')
    parser.add_argument('evalue_cut_off', metavar='evalue_cut_off', type=float,
                        help='maximum e-value for match to be included')
    parser.add_argument('email', metavar='email',
                        help='email address to be used for NCBI connection')
    args = parser.parse_args()
    main(args.blastn_filename, args.evalue_cut_off, args.email)
