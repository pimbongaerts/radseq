#!/usr/bin/env python
"""
Creates a list of GO terms for each annotated gene.
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

def create_uniprot2go_dict(uniprot2go_filename):
    """" Creates look-up dict with uniprot ID (key) and GO terms (value) """ 
    uniprot2go = {}
    uniprot2go_file = open(uniprot2go_filename, 'r')
    for line in uniprot2go_file:
        cols = line.rstrip().split('\t')
        if len(cols) > 1:
            uniprot2go[cols[0]] = cols[1]
    uniprot2go_file.close()
    return uniprot2go

def main(gene2uniprot_filename, uniprot2go_filename):
    # Create look-up table for uniprot -> go_terms
    uniprot2go = create_uniprot2go_dict(uniprot2go_filename)
    # Iterate over gene2uniprot file and add output gene2go
    gene2uniprot_file = open(gene2uniprot_filename, 'r')
    last_outputted_gene_id = ''
    match_found = False
    for line in gene2uniprot_file:
        cols = line.rstrip().split('\t')
        gene_id, uniprot_id = cols[0:2]
        if  ((gene_id != last_outputted_gene_id) and 
             (uniprot_id in uniprot2go.keys())):
            print('{0}\t{1}'.format(gene_id, uniprot2go[uniprot_id]))
            last_outputted_gene_id = gene_id
    gene2uniprot_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene2uniprot_filename', 
                        metavar='gene2uniprot_filename',
                        help='input file (tsv) with the custom gene ids (first\
                        column) and the corresponding UniProt IDs (second\
                        column); when multiple UniProt IDs are given for each\
                        gene, they should be sorted by highest match)')
    parser.add_argument('uniprot2go_filename', 
                        metavar='uniprot2go_filename',
                        help='input file (tsv) with UniProt gene ids (first\
                        column) and the corresponding GO terms (second\
                        column) separated by semi-colons (;)')
    args = parser.parse_args()
    main(args.gene2uniprot_filename, args.uniprot2go_filename)
