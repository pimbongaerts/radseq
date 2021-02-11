#!/usr/bin/env python
"""
Creates a gene set file for Gene Ontology (GO) enrichment analysis in Gowinda.
Input file is a tsv with gene_id as first column, and the list of corresponding
GO terms (with corresponding IDs in square brackets) as second column. The
input file can be created from a tsv downloaded from the UniProt portal:
https://www.uniprot.org/uniprot/
"""
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

def main(input_filename):
    # Dicts to save GO descriptions and assigned genes
    GO_genes = {}
    GO_descriptions = {}
    # Read input file
    input_file = open(input_filename, 'r')
    for line in input_file:
        # Iterate over each gene_id
        gene_id = line.split('\t')[0]
        GO_terms = line.split('\t')[1].split(';')
        # Extract GO ids and GO terms
        for GO_term in GO_terms:
            if GO_term.strip() != '':
                GO_description = GO_term.strip().split('[GO:')[0]
                GO_id_extract = GO_term.split('[GO:', 1)[1].split(']')[0]
                GO_id = 'GO:{0}'.format(GO_id_extract)
                if not GO_id in GO_genes:
                    GO_genes[GO_id] = []
                    GO_descriptions[GO_id] = GO_description.strip()
                GO_genes[GO_id].append(gene_id)
    input_file.close()
    # Output results by GO ID
    for GO_id in sorted(GO_descriptions.keys()):
        print('{0}\t{1}\t{2}'.format(GO_id,
                                     GO_descriptions[GO_id],
                                     ' '.join(GO_genes[GO_id])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_filename', metavar='input_file',
                        help='input file with gene_ids, GO IDs and GO terms')
    args = parser.parse_args()
    main(args.input_filename)
