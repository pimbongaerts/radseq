#!/usr/bin/env python
"""
Creates NJ tree from a genetic distance matrix. Outputs ASCII format to
STDOUT, a nexus-formatted tree to output file, and a PDF visualization of the
tree (same name as the output file but with a `.pdf` extension). Note:
distance matrix can be created from `vcf` using `vcf_gdmatrix.py`.
"""
import sys
import os
import argparse
import numpy as np
import Bio.Phylo
import Bio.Phylo.TreeConstruction
import tre2pdf

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def get_matrix_from_file(matrix_filename):
    """ Get genetic distance matrix from file """
    matrix_file = open(matrix_filename, 'r')
    names = []
    matrix = []
    for index, line in enumerate(matrix_file, start=1):
        if index > 1:      # Skip header
            cols = line.strip().split()
            names.append(cols[0])
            matrix.append([float(i) for i in cols[1:index]])
    matrix_file.close()
    distance_matrix = Bio.Phylo.TreeConstruction._DistanceMatrix(names=names,
                                                                 matrix=matrix)
    return distance_matrix


def main(matrix_filename, tree_output_filename):
    # Read matrix and labels from file
    distance_matrix = get_matrix_from_file(matrix_filename)

    # Generate tree from distance matrix
    constructor = Bio.Phylo.TreeConstruction.DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    tree.ladderize()

    # Output to screen
    Bio.Phylo.draw_ascii(tree)
    Bio.Phylo.write(tree, tree_output_filename, 'nexus')

    # Output PDF visualization (same name as tree file but with .pdf extension)
    pdf_output_filename = os.path.splitext(tree_output_filename)[0] + '.pdf'
    tre2pdf.main(tree_output_filename, output_filename=pdf_output_filename,
                 tree_format='nexus')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix_filename', metavar='matrix_file',
                        help='text file (tsv or csv) with genetic distance \
                        matrix')
    parser.add_argument('tree_output_filename', metavar='tree_output_file',
                        help='nexus file with output tree')
    args = parser.parse_args()
    main(args.matrix_filename, args.tree_output_filename)
