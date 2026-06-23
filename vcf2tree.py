#!/usr/bin/env python
"""
vcf2tree.py: Wrapper that runs vcf_gdmatrix.py and gdmatrix2tree.py to produce
a NJ tree directly from a VCF file.
"""
import os
import sys
import argparse
import subprocess

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def main(vcf_filename, pop_filename=None, min_threshold=100, force=False,
         color_field=-1):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    basename = os.path.splitext(vcf_filename)[0]
    matrix_filename = '{}_gd.txt'.format(basename)
    tree_filename = '{}_gd.tre'.format(basename)
    pdf_filename = '{}_gd.pdf'.format(basename)

    # Step 1: VCF -> distance matrix
    gdmatrix_cmd = [sys.executable,
                    os.path.join(script_dir, 'vcf_gdmatrix.py'),
                    vcf_filename]
    if pop_filename:
        gdmatrix_cmd.append(pop_filename)
    gdmatrix_cmd.extend(['--min-threshold', str(min_threshold)])
    if force:
        gdmatrix_cmd.append('--force')

    sys.stderr.write('Calculating genetic distances...\n')
    with open(matrix_filename, 'w') as f:
        result = subprocess.run(gdmatrix_cmd, stdout=f, stderr=sys.stderr)
    if result.returncode != 0:
        sys.exit('Error: vcf_gdmatrix.py failed')
    sys.stderr.write('Distance matrix written to {}\n'.format(matrix_filename))

    # Step 2: distance matrix -> NJ tree
    tree_cmd = [sys.executable,
                os.path.join(script_dir, 'gdmatrix2tree.py'),
                matrix_filename, tree_filename]

    sys.stderr.write('Building NJ tree...\n')
    result = subprocess.run(tree_cmd)
    if result.returncode != 0:
        sys.exit('Error: gdmatrix2tree.py failed')
    sys.stderr.write('Tree written to {}\n'.format(tree_filename))

    # Step 3: tree -> PDF
    pdf_cmd = [sys.executable,
               os.path.join(script_dir, 'tre2pdf.py'),
               tree_filename, '-o', pdf_filename,
               '--color_field', str(color_field)]

    sys.stderr.write('Generating tree PDF...\n')
    result = subprocess.run(pdf_cmd, stderr=sys.stderr)
    if result.returncode != 0:
        sys.exit('Error: tre2pdf.py failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf_filename', metavar='vcf_file',
                        help='input file with SNP data (`.vcf`)')
    parser.add_argument('pop_filename', metavar='pop_file',
                        nargs='?', default=None,
                        help='optional text file (tsv or csv) with individuals '
                        'and populations (controls output order)')
    parser.add_argument('--min-threshold', type=int, default=100,
                        help='minimum number of shared SNPs without missing '
                        'data (default: 100)')
    parser.add_argument('--force', action='store_true',
                        help='iteratively remove samples with fewest genotyped '
                        'loci when pairs fall below threshold')
    parser.add_argument('--color_field', type=int, default=-1,
                        help='which underscore-separated field in the sample '
                        'name to use for color coding in the PDF '
                        '(default: -1, i.e. last field)')
    args = parser.parse_args()
    main(args.vcf_filename, args.pop_filename,
         min_threshold=args.min_threshold, force=args.force,
         color_field=args.color_field)
