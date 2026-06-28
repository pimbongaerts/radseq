#!/usr/bin/env python
"""
Batch filter VCFs with vcftools. Recursively finds every `.vcf` file (in a
directory and its subdirectories) that has a matching `_samples_to_remove.txt`
file (as produced by `vcf_missing_data.py -s`) and runs vcftools to:
  1. remove the listed individuals (--remove),
  2. keep only sites with >= 2 alleles (--min-alleles 2),
  3. remove singleton alleles (--mac 2),
  4. optionally apply a --max-missing site filter.

vcftools computes the allele-based site statistics (--mac / --max-missing) over
only the retained individuals, so the allele filters are inherently applied
after the individual removal within the single command.

For each input, output is written to `<base><postfix>.vcf` (default postfix
`_filtered`), e.g. `STEPHANOCOENIA.vcf` -> `STEPHANOCOENIA_filtered.vcf`.
"""
import os
import sys
import glob
import shutil
import argparse
import subprocess

import vcf_missing_data

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

VCF_PATTERN = '*.vcf'
VCFTOOLS_CMD = 'vcftools'
DEFAULT_POSTFIX = '_filtered'
MIN_ALLELES = '2'
MIN_MAC = '2'


def find_vcf_files(directory):
    """Recursively find all `.vcf` files in directory and subdirectories"""
    pattern = os.path.join(directory, '**', VCF_PATTERN)
    return sorted(glob.glob(pattern, recursive=True))


def filter_vcf(vcf_filename, remove_filename, output_filename, max_missing):
    """Run vcftools to remove individuals and apply allele-based filters.

    The allele filters (--min-alleles, --mac) and --max-missing are computed by
    vcftools over the retained individuals only, i.e. after the --remove step.
    """
    cmd = [VCFTOOLS_CMD,
           '--vcf', vcf_filename,
           '--remove', remove_filename,
           '--min-alleles', MIN_ALLELES,
           '--mac', MIN_MAC]
    if max_missing is not None:
        cmd += ['--max-missing', str(max_missing)]
    cmd += ['--recode', '--recode-INFO-all', '--stdout']

    with open(output_filename, 'w') as output_file:
        result = subprocess.run(cmd, stdout=output_file, stderr=sys.stderr)
    return result.returncode == 0


def main(directory, max_missing=None, postfix=DEFAULT_POSTFIX):
    # Ensure vcftools is available
    if shutil.which(VCFTOOLS_CMD) is None:
        sys.exit('Error: `{0}` not found in PATH'.format(VCFTOOLS_CMD))

    vcf_files = find_vcf_files(directory)
    if not vcf_files:
        sys.exit('No `.vcf` files found in {0}'.format(directory))

    filtered = skipped = failed = 0
    for vcf_filename in vcf_files:
        base = os.path.splitext(vcf_filename)[0]
        remove_filename = base + vcf_missing_data.OUTPUT_SUFFIX

        # Only filter VCFs that have a matching `_samples_to_remove.txt` file
        if not os.path.isfile(remove_filename):
            skipped += 1
            continue

        output_filename = base + postfix + '.vcf'
        sys.stderr.write('Filtering {0} -> {1}\n'.format(vcf_filename,
                                                         output_filename))
        if filter_vcf(vcf_filename, remove_filename, output_filename,
                      max_missing):
            filtered += 1
        else:
            failed += 1
            sys.stderr.write('Error: vcftools failed on {0}\n'.format(
                vcf_filename))
            # Remove the (likely incomplete) output file
            if os.path.isfile(output_filename):
                os.remove(output_filename)

    sys.stderr.write('Filtered {0} file(s); skipped {1} without remove file; '
                     '{2} failed\n'.format(filtered, skipped, failed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse
                                     .RawDescriptionHelpFormatter)
    parser.add_argument('directory', metavar='directory', nargs='?',
                        default=os.getcwd(),
                        help='directory to search for `.vcf` files '
                             '(default: current working directory)')
    parser.add_argument('--max-missing', type=float, default=None,
                        metavar='VALUE',
                        help='optional vcftools `--max-missing` site filter '
                             '(0-1, where 1 allows no missing data and 0 '
                             'allows completely missing sites)')
    parser.add_argument('--postfix', default=DEFAULT_POSTFIX, metavar='POSTFIX',
                        help='postfix for output filenames (default: {0}); '
                             'e.g. `--postfix _2b` turns `STEPHANOCOENIA.vcf` '
                             'into `STEPHANOCOENIA_2b.vcf`'.format(
                                 DEFAULT_POSTFIX))
    args = parser.parse_args()
    main(args.directory, args.max_missing, args.postfix)
