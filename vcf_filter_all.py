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
`_filtered`), e.g. `STEPHANOCOENIA.vcf` -> `STEPHANOCOENIA_filtered.vcf`, with a
matching `<base><postfix>.log` capturing vcftools' summary (number of
individuals and sites kept).
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


def filter_vcf(vcf_filename, remove_filename, output_filename, log_filename,
               max_missing, dry_run=False):
    """Run vcftools to remove individuals and apply allele-based filters.

    The allele filters (--min-alleles, --mac) and --max-missing are computed by
    vcftools over the retained individuals only, i.e. after the --remove step.
    vcftools' summary output (parameters and counts of individuals/sites kept)
    is captured to `log_filename`. When `dry_run` is True the vcftools command
    is printed instead of executed.
    """
    cmd = [VCFTOOLS_CMD,
           '--vcf', vcf_filename,
           '--remove', remove_filename,
           '--min-alleles', MIN_ALLELES,
           '--mac', MIN_MAC]
    if max_missing is not None:
        cmd += ['--max-missing', str(max_missing)]
    cmd += ['--recode', '--recode-INFO-all', '--stdout']

    if dry_run:
        print('{0} > {1} 2> {2}'.format(' '.join(cmd), output_filename,
                                        log_filename))
        return True

    with open(output_filename, 'w') as output_file, \
            open(log_filename, 'w') as log_file:
        result = subprocess.run(cmd, stdout=output_file, stderr=log_file)
    return result.returncode == 0


def main(directory, max_missing=None, postfix=DEFAULT_POSTFIX, dry_run=False):
    # Ensure vcftools is available (not required for a dry run)
    if not dry_run and shutil.which(VCFTOOLS_CMD) is None:
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
        log_filename = base + postfix + '.log'
        if not dry_run:
            sys.stderr.write('Filtering {0} -> {1}\n'.format(vcf_filename,
                                                             output_filename))
        if filter_vcf(vcf_filename, remove_filename, output_filename,
                      log_filename, max_missing, dry_run):
            filtered += 1
        else:
            failed += 1
            sys.stderr.write('Error: vcftools failed on {0} (see {1})\n'.format(
                vcf_filename, log_filename))
            # Remove the (likely incomplete) output VCF; keep the log for
            # debugging
            if os.path.isfile(output_filename):
                os.remove(output_filename)

    verb = 'Would filter' if dry_run else 'Filtered'
    sys.stderr.write('{0} {1} file(s); skipped {2} without remove file; '
                     '{3} failed\n'.format(verb, filtered, skipped, failed))


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
    parser.add_argument('--dry-run', action='store_true',
                        help='only print the vcftools commands that would be '
                             'run (to STDOUT) without executing them')
    args = parser.parse_args()
    main(args.directory, args.max_missing, args.postfix, args.dry_run)
