#!/usr/bin/env python
"""
Reduce an ipyrad `.loci` file to just the loci that are present in a `.vcf`
(i.e. the loci that still contributed at least one SNP to the VCF, e.g. after
filtering). Each locus block in the `.loci` file is kept or dropped as a whole,
so the output is a valid `.loci` file containing only the retained loci.

Loci are matched by their integer locus id. ipyrad writes that id into the VCF
`ID` column as `loc<N>_pos<M>` (and into each `.loci` separator line as
`// ... |<N>...|`). When the VCF has no `loc<N>` ids (e.g. a de novo VCF whose
`CHROM` is the locus number) the `CHROM` column is used as the fallback key.

By default the reduced file is named after the VCF with a `.loci` extension (so
`--vcf FAVIINAE_filtered.vcf` writes `FAVIINAE_filtered.loci`, ready to pair with
that VCF); if that name already exists it falls back to `<vcf>_invcf.loci`.

Example:
  python3 ipyrad_reduce_loci_to_match_vcf.py --vcf FAVIINAE_filtered.vcf \
      --loci FAVIINAE.loci
"""
import os
import re
import sys
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2026 Pim Bongaerts'
__license__ = 'GPL'

LOC_RE = re.compile(r'loc(\d+)')       # ipyrad VCF ID column: loc<N>_pos<M>
INT_RE = re.compile(r'(\d+)')          # fallback: first integer in CHROM
SEP_ID_RE = re.compile(r'\s*(\d+)')    # leading id in the `.loci` separator id


def vcf_locus_ids(vcf_filename):
    """ Return (set_of_locus_ids, source_description) present in the VCF. Uses
    the `loc<N>` id in the ID column; falls back to the CHROM integer if the VCF
    carries no `loc<N>` ids. """
    ids = set()
    with open(vcf_filename, 'r') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            cols = line.split('\t', 3)
            if len(cols) >= 3:
                match = LOC_RE.search(cols[2])
                if match:
                    ids.add(int(match.group(1)))
    if ids:
        return ids, 'ID column (loc<N>)'
    with open(vcf_filename, 'r') as handle:      # de novo fallback: CHROM = id
        for line in handle:
            if not line or line[0] == '#':
                continue
            match = INT_RE.search(line.split('\t', 1)[0])
            if match:
                ids.add(int(match.group(1)))
    return ids, 'CHROM column (integer)'


def locus_id_from_separator(line):
    """ Parse the integer locus id from a `.loci` `//` separator line. The id is
    the first field after the SNP-marker section (which never contains `|`). """
    parts = line.split('|', 2)
    if len(parts) < 2:
        return None
    match = SEP_ID_RE.match(parts[1])
    return int(match.group(1)) if match else None


def reduce_loci(loci_filename, keep_ids, out_filename):
    """ Stream the `.loci` file, writing only the locus blocks whose id is in
    `keep_ids`. Returns (kept, total). """
    kept = total = 0
    block = []
    with open(loci_filename, 'r') as fin, open(out_filename, 'w') as fout:
        for line in fin:
            block.append(line)
            if line.startswith('//'):
                total += 1
                locus_id = locus_id_from_separator(line)
                if locus_id is not None and locus_id in keep_ids:
                    fout.writelines(block)
                    kept += 1
                block = []
    return kept, total


def main(vcf_filename, loci_filename, out_filename):
    if not os.path.isfile(vcf_filename):
        sys.exit('Error: vcf file `{0}` not found.'.format(vcf_filename))
    if not os.path.isfile(loci_filename):
        sys.exit('Error: loci file `{0}` not found.'.format(loci_filename))
    if not out_filename:
        base = os.path.splitext(vcf_filename)[0]
        out_filename = base + '.loci'
        if os.path.exists(out_filename):        # don't clobber an existing file
            out_filename = base + '_invcf.loci'
    if os.path.abspath(out_filename) == os.path.abspath(loci_filename):
        sys.exit('Error: output would overwrite the input loci file; use -o.')

    keep_ids, source = vcf_locus_ids(vcf_filename)
    if not keep_ids:
        sys.exit('Error: no locus ids found in `{0}`.'.format(vcf_filename))
    print('VCF loci: {0} distinct locus id(s) [{1}]'.format(
        len(keep_ids), source))
    kept, total = reduce_loci(loci_filename, keep_ids, out_filename)
    print('Loci file: {0} of {1} loci kept ({2} dropped).'.format(
        kept, total, total - kept))
    print('Reduced loci written to `{0}`.'.format(out_filename))
    if kept == 0:
        sys.stderr.write('Warning: no loci matched; check that the VCF and loci '
                         'file come from the same ipyrad assembly.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='vcf_file',
                        required=True, help='input `.vcf` file (defines which '
                        'loci to keep)')
    parser.add_argument('-l', '--loci', dest='loci_filename',
                        metavar='loci_file', required=True,
                        help='input ipyrad `.loci` file to reduce')
    parser.add_argument('-o', '--output', dest='out_filename',
                        metavar='loci_file', default=None,
                        help='output `.loci` file (default: the vcf name with a '
                             '`.loci` extension, e.g. FAVIINAE_filtered.vcf -> '
                             'FAVIINAE_filtered.loci; falls back to '
                             '`<vcf>_invcf.loci` if that already exists)')
    args = parser.parse_args()
    main(args.vcf_filename, args.loci_filename, args.out_filename)
