#!/usr/bin/env python
"""
Batch-run `vcf_cluster_explore.py` over every `.vcf` found recursively under a
directory. For each VCF a sibling `.loci` file with the same basename is used
automatically when present, and the CSV + PDF outputs are written with their
default names into the folder that contains the VCF (so each dataset's results
sit next to its input).

Any `vcf_cluster_explore.py` option that is not file-specific is exposed here and
applied to every dataset (e.g. `--pops-from-sample-id`, `--fields`, `--method`,
`--max-k`, `--tree`, `--ordination`, `--auto-clone`, `--no-pdf`). The per-file
options (`--vcf`, `--loci`, `--output`, `--pop`, `--pdf-output`) are handled
per-VCF by this wrapper and are intentionally not forwarded.

Each VCF is run in its own subprocess, so one failing dataset does not stop the
batch; a summary of successes/failures is printed at the end.

Example:
  python3 vcf_cluster_explore_all.py datasets/ --pops-from-sample-id \
      --fields species,location,depth
"""
import os
import sys
import glob
import argparse
import subprocess

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2026 Pim Bongaerts'
__license__ = 'GPL'

SCRIPT = 'vcf_cluster_explore.py'
LOCI_EXT = '.loci'


def find_vcfs(root, pattern):
    """ Recursively find VCFs matching `pattern` under `root` (sorted). """
    return sorted(glob.glob(os.path.join(root, '**', pattern), recursive=True))


def loci_for(vcf_filename):
    """ Return the sibling `.loci` file (same basename) if it exists, else None. """
    candidate = os.path.splitext(vcf_filename)[0] + LOCI_EXT
    return candidate if os.path.isfile(candidate) else None


def build_forwarded_args(args):
    """ Reconstruct the cross-dataset `vcf_cluster_explore.py` flags to forward.
    Only options that are meaningful across every dataset are passed on. """
    fwd = ['--method', args.method,
           '--max-k', str(args.max_k),
           '--min-cluster-size', str(args.min_cluster_size),
           '--tree', args.tree_mode,
           '--ordination', args.ordination]
    if args.auto_clone:
        fwd.append('--auto-clone')
    if args.clone_threshold is not None:
        fwd += ['--clone-threshold', str(args.clone_threshold)]
    if args.force_k is not None:
        fwd += ['--k', str(args.force_k)]
    if args.pops_from_sample_id:
        fwd.append('--pops-from-sample-id')
    if args.field_names:
        fwd += ['--fields', args.field_names]
    if args.no_pdf:
        fwd.append('--no-pdf')
    return fwd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('root', nargs='?', default='.',
                        help='directory to search recursively for VCFs '
                             '(default: current directory)')
    parser.add_argument('--pattern', dest='pattern', default='*.vcf',
                        metavar='glob',
                        help='filename glob for VCFs (default: *.vcf)')
    parser.add_argument('--script', dest='script', default=None,
                        metavar='path',
                        help='path to vcf_cluster_explore.py (default: the copy '
                             'alongside this wrapper)')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
                        help='list the VCFs and the commands without running '
                             'them')
    # ---- forwarded (cross-dataset) vcf_cluster_explore.py options ----
    parser.add_argument('-m', '--method', dest='method', default='ibs',
                        choices=('ibs', 'het-masked', 'dosage', 'single-read'),
                        help='similarity measure (default: ibs)')
    parser.add_argument('--auto-clone', dest='auto_clone', action='store_true',
                        help='clone-correct using an auto-inferred threshold')
    parser.add_argument('--clone-threshold', dest='clone_threshold',
                        default=None, metavar='pct',
                        help='clone-correct using this manual similarity %% '
                             'threshold')
    parser.add_argument('--max-k', dest='max_k', type=int, default=10,
                        metavar='K', help='maximum K to evaluate (default: 10)')
    parser.add_argument('--min-cluster-size', dest='min_cluster_size', type=int,
                        default=2, metavar='N',
                        help='clusters smaller than this cannot be chosen as '
                             'best K (default: 2)')
    parser.add_argument('--k', dest='force_k', type=int, default=None,
                        metavar='K', help='force which K drives the '
                        'differentiation panels (default: best-supported K)')
    parser.add_argument('--tree', dest='tree_mode', default='upgma',
                        choices=('upgma', 'nj'),
                        help='tree to draw (default: upgma)')
    parser.add_argument('--ordination', dest='ordination', default='pca',
                        choices=('pca', 'pcoa'),
                        help='ordination for the scatter panels (default: pca)')
    parser.add_argument('--pops-from-sample-id', dest='pops_from_sample_id',
                        action='store_true',
                        help='derive annotation tracks from the sample name '
                             'fields (see vcf_cluster_explore.py)')
    parser.add_argument('--fields', dest='field_names', metavar='names',
                        default=None,
                        help='comma-separated names for the annotation tracks, '
                             'e.g. "species,location,depth"')
    parser.add_argument('--no-pdf', dest='no_pdf', action='store_true',
                        help='do not generate the PDF reports (text only)')
    args = parser.parse_args()

    script = args.script or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), SCRIPT)
    if not os.path.isfile(script):
        sys.exit('Error: cannot find `{0}` (use --script).'.format(script))

    vcfs = find_vcfs(args.root, args.pattern)
    if not vcfs:
        sys.exit('No VCFs matching `{0}` under `{1}`.'.format(
            args.pattern, args.root))
    forwarded = build_forwarded_args(args)
    print('Found {0} VCF(s) under `{1}`.'.format(len(vcfs), args.root))

    n_ok = n_fail = 0
    for i, vcf in enumerate(vcfs, 1):
        loci = loci_for(vcf)
        cmd = [sys.executable, script, '--vcf', vcf]
        if loci:
            cmd += ['--loci', loci]
        cmd += forwarded
        print('\n[{0}/{1}] {2}{3}'.format(
            i, len(vcfs), vcf, ' (+{0})'.format(os.path.basename(loci))
            if loci else ''))
        if args.dry_run:
            print('  ' + ' '.join(cmd))
            continue
        result = subprocess.run(cmd)
        if result.returncode == 0:
            n_ok += 1
        else:
            n_fail += 1
            sys.stderr.write('  FAILED (exit {0}): {1}\n'.format(
                result.returncode, vcf))

    if not args.dry_run:
        print('\nDone: {0} succeeded, {1} failed, {2} total.'.format(
            n_ok, n_fail, len(vcfs)))
        if n_fail:
            sys.exit(1)


if __name__ == '__main__':
    main()
