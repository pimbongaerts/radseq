#!/usr/bin/env python
"""
Calculates Genetic Distance (Hamming / p-distance) for each pair of individuals
in a `.vcf` file and outputs as matrix. Popfile can optionally be supplied to
indicate order in matrix (otherwise VCF sample order is used).

The distance is the allele-sharing `1 - IBS` over sites genotyped in both
individuals (pairwise-complete), where heterozygote-vs-heterozygote counts as a
full match. This is the L1 / allele-count family and is equivalent to PLINK
`--distance 1-ibs flat-missing` (verified identical to rounding) and to PLINK's
default `--distance` / `allele-ct` output (same metric, rescaled). It is the
same measure as the default `ibs` method in `vcf_clone_detect.py`; see that
script for alternative measures (dosage = PLINK `--make-rel cov`; single-read =
ANGSD `-doIBS 1`). Note PLINK/ANGSD apply their own missing-data scaling, so
values can drift slightly on data with missing genotypes.
"""
import sys
import argparse
import operator
import itertools
from collections import Counter

__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2016 Pim Bongaerts"
__license__ = "GPL"


VCF_HEADER_CHAR = "#"
VCF_INDIV_HEADER = "#CHROM"
ERROR_CODE = -1
SELF_MATCH = 0

SEPARATOR = "\t"
DEFAULT_MIN_THRESHOLD = 100


def get_genotypes_from_vcf(vcf_filename):
    """Read genotypes from VCF into dict of lists (indivs: genotypes)"""
    indivs_gts = {}
    colnrs_indivs = {}
    # Iterate through VCF and store genotypes
    vcf_file = open(vcf_filename, "r")
    for line in vcf_file:
        cols = line.rstrip().split()
        if line[0 : len(VCF_INDIV_HEADER)] == VCF_INDIV_HEADER:
            # Extract indiv names from header
            colnrs_indivs = {i: indiv for (i, indiv) in enumerate(cols[9:], start=9)}
            # Create dict item with list for indiv
            indivs_gts = {indiv: [] for indiv in cols[9:]}
        elif line[0] != VCF_HEADER_CHAR:
            # Extract genotypes for all individuals (current SNP)
            for index, genotype in enumerate(cols[9:], start=9):
                indivs_gts[colnrs_indivs[index]].append(genotype)
    vcf_file.close()
    return indivs_gts


def count_genotyped_loci(indivs_gts):
    """Count non-missing genotypes per individual"""
    return {
        indiv: sum(1 for gt in gts if gt[:1] != ".")
        for indiv, gts in indivs_gts.items()
    }


def compare_indivs(indiv1_gts, indiv2_gts, min_threshold):
    """Calculate genetic distance between 2 individuals"""
    matches_count = total_count = 0

    # Evaluate match across each individual SNP
    for x in range(0, len(indiv1_gts)):
        genotype1 = indiv1_gts[x][:3]
        genotype2 = indiv2_gts[x][:3]

        # Only consider SNPs that are genotyped for both individuals
        if genotype1[0] != "." and genotype2[0] != ".":
            if genotype1 == genotype2:
                matches_count += 1
            elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2]:
                matches_count += 0.5
            total_count += 1

    if total_count >= min_threshold:
        return round(1 - (matches_count / total_count), 4)
    else:
        return ERROR_CODE


def get_comparison_key(pairs_gds, indiv1, indiv2):
    """Obtain comparison key from dict for pair of individuals"""
    comparison_key1 = "{0}:{1}".format(indiv1, indiv2)
    comparison_key2 = "{0}:{1}".format(indiv2, indiv1)
    if comparison_key1 in pairs_gds:
        return comparison_key1
    elif comparison_key2 in pairs_gds:
        return comparison_key2
    else:
        return sys.exit("Unexpected error")


def get_match_score(pairs_gds, indiv1, indiv2):
    """Obtain comparison key from dict for pair of individuals"""
    if indiv1 == indiv2:
        return str(SELF_MATCH)
    else:
        comparison_key = get_comparison_key(pairs_gds, indiv1, indiv2)
        return str(pairs_gds[comparison_key])


def main(
    vcf_filename,
    pop_filename=None,
    min_threshold=DEFAULT_MIN_THRESHOLD,
    force=False,
    output=None,
):
    # Read genotypes for all individual and loci into a dict
    indivs_gts = get_genotypes_from_vcf(vcf_filename)

    # Read order of individuals from popfile, or use VCF order
    if pop_filename:
        with open(pop_filename) as file:
            indiv_order = [
                line.replace(",", " ").split()[0] for line in file if line.strip()
            ]
        # Ensure all popfile samples are present in the VCF
        missing = [indiv for indiv in indiv_order if indiv not in indivs_gts]
        if missing:
            sys.exit(
                "Error: sample(s) in pop_file not found in VCF: "
                "{0}".format(", ".join(missing))
            )
    else:
        indiv_order = list(indivs_gts.keys())

    removed_samples = []

    # Iteratively compute pairwise distances, removing low-data samples
    # when --force is enabled and a pair falls below threshold
    while True:
        n = len(indivs_gts)
        total_pairs = n * (n - 1) // 2
        pairs_gds = {}
        restart = False

        for i, (indiv1, indiv2) in enumerate(itertools.combinations(indivs_gts, 2), 1):
            comparison_key = "{0}:{1}".format(indiv1, indiv2)
            pairs_gds[comparison_key] = compare_indivs(
                indivs_gts[indiv1], indivs_gts[indiv2], min_threshold
            )

            if pairs_gds[comparison_key] == ERROR_CODE:
                if not force:
                    sys.exit(
                        "Error: less than {0} SNPs shared between "
                        "{1} and {2}".format(min_threshold, indiv1, indiv2)
                    )

                # Remove the worse-genotyped sample of the failing pair, so the
                # removal actually targets the comparison that fell below threshold
                gt_counts = count_genotyped_loci(indivs_gts)
                worst = min((indiv1, indiv2), key=lambda x: gt_counts[x])
                sys.stderr.write(
                    "\rRemoving {0} (genotyped loci: {1}/{2})\n".format(
                        worst, gt_counts[worst], len(indivs_gts[worst])
                    )
                )
                removed_samples.append(worst)
                del indivs_gts[worst]
                indiv_order = [s for s in indiv_order if s in indivs_gts]
                restart = True
                break

            pct = i * 100 // total_pairs
            bar = "#" * (pct // 2) + "-" * (50 - pct // 2)
            sys.stderr.write("\r[{0}] {1}% ({2}/{3})".format(bar, pct, i, total_pairs))

        if not restart:
            sys.stderr.write("\n")
            break

    if removed_samples:
        sys.stderr.write(
            "Excluded {0} sample(s): {1}\n".format(
                len(removed_samples), ", ".join(removed_samples)
            )
        )

    # Output distance matrix
    out = open(output, "w") if output else sys.stdout

    out.write("\t{0}\n".format(SEPARATOR.join(indiv_order)))
    for indiv1 in indiv_order:
        matrix_row = [indiv1]
        for indiv2 in indiv_order:
            matrix_row.append(get_match_score(pairs_gds, indiv1, indiv2))
        out.write(SEPARATOR.join(matrix_row) + "\n")

    if output:
        out.close()
        sys.stderr.write("Matrix written to {}\n".format(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "vcf_filename", metavar="vcf_file", help="input file with SNP data (`.vcf`)"
    )
    parser.add_argument(
        "pop_filename",
        metavar="pop_file",
        nargs="?",
        default=None,
        help="optional text file (tsv or csv) with individuals "
        "and populations (controls output order; if omitted, "
        "VCF sample order is used)",
    )
    parser.add_argument(
        "--min-threshold",
        type=int,
        default=DEFAULT_MIN_THRESHOLD,
        help="minimum number of shared SNPs without missing data "
        "(default: {})".format(DEFAULT_MIN_THRESHOLD),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="instead of exiting on below-threshold pairs, iteratively "
        "remove the sample with the fewest genotyped loci and retry",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="output file for distance matrix (default: stdout)",
    )
    args = parser.parse_args()
    main(
        args.vcf_filename,
        args.pop_filename,
        min_threshold=args.min_threshold,
        force=args.force,
        output=args.output,
    )
