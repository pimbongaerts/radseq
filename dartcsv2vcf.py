#!/usr/bin/env python
"""
dartcsv2vcf.py: Converts DArT-seq CSV file (Report_*_SNP_mapping_1.csv) to VCF
"""
import sys
import math
import argparse

CSV_COMMENT = "*"
CSV_HEADER_START = "AlleleID"

KNOWN_METADATA_COLS = {
    "AlleleID",
    "CloneID",
    "AlleleSequence",
    "AlleleSequenceRef",
    "AlleleSequenceSnp",
    "TrimmedSequence",
    "TrimmedSequenceRef",
    "TrimmedSequenceSnp",
    "SNP",
    "SnpPosition",
    "CallRate",
    "OneRatioRef",
    "OneRatioSnp",
    "FreqHomRef",
    "FreqHomSnp",
    "FreqHets",
    "PICRef",
    "PICSnp",
    "AvgPIC",
    "AvgCountRef",
    "AvgCountSnp",
    "RepAvg",
}

VCF_META_HEADER = """\
##fileformat=VCFv4.2
##source=dartcsv2vcf.py
##INFO=<ID=DP,Number=1,Type=Float,Description="Average total read depth (AvgCountRef + AvgCountSnp)">
##INFO=<ID=CR,Number=1,Type=Float,Description="Call rate (fraction of samples with a genotype call)">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency (FreqHomSnp + 0.5 * FreqHets)">
##INFO=<ID=HET,Number=1,Type=Float,Description="Observed heterozygote frequency">
##INFO=<ID=PIC,Number=1,Type=Float,Description="Average Polymorphism Information Content">
##INFO=<ID=REP,Number=1,Type=Float,Description="Reproducibility average">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
VCF_COL_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2021 Pim Bongaerts"
__license__ = "GPL"


def get_output_filename(csv_filename, extension):
    return "{0}.{1}".format(csv_filename.replace(".csv", ""), extension)


def get_dict_of_name_changes(samplenames_filename):
    """Initialise dictionary with old (keys) and new (values)"""
    rename_file = open(samplenames_filename, "r")
    name_changes = {}
    for line in rename_file:
        cols = line.replace(",", "\t").split()
        name_changes[cols[0].strip()] = cols[1].strip()
    rename_file.close()
    return name_changes


def get_new_sample_names(old_samples, name_changes):
    """Change samples to new sample names"""
    sample_names = []
    for sample in old_samples:
        if sample.strip() in name_changes:
            sample_names.append(name_changes[sample.strip()])
        else:
            print("Not renamed: {0}".format(sample))
            sample_names.append(sample.strip())
    return sample_names


def find_column_indices(header_cols):
    """Detect column positions from header names, returns dict"""
    stripped = [c.strip() for c in header_cols]
    col_map = {name: i for i, name in enumerate(stripped)}

    idx = {}
    idx["CloneID"] = col_map.get("CloneID")
    idx["SNP"] = col_map.get("SNP")
    idx["SnpPosition"] = col_map.get("SnpPosition")
    idx["AlleleSeq"] = col_map.get("AlleleSequenceRef", col_map.get("AlleleSequence"))

    missing = [k for k, v in idx.items() if v is None]
    if missing:
        sys.exit("Error: missing required columns: {}".format(", ".join(missing)))

    for name in (
        "CallRate",
        "AvgCountRef",
        "AvgCountSnp",
        "RepAvg",
        "FreqHomSnp",
        "FreqHets",
        "AvgPIC",
    ):
        idx[name] = col_map.get(name)

    first_sample_col = None
    for i, name in enumerate(stripped):
        if name not in KNOWN_METADATA_COLS:
            first_sample_col = i
            break
    if first_sample_col is None:
        sys.exit("Error: no sample columns found in header")
    idx["first_sample"] = first_sample_col

    return idx


def convert_to_vcf_genotype(genotype):
    """Convert DaRT to VCF genotype"""
    dart_gt = genotype.strip()
    if dart_gt == "-":
        vcf_gt = "./."
    elif dart_gt == "0":
        vcf_gt = "0/0"
    elif dart_gt == "1":
        vcf_gt = "1/1"
    elif dart_gt == "2":
        vcf_gt = "0/1"
    else:
        print(dart_gt)
        vcf_gt = "err"
    return vcf_gt


def get_vcf_genotypes(dart_genotypes):
    """Get VCF genotypes"""
    return "\t".join([convert_to_vcf_genotype(dart_gt) for dart_gt in dart_genotypes])


def repavg_to_phred(repavg):
    """Convert reproducibility (0-1) to Phred-scaled quality score"""
    try:
        r = float(repavg)
    except (ValueError, TypeError):
        return "."
    if r >= 1.0:
        return "100"
    if r <= 0.0:
        return "0"
    return "{:.1f}".format(-10 * math.log10(1 - r))


def build_info_field(cols, idx):
    """Build VCF INFO field from DArT metadata columns"""
    info_parts = []

    if idx.get("AvgCountRef") is not None and idx.get("AvgCountSnp") is not None:
        try:
            dp = float(cols[idx["AvgCountRef"]]) + float(cols[idx["AvgCountSnp"]])
            info_parts.append("DP={:.1f}".format(dp))
        except (ValueError, IndexError):
            pass

    if idx.get("CallRate") is not None:
        try:
            info_parts.append("CR={}".format(cols[idx["CallRate"]].strip()))
        except IndexError:
            pass

    if idx.get("FreqHomSnp") is not None and idx.get("FreqHets") is not None:
        try:
            af = float(cols[idx["FreqHomSnp"]]) + 0.5 * float(cols[idx["FreqHets"]])
            info_parts.append("AF={:.4f}".format(af))
        except (ValueError, IndexError):
            pass

    if idx.get("FreqHets") is not None:
        try:
            info_parts.append("HET={}".format(cols[idx["FreqHets"]].strip()))
        except IndexError:
            pass

    if idx.get("AvgPIC") is not None:
        try:
            info_parts.append("PIC={}".format(cols[idx["AvgPIC"]].strip()))
        except IndexError:
            pass

    if idx.get("RepAvg") is not None:
        try:
            info_parts.append("REP={}".format(cols[idx["RepAvg"]].strip()))
        except IndexError:
            pass

    return ";".join(info_parts) if info_parts else "."


def output_vcf_line(vcf_file, chrom, pos, ref, alt, qual, info, vcf_genotypes_concat):
    """Output VCF line"""
    vcf_file.write(
        "{0}\t{1}\t.\t{2}\t{3}\t{4}\tPASS\t{5}\tGT\t{6}\n".format(
            chrom, pos, ref, alt, qual, info, vcf_genotypes_concat
        )
    )


def main(dartcsv_filename, samplenames_filename=None):
    name_changes = None
    if samplenames_filename:
        name_changes = get_dict_of_name_changes(samplenames_filename)

    vcf_file = open(get_output_filename(dartcsv_filename, "vcf"), "w")
    fasta_file = open(get_output_filename(dartcsv_filename, "fa"), "w")

    vcf_file.write(VCF_META_HEADER)

    idx = None

    csv_file = open(dartcsv_filename, "r")
    for line in csv_file:
        cols = line.split(",")
        if line.startswith(CSV_HEADER_START):
            idx = find_column_indices(cols)

            if name_changes:
                sample_names = get_new_sample_names(
                    cols[idx["first_sample"] :], name_changes
                )
            else:
                sample_names = [s.strip() for s in cols[idx["first_sample"] :]]

            vcf_file.write("{0}\t{1}\n".format(VCF_COL_HEADER, "\t".join(sample_names)))
        elif line[0] != CSV_COMMENT and idx is not None:
            ref_allele, divider, alt_allele = cols[idx["SNP"]].split(":")[1]
            vcf_genotypes_concat = get_vcf_genotypes(cols[idx["first_sample"] :])

            qual = "."
            if idx.get("RepAvg") is not None:
                qual = repavg_to_phred(cols[idx["RepAvg"]])
            info = build_info_field(cols, idx)

            output_vcf_line(
                vcf_file,
                cols[idx["CloneID"]],
                cols[idx["SnpPosition"]],
                ref_allele,
                alt_allele,
                qual,
                info,
                vcf_genotypes_concat,
            )
            fasta_file.write(">{0}\n".format(cols[idx["CloneID"]]))
            fasta_file.write("{0}\n".format(cols[idx["AlleleSeq"]]))

    csv_file.close()
    vcf_file.close()
    fasta_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "dartcsv_filename", metavar="dartcsv_filename", help="DArT-seq csv filename"
    )
    parser.add_argument(
        "samplenames_filename",
        metavar="samplenames_file",
        nargs="?",
        default=None,
        help="optional text file (tsv or csv) with old and "
        "new name for each sample (not all samples have to "
        "be listed)",
    )
    args = parser.parse_args()
    main(args.dartcsv_filename, args.samplenames_filename)
