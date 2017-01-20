#!/usr/bin/env python
"""
Clusters reads of paired-end RAD-seq data for downstream contig
assembly. It maps R1 reads to a reference, and then outputs those reads and
the corresponding R2 reads to a separate 'shuffled' FASTQ file per locus.

Note: when using an existing output folder, reads are being appended to
existing files (use this to append data from multiple samples). BWA needs to
be installed and accessible through `PATH` environmental variable.
"""
import os
import sys
import argparse
import subprocess

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

SAM_HEADER_ID = '@'
FASTQ_HEADER_ID = '@HWI'
COL_QNAME = 0
COL_FLAG = 1
COL_RNAME = 2
COL_POS = 3
COL_MAPQ = 4
COL_SEQ = 9


def is_header(line, header_id):
    """ Return True if line is a header """
    return line[:len(header_id)] == header_id


def execute_to_screen(cmd):
    """ Execute shell command """
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    lines_iterator = iter(popen.stdout.readline, b"")
    for line in lines_iterator:
        print(line)


def execute_to_file(cmd, output_filename):
    """ Execute shell command with stdout redirected to file """
    with open(output_filename, "w") as output_file:
        subprocess.call(cmd, stdout=output_file)


def bwa_index(index_name, ref_fasta_filename):
    """ Index sequences using bwa index """
    # Assess whether BWA index file already exists
    if os.path.isfile('{0}.bwt'.format(index_name)):
        print('(skipped: index file already exists)')
    else:
        cmd = ['bwa', 'index', '-p', index_name, ref_fasta_filename]
        execute_to_screen(cmd)


def bwa_mem(index_name, query_filename, threads, sam_output_filename):
    """ Map reads against reference using bwa mem """
    cmd = ['bwa', 'mem', '-t', threads, index_name, query_filename]
    execute_to_file(cmd, sam_output_filename)


def extract_matches_from_sam_file(sam_filename):
    """ Extract mapped reads (and corresponding ref ids) from SAM file """
    sam_file = open(sam_filename, 'r')
    reference_matches = {}
    for line in sam_file:
        if not is_header(line, SAM_HEADER_ID):
            cols = line.split('\t')
            # Only consider mapped reads (FLAG = 0) and starting at the first
            # position of the reference (POS = 1)
            if int(cols[COL_FLAG]) == 0 and int(cols[COL_POS]) == 1:
                reference_matches[cols[COL_QNAME]] = cols[COL_RNAME]
    sam_file.close()
    return reference_matches


def output_fastq_sequence(output_folder, ref_id, seq_r1, seq_r2):
    """ Store FASTQ sequences in separate files based on ref_id """
    # Append to output file with refid as filename
    output_filename = '{0}/{1}.fastq'.format(output_folder, ref_id)
    output_file = open(output_filename, 'a')
    # Output both R1 and R2 (sequentially)
    # for line in seq_r1:
    #    output_file.write(line)
    for line in seq_r2:
        output_file.write(line)
    output_file.close()


def cluster_fastq_files(reference_matches, r1_fastq_filename,
                        r2_fastq_filename, output_folder):
    """ Cluster FASTQ reads based on mapped R1 reads... """
    # Initialise temporary variables
    ref_id = False
    seq_r1 = []
    seq_r2 = []
    line_count = 0
    total_line_count = 0
    # Iterate through both R1 and R2 files at the same time
    with open(r1_fastq_filename) as r1, open(r2_fastq_filename) as r2:
        for line_r1, line_r2 in zip(r1, r2):
            # Progress update every 100,000 sequences
            if line_count >= 100000:
                print('outputted {0} read pairs...'.format(total_line_count))
                sys.stdout.flush()
                line_count = 0
            line_count += 1
            total_line_count += 1
            #
            if is_header(line_r1, FASTQ_HEADER_ID):
                # Output previous sequence reads if flagged
                if ref_id:
                    output_fastq_sequence(output_folder, ref_id, seq_r1,
                                          seq_r2)
                    # Reset temporary variables
                    del seq_r1[:]
                    del seq_r2[:]
                    seq_r1 = []
                    seq_r2 = []
                    ref_id = False
                # Evaluate current sequence
                name_r1 = line_r1.split()[0][1:]
                name_r2 = line_r2.split()[0][1:]
                # Check if FASTQ files are in sync
                if name_r1 != name_r2:
                    print('FASTQ files of R1 and R2 not in sync')
                    sys.exit('FASTQ files of R1 and R2 not in sync')
                # If mapped read then store sequence info in temporary list
                if name_r1 in reference_matches:
                    ref_id = reference_matches[name_r1]
                    seq_r1.append(line_r1)
                    seq_r2.append(line_r2)
            elif ref_id:
                seq_r1.append(line_r1)
                seq_r2.append(line_r2)


def main(r1_fastq_filename, r2_fastq_filename, ref_fasta_filename, threads,
         output_folder):
    # Index reference sequences in BWA
    print('*** Indexing reference sequences in BWA...')
    index_name = os.path.splitext(ref_fasta_filename)[0]
    bwa_index(index_name, ref_fasta_filename)

    # Map R1 reads against reference (outputted to '[r1]_vs_[ref].sam')
    print('*** Mapping R1 reads against reference sequences with BWA-MEM...')
    sam_output_filename = '{0}_vs_{1}.sam'.format(
        os.path.splitext(r1_fastq_filename)[0], index_name)
    bwa_mem(index_name, r1_fastq_filename, threads, sam_output_filename)

    # Identify mapped reads and corresponding references and store in dict
    print('*** Storing mapped R1 ids and corresponding ref ids into memory...')
    reference_matches = extract_matches_from_sam_file(sam_output_filename)

    # Cluster FASTQ reads based on mapped R1 reads
    print('*** Clustering FASTQ reads based on mapped R1 reads...')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    cluster_fastq_files(reference_matches, r1_fastq_filename,
                        r2_fastq_filename, output_folder)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('r1_fastq_filename', metavar='r1_fastq_file',
                        help='file in FASTQ format with R1 reads')
    parser.add_argument('r2_fastq_filename', metavar='r2_fastq_file',
                        help='file in FASTQ format with R2 reads')
    parser.add_argument('ref_fasta_filename', metavar='ref_fasta_file',
                        help='file in FASTA format with reference contigs')
    parser.add_argument('threads', metavar='threads',
                        help='number of threads to be used by BWA')
    parser.add_argument('output_folder', metavar='output_folder',
                        help='name of output folder')
    args = parser.parse_args()
    main(args.r1_fastq_filename, args.r2_fastq_filename,
         args.ref_fasta_filename, args.threads, args.output_folder)
