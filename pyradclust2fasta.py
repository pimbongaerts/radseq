#!/usr/bin/env python
"""
Creates one large FASTA from all PyRAD clustS files in directory. Only
outputs clusters that exceed size threshold (min. number of sequences in
cluster). First sequence of each cluster is outputted (together with size
of overall cluster - note: not of that specific sequence). Prints the outputted
and total number of clusters to STDOUT.
"""
import sys
import os
import gzip
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'

FASTA_HEADER_CHAR = '>'
COMMENT_CHAR = '/'
CLUSTS_EXT = '.clustS.gz'


def output_clusters_to_fasta(clustS_filename, threshold, output_file):
    """ Iterate over clustS file and output clusters that exceed threshold """
    current_locus = current_sequence_line = ""
    cluster_size = clusters_tot = clusters_out = 0
    clustS_file = gzip.open(clustS_filename, 'rb')
    for gzip_line in clustS_file:
        line = gzip_line.decode("utf-8")
        if line[0] == FASTA_HEADER_CHAR:
            # Extract locus name if first sample in cluster
            if not current_locus:
                current_locus = line[1:].split(';')[0]
            # Extract size from sequence header and add to cluster size
            size = int(line[1:].split(';')[1].split('=')[1])
            cluster_size += size
        elif line[0] not in (COMMENT_CHAR, FASTA_HEADER_CHAR):
            # Store sequence line
            current_sequence_line = line
        elif line[0] == COMMENT_CHAR:
            # Update cluster stats
            clusters_tot += 1
            # Output cluster when above threshold
            if cluster_size >= threshold:
                fasta_header = '>{0};cluster_size={1}\n'.format(current_locus,
                                                                cluster_size)
                output_file.write(fasta_header)
                output_file.write(current_sequence_line)
                clusters_out += 1
            # Reset cluster variables
            current_locus = current_sequence_line = ""
            cluster_size = 0
    clustS_file.close()
    return clusters_out, clusters_tot


def main(path, threshold, output_filename):
    # Open output file and print header
    output_file = open(output_filename, 'w')

    # Iterate over all files
    clusters_out_tot = 0
    for file in os.listdir(path):
        if file.endswith(CLUSTS_EXT):
            clusters_out, clusters_tot = output_clusters_to_fasta(file,
                                                                  threshold,
                                                                  output_file)
            print('{0}\t{1:,}\t{2:,}'.format(file.replace(CLUSTS_EXT, ''),
                                             clusters_out,
                                             clusters_tot))
            clusters_out_tot += clusters_out

    # Close output file
    output_file.close()
    print('Total seqs in {0}: {1:,}'.format(output_filename, clusters_out_tot))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', metavar='path',
                        help='path that contains PyRAD `.clustS` files')
    parser.add_argument('threshold', type=int, metavar='cluster_threshold',
                        help='minimum size of cluster to be included')
    parser.add_argument('output_filename', metavar='output_file',
                        help='name of output FASTA file')
    args = parser.parse_args()
    main(args.path, args.threshold, args.output_filename)
