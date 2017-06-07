#!/usr/bin/env python
"""
Plots all the results from a structure_mp run to a multi-page PDF.
"""

import argparse
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2017 Pim Bongaerts"
__license__ = "GPL"

CLUMPP_PREFIX = 'clumpp_K'
CLUMPP_NAME = '.out.csv'
COLORS = ['#9e0142', '#5e4fa2', '#d53e4f', '#3288bd', '#f46d43', '#66c2a5',
          '#fdae61', '#abdda4', '#fee08b', '#e6f598', '#ffffbf']


def get_csv_files_in_path(path):
    """ Return list with all csv files in path """
    return [filename for filename in glob.glob(os.path.join(path, '*.csv'))]


def group_csv_files_by_K(csv_files):
    """ Group csv files by the number of clusters """
    csv_files_by_K = {}
    for filename in csv_files:
        K_value = filename.split(CLUMPP_PREFIX)[1].split('_')[0].split('.')[0]
        if K_value not in csv_files_by_K:
            csv_files_by_K[K_value] = []
        csv_files_by_K[K_value].append(filename)
    return csv_files_by_K


def plot_from_csv_file(csv_file, popnames):
    """ Generate plot from data in csv file """
    df = pd.read_csv(csv_file, header=None)

    # Plot bar settings
    bar_width = 1
    bar_left_pos = range(0, df.shape[0])
    tick_pos = [i + (bar_width / 2) for i in bar_left_pos]
    bottom_pos = np.zeros_like(bar_left_pos).astype('float')
    if popnames:
        names = df[df.columns[1]]
    else:
        names = df[df.columns[0]]

    # Generate plot and1 bars
    f, ax = plt.subplots(1, figsize=(10, 5))
    for i in range(2, len(df.columns)):
        ax.bar(bar_left_pos, df[df.columns[i]], width=bar_width,
               bottom=bottom_pos, color=COLORS[i - 2])
        bottom_pos += df[df.columns[i]]

    # Plot properties
    plt.title('STRUCTURE {0}'.format(os.path.basename(csv_file)))
    plt.xticks(tick_pos, names)
    plt.xlim([min(tick_pos) - bar_width, max(tick_pos) + bar_width])
    plt.ylim(0, 1)
    plt.setp(plt.gca().get_xticklabels(), rotation=90,
             horizontalalignment='right', size='xx-small')
    ax.set_xlabel('Individuals')
    ax.set_ylabel('Assignment')
    f.subplots_adjust(right=0.9, bottom=0.2)
    return plt


def output_plots_to_pdf(pdf_output_filename, csv_files_by_K, clumpp_only,
                        popnames):
    """ Output plots (sorted by K) to multi-page PDF """
    with PdfPages(pdf_output_filename) as pdf:
        for K_value in sorted(csv_files_by_K):
            print('Outputting plots for K = {0}'.format(K_value))
            for csv_file in csv_files_by_K[K_value]:
                # If clumpp flag set  then only output summary
                if (not clumpp_only) or (CLUMPP_NAME in csv_file):
                    plot = plot_from_csv_file(csv_file, popnames)
                    pdf.savefig()
                    plot.close()


def main(path, clumpp_only, popnames):
    csv_files = get_csv_files_in_path(path)
    csv_files_by_K = group_csv_files_by_K(csv_files)
    pdf_output_filename = '{0}.pdf'.format(path.replace('/', ''))
    output_plots_to_pdf(pdf_output_filename, csv_files_by_K, clumpp_only,
                        popnames)
    os.system('open {0}'.format(pdf_output_filename))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', metavar='path',
                        help='path to structure_mp results')
    parser.add_argument('-p', '--popnames', action='store_true',
                        help='set flag to output population names')
    parser.add_argument('-c', '--clumpp_only', action='store_true',
                        help='set flag to only plot CLUMPP summary')
    args = parser.parse_args()
    main(args.path, args.clumpp_only, args.popnames)
