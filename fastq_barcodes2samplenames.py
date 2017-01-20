#!/usr/bin/env python
"""
Renames barcoded .fastq.gz files in a folder to sample/individual names.
Takes relative or absolute path as first argument (e.g. 'samples'; without
trailing slash) and a text file as second argument. The latter should be a
tab-separated text file, with the barcode in the first column, and the new
sample/individual name in the second column. The script conducts a trial run
first, listing the files to be renamed, and then asks for confirmation before
doing the actual renaming.
"""
import sys
import os
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


def human_readable_size(number):
    """ Convert byte size to human readable format """
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(number) < 1024.0:
            return '%3.1f%s' % (number, unit)
        number /= 1024.0
    return '%.1f%s' % (number, 'Yi')


def dictionary_from_barcodes(barcode_filename):
    """ Create dictionary with barcodes as key and individual as value """
    new_filenames = {}
    barcode_file = open(barcode_filename, 'r')
    for line in barcode_file:
        line = line.replace(',', ' ').rstrip()
        cols = line.split('\t')
        new_filenames[cols[0]] = cols[1]
    return new_filenames


def rename_files(path, new_filenames, perform_rename=False):
    """ Rename files according to barcode dictionary """
    # Create renaming log file
    rename_log = open('rename_log.txt', 'w')
    # Loop through all files in specified path
    file_count = 0
    for filename in os.listdir(path):
        # Loop through barcodes to see if filename contains barcode
        for barcode in new_filenames.keys():
            if barcode in filename:
                new_filename = '{0}.fastq.gz'.format(new_filenames[barcode])
                new_filename_with_path = '{0}/{1}'.format(path, new_filename)
                filename_with_path = '{0}/{1}'.format(path, filename)
                file_size = os.path.getsize(filename_with_path)
                human_file_size = human_readable_size(file_size)
                try:
                    # Only rename file if perform_rename = True
                    if perform_rename:
                        os.rename(filename_with_path, new_filename_with_path)
                    msg = '{0}\t{1}\t{2}\t{3}\t{4}'.format(filename, barcode,
                                                           new_filename,
                                                           file_size,
                                                           human_file_size)
                    file_count += 1
                except:
                    msg = 'Error: could not rename {0}'.format(filename)
                print(msg)
                rename_log.write('{0}\n'.format(msg))
    rename_log.close()
    return file_count


def ask_for_permission(question):
    """ Ask for user confirmation using input """
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    while True:
        sys.stdout.write('{0} [y/n] '.format(question))
        choice = input().lower()
        if choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with yes or no (or y or n). \n ')


def main(path, barcode_filename):
    # Initialise dictionary with barcodes
    new_filenames = dictionary_from_barcodes(barcode_filename)
    # Perform a test rename to show user (perform_rename = False)
    file_count = rename_files(path, new_filenames, False)
    # Summarise output of test rename
    barcode_count = len(new_filenames)
    print(('\nFound {0} out of {1} barcodes'
           ' listed in {2}').format(file_count, barcode_count,
                                    barcode_filename))
    # Ask for permission to conduct actual rename
    if ask_for_permission('Based on the above - would you like to continue '
                          'with conducting the actual renaming?'):
        file_count = rename_files(path, new_filenames, True)
        # Summarise output of actual rename
        print(('\nRenamed {0} files - {1} barcodes'
               ' listed').format(file_count, barcode_count))
    else:
        print('\nNo files were renamed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', metavar='path',
                        help='path (for current directory use `.`)')
    parser.add_argument('barcode_filename', metavar='barcode_file',
                        help='text file (tsv or csv) with barcodes and sample \
                        names')
    args = parser.parse_args()
    main(args.path, args.output_filename)
