#!/usr/bin/env python
"""
Compiles README markdown file for this repository
(https://github.com/pimbongaerts/radseq).
Categories are assigned based on prefix, usage information is extracted
from argparse, and example input files are assigned based on argument names.
"""
import sys
import os
import argparse
import glob
import subprocess

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_FILE = 'README_header.md'
PATH = '.'
DOCSTRING = '"""'
GROUP_TAGS = ['vcf', 'pyrad', 'fastq', 'mapping', 'popfile', 'other']
GROUP_OTHER = 'other'
EXAMPLES_FOLDER = 'input_examples/'
EXAMPLES_DESCR = 'Example input file(s): '
PEP8_FAIL = ' *[File did not pass PEP8 check]*'


def get_python_help_message(filename):
    """ Capture help mesage of python script """
    help_msg = subprocess.Popen(['python3', filename, '-h'],
                                stdout=subprocess.PIPE).communicate()[0]
    return help_msg


def check_PEP8_compliance(filename):
    """ Capture help mesage of python script """
    pep8_output = subprocess.Popen(['pycodestyle', filename, ],
                                   stdout=subprocess.PIPE).communicate()[0]
    if pep8_output.decode("utf-8") == '':
        return ''
    else:
        return PEP8_FAIL


def extract_metadata_from_script(filename):
    """ Extract metadata from script """
    metadata = {}

    # Extract information from argparse help output
    help_msg = get_python_help_message(filename)
    help_metadata = help_msg.decode("utf-8").split('\n\n')
    metadata['usage'] = help_metadata[0]
    metadata['title'] = help_metadata[1]
    metadata['arguments'] = help_metadata[2:]
    # Extract group from filename
    for group_tag in GROUP_TAGS:
        if group_tag in filename:
            metadata['group'] = group_tag
            break
        metadata['group'] = GROUP_OTHER
    # Evaluate if file passes PEP8 check
    metadata['pep8'] = check_PEP8_compliance(filename)

    return metadata


def get_example_hyperlinks(usage):
    """ Return hyperlinks to example input files if they exist """
    hyper_links = []
    usage_split = usage.replace('\n', '').split(']')
    potential_files = usage_split[len(usage_split) - 1].split()

    for potential_file in potential_files:
        potential_file_path = '{0}{1}'.format(EXAMPLES_FOLDER,
                                              potential_file.strip())
        actual_file = glob.glob('{0}.*'.format(potential_file_path))
        if actual_file:
            actual_file_short = actual_file[0].replace(EXAMPLES_FOLDER, '')
            hyper_links.append('[{0}]({1})'.format(actual_file_short,
                                                   actual_file[0]))
    if len(hyper_links) > 0:
        return '{0} {1}.'.format(EXAMPLES_DESCR, ', '.join(hyper_links))
    else:
        return ''


def get_formatted_metadata(script, script_metadata):
    """ Return formatted metadata for script """
    metadata = ''
    # Title and description
    metadata += '**[{0}]({0})** - {1}'.format(script,
                                              script_metadata[script]['title'])
    # Check for PEP8 compliance
    metadata += script_metadata[script]['pep8']
    # Print code block with usage and argument info
    metadata += '\n\n\t{0}\n'.format(script_metadata[script]['usage'])
    arguments = '\n\n'.join(script_metadata[script]['arguments'])
    metadata += '\n\t{0}'.format(arguments.replace('\n', '\n\t'))
    metadata += '\n\n'
    # Print input file example hyperlinks
    metadata += get_example_hyperlinks(script_metadata[script]['usage'])
    metadata += '\n\n'
    return metadata


def main():
    # Output markdown header (from separate file)
    header_file = open(HEADER_FILE, 'r')
    print(header_file.read())
    header_file.close()

    # Loop over all python scripts in path and extract metadata
    grouped_scripts = {}
    script_metadata = {}
    for filename in os.listdir(PATH):
        if filename.endswith('.py'):
            metadata = extract_metadata_from_script(filename)
            # Add script to dict with grouping
            if not metadata['group'] in grouped_scripts:
                grouped_scripts[metadata['group']] = []
            grouped_scripts[metadata['group']].append(filename)
            # Store metadata in dict
            script_metadata[filename] = metadata

    # Output script summary by group:
    # print('\n|script|description (truncated)|\n|---|---|')
    # for group in grouped_scripts:
    #    # Output name of each group
    #    print('|**{0}**| |'.format(group))
    #    # Output hyperlinks to each script in group
    #    for script in grouped_scripts[group]:
    #        print('{0}|'.format(script), end='')
    #        print('{0}...|'.format(script_metadata[script]['title'][0:40]))

    # Output script metadata by group:
    for group in GROUP_TAGS:
        # Output header for group
        print('\n## {0}\n'.format(group))

        # Output formatted metadata for each script
        for script in grouped_scripts[group]:
            print(get_formatted_metadata(script, script_metadata))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    main()
