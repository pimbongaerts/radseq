#!/usr/bin/env python
"""
Set the color of each label in a NEXUS tree file.
"""
import argparse

TAXLABELS = "taxlabels"
SEMICOLON = ";"


def get_label_colors(color_filename):
    """ Get list of indvs and colors from tsv or csv """
    label_colors = {}
    color_file = open(color_filename, 'r')
    for line in color_file:
        rows = line.replace(',', '\t').split()
        label_colors[rows[0]] = rows[1]
    return label_colors


def print_new_taxlabels_multi_line(taxlabel_list, label_colors):
    """ Print new tax labels with colors """
    print(TAXLABELS)
    for taxlabel in taxlabel_list:
        taxlabel_only = taxlabel.split('[')[0]
        if taxlabel_only in label_colors:
            print('{0}[&!color={1}]'.format(taxlabel_only,
                                            label_colors[taxlabel_only]))
        else:
            print('{0}[&!color=#000000]'.format(taxlabel_only))
    print(SEMICOLON)


def main(nexus_filename, color_filename):
    label_colors = get_label_colors(color_filename)

    taxlabels_flag = False
    taxlabels = []
    taxlabels_line = ''

    nexus_file = open(nexus_filename, 'r')
    for line in nexus_file:
        if line.strip()[0:len(TAXLABELS)].lower() == TAXLABELS:
            taxlabels_flag = True
            if line.rstrip()[len(line.rstrip()) - 1] == SEMICOLON:
                print_new_taxlabels_multi_line(line.split()[1:],
                                               label_colors)
                taxlabels_flag = False
        elif taxlabels_flag:
            if line.strip()[0] == SEMICOLON:
                print_new_taxlabels_multi_line(taxlabels, label_colors)
                taxlabels_flag = False
            else:
                taxlabels.append(line.strip())
        else:
            print(line, end='')


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nexus_filename', metavar='nexus_filename',
                        help='nexus input file)')
    parser.add_argument('color_filename', metavar='color_filename',
                        help='file with samples and corresponding colors')
    args = parser.parse_args()

    main(args.nexus_filename, args.color_filename)
