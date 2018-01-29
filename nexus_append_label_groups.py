#!/usr/bin/env python
"""
Appends group to each label in a NEXUS tree file.
"""
import argparse

TAXLABELS = "taxlabels"
SEMICOLON = ";"
UNDEFINED = "UNDEFINED"

def get_label_groups(group_filename):
    """ Get list of indvs and groups from tsv or csv """
    label_groups = {}
    group_file = open(group_filename, 'r')
    for line in group_file:
        rows = line.replace(',', '\t').split()
        label_groups[rows[0]] = rows[1]
    return label_groups


def print_new_taxlabels_multi_line(taxlabel_list, label_groups):
    """ Print new tax labels with groups """
    print(TAXLABELS)
    for taxlabel in taxlabel_list:
        taxlabel_only = taxlabel.split('[')[0]
        if taxlabel_only in label_groups:
            print('{0}_{1}'.format(taxlabel_only,
                                            label_groups[taxlabel_only]))
        else:
            print('{0}_{1}'.format(taxlabel_only, UNDEFINED))
    print(SEMICOLON)


def main(nexus_filename, group_filename):
    label_groups = get_label_groups(group_filename)

    taxlabels_flag = False
    taxlabels = []
    taxlabels_line = ''

    nexus_file = open(nexus_filename, 'r')
    for line in nexus_file:
        if line.strip()[0:len(TAXLABELS)].lower() == TAXLABELS:
            taxlabels_flag = True
            if line.rstrip()[len(line.rstrip()) - 1] == SEMICOLON:
                print_new_taxlabels_multi_line(line.split()[1:],
                                               label_groups)
                taxlabels_flag = False
        elif taxlabels_flag:
            if line.strip()[0] == SEMICOLON:
                print_new_taxlabels_multi_line(taxlabels, label_groups)
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
    parser.add_argument('group_filename', metavar='group_filename',
                        help='file with samples and corresponding groups')
    args = parser.parse_args()

    main(args.nexus_filename, args.group_filename)
