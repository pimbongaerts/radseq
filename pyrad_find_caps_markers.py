#!/usr/bin/env python
"""
Search PyRAD output file for diagnostic CAPS loci that can distinguish two
groups (or one group and all other samples).
"""
import sys
import os
import re
import errno
import itertools
import argparse

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2017 Pim Bongaerts'
__license__ = 'GPL'

CONSIDERED_CHARS = ['A', 'C', 'G', 'T']
UNASSIGNED = 'unassigned'
OTHER_SAMPLES = '#others'
OUTPUT_HIGHLIGHT_START = '\x1b[6;30;42m'
OUTPUT_HIGHLIGHT_END = '\x1b[0m'


class SampleGroups(object):
    """ Groups of samples that are to be compared """

    def __init__(self, group1_name, group2_name):
        self.group1_name = group1_name
        if group2_name:
            self.group2_name = group2_name
        else:
            # If not specified, Group2 contains indivs from all other groups
            self.group2_name = OTHER_SAMPLES

        self.group1_samples = []   # Indivs in Group1
        self.group2_samples = []   # Indivs in Group2
        self.other_groups = {}     # Other groups (keys) & their indivs (lists)
        self.sample_groups = {}    # All indivs and their group assignment

    def add_sample(self, sample_name, group_name):
        """ Add sample and its group assignment """
        self.sample_groups[sample_name] = group_name

        if group_name == self.group1_name:      # Group 1
            self.group1_samples.append(sample_name)
        elif group_name == self.group2_name:    # Group 2
            self.group2_samples.append(sample_name)
        else:
            # Other groups (create if not yet present)
            if group_name not in self.other_groups.keys():
                self.other_groups[group_name] = []
            self.other_groups[group_name].append(sample_name)

    def get_sample_group(self, sample_name):
        """ Return group to which sample belongs """
        if sample_name in self.sample_groups:
            return self.sample_groups[sample_name]
        else:
            return UNASSIGNED      # Indiv in pyrad_file but not group_file


class SeqGroups(object):
    """ Sequences of groups of samples that are to be compared """

    def __init__(self, group1_name, group2_name):
        self.locus_name = ''
        self.group1_name = group1_name
        if group2_name:
            self.group2_name = group2_name
        else:
            # If not specified, Group2 contains indivs from all other groups
            self.group2_name = OTHER_SAMPLES

        self.group1_seqs = {}       # Seqs of indivs in Group1
        self.group2_seqs = {}       # Seqs of indivs in Group2
        self.other_groups_seqs = {}  # Other groups (keys) & their seqs (dicts)
        self.indiv_seq = {}         # All indivs and their seq

    def set_locus_name(self, locus_name):
        """ Assign a name to current locus """
        self.locus_name = locus_name

    def add_sample(self, sample_name, group_name, seq):
        """ Add sample and sequence to its corresponding group """
        self.indiv_seq[sample_name] = seq

        if group_name == UNASSIGNED:
            return      # Skip indiv that is in pyrad_file but not group_file
        elif group_name == self.group1_name:
            self.group1_seqs[sample_name] = seq
        elif group_name == self.group2_name:
            self.group2_seqs[sample_name] = seq
        else:
            # Other groups (create if not yet present)
            if group_name not in self.other_groups_seqs.keys():
                self.other_groups_seqs[group_name] = {}
            self.other_groups_seqs[group_name][sample_name] = seq
            # Also add to Group 2 (so Group 1 is compared to everything else)
            if self.group2_name == OTHER_SAMPLES:
                self.group2_seqs[sample_name] = seq

    def get_seq(self, indiv):
        """ Retrieve sequence for individual """
        if indiv in self.indiv_seq:
            return self.indiv_seq[indiv]
        else:
            return False

    def reinit(self):
        """ Clear all sequence dictionaries """
        self.group1_seqs.clear()
        self.group2_seqs.clear()
        self.other_groups_seqs.clear()
        self.indiv_seq.clear()

    def sufficient_representation(self, min_samples):
        """ Assess if Group 1 and 2 both have at least min_samples """
        return (len(self.group1_seqs) >= min_samples and
                len(self.group2_seqs) >= min_samples)

    def get_potential_caps_markers(self, re_sites):
        """ Assess potential of each re_site as CAPS marker to distinguish
            sequence groups 1 and 2 """

        # Create consensus sequences of both groups to evaluate for CAPS
        group1_cons_seq = self.__create_consensus_seq(self.group1_seqs)
        group2_cons_seq = self.__create_consensus_seq(self.group2_seqs)

        # Evaluate each individual restriction site as a CAPS marker to
        # distinguish both consensus seqs
        caps_results = []
        for re_name, re_seq in re_sites.items():
            caps_results += self.__evaluate_caps_site(re_name, re_seq,
                                                      group1_cons_seq,
                                                      group2_cons_seq)
        return caps_results

    def __evaluate_caps_site(self, re_name, re_seq, seq1, seq2):
        """ Evaluate specific re_site as CAPS marker to distinguish sequence
            groups 1 and 2 """

        # Find occurrences of restriction site in both sequences
        rs_pos_in_seq1 = [x.start() for x in re.finditer(re_seq, seq1)]
        rs_pos_in_seq2 = [x.start() for x in re.finditer(re_seq, seq2)]
        rs_pos_combined = set(rs_pos_in_seq1 + rs_pos_in_seq2)

        # Evaluate potential caps sites (present in one, mutated in other)
        caps_pos_in_seq = []
        for rs_pos in rs_pos_combined:
            match_seq1 = seq1[rs_pos:rs_pos + len(re_seq)]  # Seq for Group 1
            match_seq2 = seq2[rs_pos:rs_pos + len(re_seq)]  # Seq for Group 2
            match_combined = match_seq1 + match_seq2       # All nucleotides

            # Identify as CAPS when (1) re_seq present in only one group, and
            #                       (2) there are no ambiguous sites ('N')
            if match_seq1 != match_seq2 and 'N' not in (match_combined):
                caps_pos_in_seq.append(rs_pos)

        if len(caps_pos_in_seq) > 0:
            return self.__format_caps_result(re_name, re_seq, seq1,
                                             seq2, rs_pos_combined)
        else:
            return []  # TODO: or return False?

    def __format_caps_result(self, re_name, re_seq, group1_cons_seq,
                             group2_cons_seq, rs_pos_combined):
        """ Format successful CAPS match lines for printing to stdout """
        caps_result = []
        group1_seq_count = len(self.group1_seq)
        group2_seq_count = len(self.group2_seqs)

        caps_result.append(self.__format_caps_header(self.locus_name,
                                                     re_name, re_seq))
        caps_result.append(self.__format_caps_match(self.group1_name,
                                                    group1_seq_count,
                                                    group1_cons_seq,
                                                    rs_pos_combined,
                                                    re_seq))
        caps_result.append(self.__format_caps_match(self.group2_name,
                                                    group2_seq_count,
                                                    group2_cons_seq,
                                                    rs_pos_combined,
                                                    re_seq))

        for group in self.other_groups_seqs.keys():
            group_cons_seq = self.__create_consensus_seq(
                self.other_groups_seqs[group])
            group_info = self.__format_caps_match(group,
                                                  len(self.other_groups_seqs[
                                                      group]),
                                                  group_cons_seq,
                                                  rs_pos_combined,
                                                  re_seq)
            caps_result.append(group_info)

        return caps_result

    @classmethod
    def __format_caps_match(cls, group_name, group_size, consensus_seq,
                            caps_site_pos, re_seq):
        """ Format sequence match for __format_caps_result  """
        seq_formatted = cls.__highlight_matches(consensus_seq, caps_site_pos,
                                                len(re_seq))
        return '{0} ({1}):\t\t{2}'.format(group_name, group_size,
                                          seq_formatted)

    @classmethod
    def __create_consensus_seq(cls, group_seqs):
        """ Create the consensus sequence for a group of sequences -
            note: any variable sites will be outputted as 'N' """
        pos_bases = {}
        seq_length = 0

        # Create set for each sequence position with different base options
        for sample, seq in group_seqs.items():
            seq_length = len(seq)
            for pos in range(0, seq_length - 1):
                if pos not in pos_bases:
                    pos_bases[pos] = set()
                pos_bases[pos].update(seq[pos])

        # Create consensus sequence
        consensus_seq = []
        for pos in range(0, seq_length - 1):
            consensus_seq.append(cls.__get_consensus_base(pos_bases[pos]))
        return ''.join(consensus_seq)

    @staticmethod
    def __format_caps_header(locus_name, re_name, re_seq):
        """ Format sequence header for __format_caps_result  """
        return 'locus_{0}; {1} ({2})'.format(locus_name, re_name, re_seq)

    @staticmethod
    def __get_consensus_base(bases):
        """ Determine consensus base from a set of bases -
            note: all variable sites are marked as ambiguous 'N'  """
        if len(bases) == 1:
            return bases.pop()
        else:
            return 'N'

    @staticmethod
    def __highlight_matches(seq, pos_in_seq, match_length):
        """ Color-code matching positions in sequence for stdout """
        end_pos_in_seq = [x + match_length for x in pos_in_seq]
        formatted_seq = []
        for pos in range(0, len(seq) - 1):
            if pos in pos_in_seq:
                formatted_seq.append(OUTPUT_HIGHLIGHT_START)
            elif pos in end_pos_in_seq:
                formatted_seq.append(OUTPUT_HIGHLIGHT_END)
            formatted_seq.append(seq[pos])
        return ''.join(formatted_seq)


def get_sample_groups(group_filename, group1, group2):
    """ Initialize sample groups from group_file  """
    group_file = open(group_filename, 'r')
    groups = SampleGroups(group1, group2)
    for line in group_file:
        cols = line.replace(',', ' ').split()
        sample_name = cols[0]
        group_name = cols[1]
        groups.add_sample(sample_name, group_name)
    group_file.close()
    return groups


def get_re_sites(re_site_filename):
    """ Initialize sample groups from group_file """
    re_site_file = open(re_site_filename, 'r')
    re_sites = {}
    for line in re_site_file:
        cols = line.replace(',', ' ').split()
        re_site_name = cols[0]
        re_site_seq = cols[1]
        re_sites[re_site_name] = re_site_seq
    re_site_file.close()
    return re_sites


def get_sample_details(line, groups):
    """ Extract sample details from line in pyrad file """
    cols = line.split()
    sample_name = cols[0][1:]
    seq = cols[1].strip()
    group = groups.get_sample_group(sample_name)
    return sample_name, group, seq


def get_locus_name(line):
    """ Extract locus name from line in pyrad file """
    return line.split('|')[1].strip()


def output_to_screen(caps_results):
    """ Output CAPS results to screen """
    for line in caps_results:
        print(line)
    print()


def output_as_fasta(output_folder, header_info, groups, seqgroups):
    """ Output individual sequences to output file """
    locus = header_info.split(';')[0]
    output_file = open('{0}/{1}.fa'.format(output_folder, locus), 'w')
    for indiv in groups.group1_samples:
        seq = seqgroups.get_seq(indiv)
        if seq:
            output_file.write(format_seq(groups.group1_name, indiv,
                                         header_info, seq))
    if groups.group2_name != OTHER_SAMPLES:
        for indiv in groups.group2_samples:
            seq = seqgroups.get_seq(indiv)
            if seq:
                output_file.write(format_seq(groups.group2_name, indiv,
                                             header_info, seq))
    for group in groups.other_groups:
        for indiv in groups.other_groups[group]:
            seq = seqgroups.get_seq(indiv)
            if seq:
                output_file.write(format_seq(group, indiv,
                                             header_info, seq))
    output_file.close()


def format_seq(group, indiv, header, seq):
    """ Format sequence and header into FASTA format """
    if seq:
        return '>{0}_{1}; {2}\n{3}\n'.format(group, indiv, header, seq)


def make_sure_path_exists(path):
    """ Ensure folder exists and raise exception if cannot be created """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def evaluate_loci_in_pyrad_file(pyrad_filename, groups, re_sites,
                                min_samples, output_folder):
    """ Evaluate all loci in PyRAD file for potential CAPS markers """
    seqgroups = SeqGroups(groups.group1_name, groups.group2_name)

    pyrad_file = open(pyrad_filename, 'r')
    for line in pyrad_file:
        # Line with sequence information -> store
        if line[0] == '>':
            seqgroups.add_sample(*get_sample_details(line, groups))

        # Line with locus name (end of sequence alignment) -> evaluate
        if line[0] == '/':
            if seqgroups.sufficient_representation(min_samples):

                # Evaluate for potential CAPS markers
                seqgroups.set_locus_name(get_locus_name(line))
                caps_results = seqgroups.get_potential_caps_markers(re_sites)
                # Output matches to screen (and optional output files)
                if len(caps_results) > 0:
                    output_to_screen(caps_results)
                    if output_folder:
                        output_as_fasta(output_folder,
                                        caps_results[0],
                                        groups, seqgroups)

            # Empty sequences / information in seqgroups object
            seqgroups.reinit()
    pyrad_file.close()


def main(pyrad_filename, group_filename, re_site_filename,
         group1, group2, min_samples, output_folder):

    if output_folder:
        make_sure_path_exists(output_folder)

    groups = get_sample_groups(group_filename, group1, group2)

    re_sites = get_re_sites(re_site_filename)

    evaluate_loci_in_pyrad_file(pyrad_filename,
                                groups,
                                re_sites,
                                min_samples,
                                output_folder)

if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', dest='pyrad_filename',
                        metavar='pyrad_file',
                        help='input pyrad .loci or .alleles file')
    parser.add_argument('-g', '--groups', dest='group_filename',
                        metavar='group_file',
                        help='text file (tsv or csv) separating individuals \
                              (first column) into groups (second column)')
    parser.add_argument('-r', '--re', dest='re_site_filename',
                        metavar='re_site_file',
                        help='text file (tsv or csv) listing restriction \
                              site names (first column) and their recognition \
                              sequences (second column)')
    parser.add_argument('-g1', '--group1', dest='group1',
                        metavar='group1',
                        help='first group that is targeted in marker search')
    parser.add_argument('-g2', '--group2', dest='group2',
                        metavar='group2',
                        help='second group (optional) that is targeted in \
                        marker search (if none given, group1 is contrasted \
                        against all other samples in group_file')
    parser.add_argument('-m', '--min_samples', dest='min_samples',
                        metavar='min_samples', default=5,
                        help='minimum number of genotyped samples in each \
                              group  for a marker to be considered')
    parser.add_argument('-o', '--output', dest='output_folder',
                        metavar='output_folder',
                        help='name of output folder for individual seqs of \
                              each diagnostic locus')
    args = parser.parse_args()

    # Assess whether all required arguments are given
    if not (args.pyrad_filename or args.group_filename or
            args.re_site_filename or args.group1):
        sys.exit(('Script requires the following arguments: pyrad_file,'
                  ' group_file, restriction_site_file, and group1'))

    main(args.pyrad_filename, args.group_filename, args.re_site_filename,
         args.group1, args.group2, int(args.min_samples), args.output_folder)
