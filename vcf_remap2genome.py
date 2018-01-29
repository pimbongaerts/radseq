#!/usr/bin/env python
"""
Remap VCF to genome.

* Currently only works with single-line FASTA files (but easy to change)
* Works best when using the optional alignment output in sam2tsv
"""
import sys
import argparse
import string

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2018 Pim Bongaerts'
__license__ = 'GPL'

# VCF file format constants
VCF_HEADER_CHAR = '#'
VCF_CHROM_COL = 0
VCF_POS_COL = 1
VCF_REF_COL = 3
VCF_ALT_COL = 4
VCF_FIRST_GT_COL = 9
# FASTA file format constants
FASTA_HEADER_CHAR = '>'
# SAMTSV file format constants
SAMTSV_HEADER_CHAR = '#'
SAMTSV_ALIGN_CHAR = ':'
SAMTSV_GAP_CHAR = '.'
SAMTSV_CHROM_COL = 0
SAMTSV_FLAG_COL = 1
SAMTSV_GENOMECHROM_COL = 2
SAMTSV_POS_COL = 3
SAMTSV_BASE_COL = 4
SAMTSV_GENOMEPOS_COL = 6
SAMTSV_GENOMEBASE_COL = 7
FORWARD_STRAND_FLAG = 0
REVERSE_STRAND_FLAG = 16
# Internal-use constants
IUPAC = {'A': ['A', 'A'], 'C': ['C', 'C'], 'G': ['G', 'G'], 'T': ['T', 'T'],
         'M': ['A', 'C'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['C', 'G'],
         'Y': ['C', 'T'], 'K': ['G', 'T']}
COMPL = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
              'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W',
              'M': 'K', 'K': 'M', 'N': 'N'}
DIRECTION = {-1: 'error', 0: 'forward', 16: 'reverse'}
ANNOTATORS = [x for x in string.ascii_lowercase]
NOT_FOUND = -1


class OriginalLoci(object):
    """ Dictionary containing the original loci (as in the vcf_file and
    fasta_file) and their SNPs as QueryLocus objects (which in turn contain Snp
     objects), with several methods go gather associated mapping data.
    """

    def __init__(self):
        self.loci = {}

    def get_snps_from_vcf(self, vcf_filename):
        """ Create the QueryLocus objects in OriginalLoci, using the locus
        names as given in the CHROM column from the VCF, and add the different
        SNPs (as indicated in the POS column from the VCF) as Snp objects.
        """
        vcf_file = open(vcf_filename, 'r')
        for line in vcf_file:
            if line[0] != VCF_HEADER_CHAR:
                cols = line.split()
                locus_name = cols[VCF_CHROM_COL]
                # Add query locus name (CHROM) if first SNP
                if locus_name not in self.loci:
                    self.loci[locus_name] = QueryLocus(locus_name)
                # Add SNP to query locus
                self.loci[locus_name].add_snp(int(cols[VCF_POS_COL]),
                                              cols[VCF_REF_COL],
                                              cols[VCF_ALT_COL])
        vcf_file.close()

    def add_seqs_from_fasta_to_original_loci(self, fasta_filename):
        """ Read sequences from the fasta file and store the in the
        corresponding QueryLocus objects
        """
        fasta_file = open(fasta_filename, 'r')
        locus_name = ''
        for line in fasta_file:
            if line[0] == FASTA_HEADER_CHAR:
                locus_name = line[1:].strip()
            else:
                if locus_name in self.loci:
                    # Add sequence
                    self.loci[locus_name].seq = line.strip()
                    # Verify for each SNP position, whether the fasta has one
                    # of the two alleles as defined in the VCF
                    if not self.loci[locus_name].verify_snps_in_fast_seq():
                        sys.exit('Error: allele mismatch between FASTA and '
                                 'VCF in locus {0}'.format(locus_name))
                else:
                    # Ignore fasta sequences with a name that is not in the
                    # VCF, as these were potentially removed during filtering
                    # of the VCF (and thus OK)
                    pass
        fasta_file.close()

    def add_mapping_info_from_samtsv(self, samtsv_filename, pid_threshold):
        """ Add mapping info from the samtsv file """
        samtsv_file = open(samtsv_filename, 'r')
        prev_locus_name = ''
        for line in samtsv_file:
            if line[0] not in [SAMTSV_HEADER_CHAR, SAMTSV_ALIGN_CHAR]:
                cols = line.split()
                locus_name = cols[SAMTSV_CHROM_COL]
                map_flag = int(cols[SAMTSV_FLAG_COL])
                # If original SNP represented a gap in the fasta: ignore,
                # (as site would not have been mapped to genome)
                if cols[SAMTSV_POS_COL] == SAMTSV_GAP_CHAR:
                    samtsv_pos = NOT_FOUND
                else:
                    samtsv_pos = int(cols[SAMTSV_POS_COL])
                # When arriving at a new locus in the samtsv (new locus_name)
                if locus_name != prev_locus_name:
                    # Create dict to cross-reference positions in the original
                    # vcf/fasta and the samtsv (different as samtsv may be
                    # based on a reverse-complement seq and without gaps)
                    cross_ref = self.loci[locus_name].crossref_samtsv(map_flag)
                    # Store orientation status of original locus and update
                    # alleles if required
                    genom_chrom = cols[SAMTSV_GENOMECHROM_COL]
                    self.loci[locus_name].set_ref_and_orientation(genom_chrom,
                                                                  map_flag)
                # When arriving at SNP save genome reference info
                if samtsv_pos in cross_ref.keys():
                    genom_base = cols[SAMTSV_GENOMEBASE_COL]
                    orig_pos = cross_ref[samtsv_pos]
                    if cols[SAMTSV_GENOMEPOS_COL] == SAMTSV_GAP_CHAR:
                        # SNP not present in genome reference (causing gap)
                        self.loci[locus_name].snps[orig_pos].delete()
                    else:
                        genom_pos = int(cols[SAMTSV_GENOMEPOS_COL])
                        (self.loci[locus_name].snps[orig_pos]
                             .set_ref(genom_pos, genom_base))
                        # Verify whether VCF alleles include allele in genome
                        if not (self.loci[locus_name]
                                    .snps[orig_pos]
                                    .verify_in_genome_seq()):
                            # If not, delete SNP from dataset
                            self.loci[locus_name].snps[orig_pos].delete()
                prev_locus_name = locus_name
            elif line[0] == SAMTSV_ALIGN_CHAR:
                # When arriving at alignment info (optional in sam2tsv output)
                # store that alignment information for current locus
                self.loci[locus_name].add_alignment_info(line, pid_threshold)
        samtsv_file.close()


class QueryLocus(object):
    """ Contains the locus name (as defined in VCF CHROM column), its sequence
    (as given in FASTA), and the corresponding SNPs (as Snp objects).
    """

    def __init__(self, name):
        self.name = name
        self.snps = {}
        self.seq = ''
        self.genom_chrom = ''
        self.orientation = -1
        self.samtsv_seq = ''
        self.alignment = []
        self.pid = 0

    def add_snp(self, orig_pos, ref_allele, alt_allele):
        """ Add SNP positions and alleles to QueryLocus """
        self.snps[orig_pos] = Snp(orig_pos,
                                  ref_allele,
                                  alt_allele,
                                  ANNOTATORS[len(self.snps)])

    def add_alignment_info(self, info, pid_threshold):
        """ Set sequence of Query Locus """
        self.alignment.append(info.strip())
        # Calculate percent ID match and evaluate when all lines are added
        if len(self.alignment) == 3:
            matches = self.alignment[1].count('|')
            self.pid = (matches / len(self.samtsv_seq)) * 100
            if pid_threshold and (self.pid < float(pid_threshold)):
                # Delete all SNPs of locus
                for snp_pos in self.snps:
                    self.snps[snp_pos].delete()

    def set_seq(self, seq):
        """ Set sequence of Query Locus """
        self.seq = seq

    def format_log(self):
        """ Return formatted string with locus info """
        info = []
        pid_info = ''
        original_snp_align, samtsv_snp_align = self.format_snp_alignment()
        # Calculate percent ID match if alignment is present in samtsv
        if len(self.alignment) == 3:
            pid_info = ' ({0:.2f} pid score)'.format(self.pid)
        # General locus and mapping information
        info.append('Locus {0} maps to genome scaffold/chrom '
                    '{1} in {2} direction{3}'.format(self.name,
                                                     self.genom_chrom,
                                                     DIRECTION[
                                                         self.orientation],
                                                     pid_info))
        # Original FASTA reference with SNPs from VCF
        info.append('*  Original FASTA:  1 {0} {1}'.format(self.seq,
                                                           len(self.seq)))
        info.append('*  original SNP pos:  {0}'.format(original_snp_align))
        # Mapping query with SNP locations
        info.append('*  Mapping query:   0 {0}'
                    ' {1}'.format(self.samtsv_seq, len(self.samtsv_seq) - 1))
        info.append('*  samtsv SNP pos:    {0}'.format(samtsv_snp_align))
        info.extend(self.alignment)
        info.append('\n')
        return '\n'.join(info)

    def set_ref_and_orientation(self, genom_chrom, map_flag):
        """ Set reference genome chrom/scaffold and orientation status of
        original locus in samtsv. Also set output alleles (i.e., complement
        if reverse mapped)
        """
        self.genom_chrom = genom_chrom
        self.orientation = map_flag
        # Infer sequence as referred to in samtsv
        seq_no_gaps = self.seq.replace('-', '')
        if map_flag == FORWARD_STRAND_FLAG:
            self.samtsv_seq = seq_no_gaps
        elif map_flag == REVERSE_STRAND_FLAG:
            self.samtsv_seq = self.__reverse_complement(seq_no_gaps)
        else:
            # Error in case of a different mapping flag
            sys.exit('Unexpected flag for locus {0}'.format(self.name))
        # Iterate over SNPs and set output_alleles
        for snp_pos in self.snps:
            if map_flag == FORWARD_STRAND_FLAG:
                new_ref = self.snps[snp_pos].ref_allele
                new_alt = self.snps[snp_pos].alt_allele
            elif map_flag == REVERSE_STRAND_FLAG:
                new_ref = COMPL[self.snps[snp_pos].ref_allele]
                new_alt = COMPL[self.snps[snp_pos].alt_allele]
            self.snps[snp_pos].output_ref_allele = new_ref
            self.snps[snp_pos].output_alt_allele = new_alt

    def verify_snps_in_fast_seq(self):
        """ Verify for each SNP position, whether the sequence (from FASTA)
        has one of the two alleles as defined in the VCF"""
        verifiable = True
        for snp_pos in self.snps:
            seq_allele = self.seq[snp_pos - 1]        # Potentially ambigious
            if seq_allele not in ('N', '-'):
                seq_alleles = IUPAC[seq_allele]     # List with non-ambiguous
                vcf_snp = self.snps[snp_pos]        # Alleles as coded in VCF
                for seq_allele in seq_alleles:
                    if seq_allele not in [vcf_snp.ref_allele,
                                          vcf_snp.alt_allele]:
                        verifiable = False      # Mismatch between VCF & fasta
                        break
        return verifiable

    def crossref_samtsv(self, mapping_flag):
        """ Returns cross-reference dict of format:
        cross_ref[samtsv_position] = original_position
        """
        cross_ref = {}
        for snp_pos in self.snps:
            if mapping_flag == FORWARD_STRAND_FLAG:
                # Count gaps from beginning of seq until SNP
                if snp_pos > 1:
                    gap_count = self.seq[0:snp_pos - 1].count('-')
                else:
                    gap_count = 0
                # [New SNP position] = [Old SNP position] - [gap count] - 1
                samtsv_pos = snp_pos - gap_count - 1
            elif mapping_flag == REVERSE_STRAND_FLAG:
                # Count gaps from SNP until end of seq
                rev_gap_count = self.seq[snp_pos:].count('-')
                # [New SNP position] = [original locus length including gaps]
                #                       - [Old SNP position]
                #                       - [reverse gap count]
                samtsv_pos = len(self.seq) - snp_pos - rev_gap_count
            else:
                # Error in case of a different mapping flag
                sys.exit('Unexpected flag for locus {0}'.format(self.name))
            self.snps[snp_pos].samtsv_pos = samtsv_pos
            cross_ref[samtsv_pos] = snp_pos
        return cross_ref

    def format_snp_alignment(self):
        """ Create a string with SNP positions indicated by asterisks """
        original_snp_align = [' ' for x in self.seq]
        samtsv_snp_align = [' ' for x in self.seq]
        for snp in self.snps.values():
            original_snp_align[snp.orig_pos - 1] = snp.alpha_id
            samtsv_snp_align[snp.samtsv_pos] = snp.alpha_id
        return ''.join(original_snp_align), ''.join(samtsv_snp_align)

    @staticmethod
    def __reverse_complement(seq):
        """ Returm reverse complement of sequence  """
        return ''.join([COMPL[base] for base in seq[::-1]])


class Snp(object):
    """ SNPs with their different corresponding coordinates """

    def __init__(self, orig_pos, ref_allele, alt_allele, alpha_id):
        self.orig_pos = orig_pos
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.alpha_id = alpha_id
        # Initialize
        self.nogap_pos = 0
        self.samtsv_pos = 0
        self.genom_pos = 0
        self.genom_base = ''
        self.output_ref_allele = ''
        self.output_alt_allele = ''
        self.deleted = False

    def set_ref(self, genom_pos, genom_base):
        """ Set sequence of Query Locus """
        self.genom_pos = genom_pos
        self.genom_base = genom_base

    def delete(self):
        """ Mark SNP as to be deleted """
        self.deleted = True

    def format_log(self):
        """ Return formatted string with SNP info for log file"""
        info = []
        additional_info = ''
        if self.deleted:
            additional_info = 'DELETED: '
        info.append('{0}SNP [{1}] in original sequence at pos {2} ({3}/{4}) '
                    'corresponds to mapping query position {5} ({6}/{7}) and '
                    'genome pos {8} ({9})\n'.format(additional_info,
                                                    self.alpha_id,
                                                    self.orig_pos,
                                                    self.ref_allele,
                                                    self.alt_allele,
                                                    self.samtsv_pos,
                                                    self.output_ref_allele,
                                                    self.output_alt_allele,
                                                    self.genom_pos,
                                                    self.genom_base))
        return '\n'.join(info)

    def format_table(self, chrom, genom_chrom):
        """ Return formatted string with SNP info for cross-ref table """
        info = []
        if self.deleted:
            status = 'deleted'
        else:
            status = 'included'
        return '{0}\n'.format('\t'.join([chrom, str(self.orig_pos),
                                         self.ref_allele, self.alt_allele,
                                         genom_chrom, str(self.genom_pos),
                                         self.genom_base,
                                         self.output_ref_allele,
                                         self.output_alt_allele, status]))

    def verify_in_genome_seq(self):
        """ Verify whether the genome reference base matches one of the two
        alleles as defined in the VCF, and change output reference allele
        to that of the genome"""
        verifiable = True
        if self.genom_base in ('N', '-', '.'):
            verifiable = False        # Mismatch between VCF & genome
        else:
            # List with non-ambiguous alleles
            seq_alleles = IUPAC[self.genom_base]
            for seq_allele in seq_alleles:
                if seq_allele not in [self.output_ref_allele,
                                      self.output_alt_allele]:
                    verifiable = False    # Mismatch between VCF & genome
        if verifiable and (self.genom_base in ['A', 'T', 'G', 'C']):
            if (self.genom_base != self.output_ref_allele and
                    self.genom_base == self.output_alt_allele):
                # Genome base corresponds to alternative allele in VCF:
                # change around reference and alternative
                self.output_alt_allele = self.output_ref_allele
                self.output_ref_allele = self.genom_base
            elif (self.genom_base == self.output_ref_allele and
                  self.genom_base != self.output_alt_allele):
                # Genome base corresponds to reference allele in VCF: leave
                pass
            else:
                sys.exit("Error: unexpected error in verify_in_genome_seq")
        return verifiable


def get_new_vcf_line(current_locus, current_snp, cols):
    """ Update line in VCF with updated information """
    if (current_locus.orientation == FORWARD_STRAND_FLAG and
        current_snp.output_ref_allele == current_snp.ref_allele and
            current_snp.output_alt_allele == current_snp.alt_allele):
        # Forward mapped & matching genome -> keep the same
        vcf_line = '\t'.join([current_locus.genom_chrom,
                              str(current_snp.genom_pos),
                              '\t'.join(cols[2:])])
    elif ((current_locus.orientation == FORWARD_STRAND_FLAG and
           current_snp.output_ref_allele == current_snp.alt_allele and
           current_snp.output_alt_allele == current_snp.ref_allele and
           current_snp.genom_base == current_snp.alt_allele) or
          # Forward mapped & not matching genome -> change both ref/alt and gts
          (current_locus.orientation == REVERSE_STRAND_FLAG and
           current_snp.genom_base == COMPL[current_snp.alt_allele] and
           current_snp.output_ref_allele == COMPL[current_snp.alt_allele] and
           current_snp.output_alt_allele == COMPL[current_snp.ref_allele])):
        # Reverse mapped & not matching genome -> change both ref/alt and gts
        genotypes = cols[VCF_FIRST_GT_COL:]
        new_genotypes = []
        for genotype in genotypes:
            if genotype[0] != '.':
                new_alleles = sorted(
                    [1 - int(genotype[0]), 1 - int(genotype[2])])
                new_genotypes.append('{0}{1}{2}{3}'.format(new_alleles[0],
                                                           genotype[1],
                                                           new_alleles[1],
                                                           genotype[3:]))
            else:
                new_genotypes.append(genotype)
        vcf_line = '\t'.join([current_locus.genom_chrom,
                              str(current_snp.genom_pos),
                              cols[2],
                              current_snp.output_ref_allele,
                              current_snp.output_alt_allele,
                              '\t'.join(cols[5:VCF_FIRST_GT_COL]),
                              '\t'.join(new_genotypes)])
    elif (current_locus.orientation == REVERSE_STRAND_FLAG and
          current_snp.genom_base == COMPL[current_snp.ref_allele] and
          current_snp.output_ref_allele == COMPL[current_snp.ref_allele] and
          current_snp.output_alt_allele == COMPL[current_snp.alt_allele]):
        # Reverse mapped & matching genome base -> ref/alt columns changed only
        vcf_line = '\t'.join([current_locus.genom_chrom,
                              str(current_snp.genom_pos),
                              cols[2],
                              current_snp.output_ref_allele,
                              current_snp.output_alt_allele,
                              '\t'.join(cols[5:])])
    else:
        print('\t'.join(cols))
        sys.exit("Error: unexpected error in new_vcf_line")
    return vcf_line


def output_new_vcf(original_loci, vcf_filename, output_vcf_filename):
    """ Output new VCF with updated positions """
    vcf_file = open(vcf_filename, 'r')
    output_vcf_file = open(output_vcf_filename, 'w')
    output_log_file = open('remap_log.txt', 'w')
    output_table_file = open('remap_table.txt', 'w')
    prev_chrom = ''
    for line in vcf_file:
        if line[0] == VCF_HEADER_CHAR:
            output_vcf_file.write(line)
        else:
            cols = line.split()
            # Current/old CHROM/POS in vcf
            vcf_chrom = cols[VCF_CHROM_COL]
            vcf_pos = int(cols[VCF_POS_COL])
            # Retrieve that locus from the original_loci dict
            current_locus = original_loci.loci[vcf_chrom]
            current_snp = current_locus.snps[vcf_pos]
            # Log: output locus information when reaching a new CHROM
            if vcf_chrom != prev_chrom:
                if prev_chrom != '':
                    output_log_file.write('\n')
                output_log_file.write(current_locus.format_log())
            # Log/table: Output SNP information for each SNP
            output_log_file.write(current_snp.format_log())
            output_table_file.write(current_snp
                                    .format_table(current_locus.name,
                                                  current_locus.genom_chrom))
            # VCF: output updated information
            if (not current_locus.snps[vcf_pos].deleted and
                    not current_locus.orientation == NOT_FOUND):
                new_vcf_line = get_new_vcf_line(current_locus,
                                                current_snp,
                                                cols)
                output_vcf_file.write('{0}\n'.format(new_vcf_line))
            prev_chrom = vcf_chrom
    vcf_file.close()
    output_vcf_file.close()
    output_log_file.close()
    output_table_file.close()


def main(vcf_filename, fasta_filename, samtsv_filename,
         output_vcf_filename, pid_threshold):
        # Initialise object holding the loci, SNPs and mapping info
    original_loci = OriginalLoci()

    # Fill original_loci with locus_names and SNP positions
    original_loci.get_snps_from_vcf(vcf_filename)

    # Add sequences to QueryLocus objects (in dict) from fasta file
    original_loci.add_seqs_from_fasta_to_original_loci(fasta_filename)

    # Get mapping info from samtsv file
    original_loci.add_mapping_info_from_samtsv(samtsv_filename, pid_threshold)

    # Output new VCF with updated positions
    output_new_vcf(original_loci, vcf_filename, output_vcf_filename)


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--vcf', dest='vcf_filename',
                        metavar='vcf_file',
                        help='original vcf file, where the CHROM correspond \
                        to the sequence name in the supplied fasta (perfect \
                        match), and the POS the position within that sequence \
                        (also counting gaps if present).')
    parser.add_argument('-f', '--fasta', dest='fasta_filename',
                        metavar='fasta_file',
                        help='fasta file representing a single sequence for \
                        each locus/CHROM in the original vcf (including gaps \
                        if present in the original aligment that the fasta \
                        is based on). Note that these gaps need to be removed \
                        before mapping these sequences back to the genome \
                        with bwa mem, but need to be present in this file.')
    parser.add_argument('-t', '--t', dest='samtsv_filename',
                        metavar='samtsv_file',
                        help='The sam file converted to tsv. The sam file \
                        represents the mapping outcome of the supplied fasta \
                        (but then without gaps) to the genome with bwa mem, \
                        e.g. through: bwa mem ref_genome \
                        fasta_file_no_gaps.fa > samfile.sam. This samfile \
                        then needs to be converted to a tsv, using sam2tsv \
                        from the jVarKit toolkit \
                        (http://lindenb.github.io/jvarkit/Sam2Tsv.html).')
    parser.add_argument('-o', '--output_vcf', dest='output_vcf_filename',
                        metavar='output_vcf_file',
                        help='remapped vcf file with genome scaffold/chroms \
                        and positions within those scaffold/chroms.')
    parser.add_argument('-pid', '--pid', dest='pid_threshold',
                        metavar='pid_threshold',
                        help='optional pid alignment threshold, to exclude \
                        loci aligning to the genome with a percent id (PID) \
                        score below the indicated value.')
    args = parser.parse_args()

    # Assess whether all arguments are given (otherwise exit):
    if not (args.vcf_filename or args.fasta_filename or
            args.samtsv_filename or args.output_vcf_filename):
        sys.exit(('Script requires the following arguments: vcf_file,'
                  ' fasta_file, samtsv_filename, and output_vcf_file'))

    main(args.vcf_filename, args.fasta_filename, args.samtsv_filename,
         args.output_vcf_filename, args.pid_threshold)
