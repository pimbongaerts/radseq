# RAD-seq script library
Collection of Python scripts for parsing/analysis of reduced representation sequencing data (e.g. RAD-seq, nextRAD). While many of the scripts are functional, some still need considerable cleaning up and more thorough testing - and this repository therefore very much represents a *work in progress*.

These scripts all require [Python 3](https://www.python.org/download/releases/3.0/), with some requiring additional packages ([BioPython](https://github.com/biopython/biopython.github.io/) and [NumPy](http://www.numpy.org/) - both of which can be easily installed using the [Miniconda](http://conda.pydata.org/miniconda.html) or [Anaconda](https://www.continuum.io/downloads) installers, or [PyVCF](https://github.com/jamescasbon/PyVCF) - which can be installed using e.g. `pip install PyVCF`). Usage information for each script can be obtained using the `-h` or `--help` flag (e.g. `python3 name_of_script.py -h`, or is also listed in this README.

This documentation is dynamically generated using the listed [README_compile.py](README_compile.py) script, extracting purpose, usage and links to example files from the [argparse](https://docs.python.org/3/library/argparse.html) information of each script.

## structure

**[structure_mp.py](structure_mp.py)** - Multi-processing STRUCTURE (Pritchard et al 2000) wrapper for RAD-seq data.
Takes a `.vcf` as input file and then creates a number of replicate datasets,
each with a different pseudo-random subsampling of one SNP per RAD contig.
Then, it runs the replicate datasets through STRUCTURE across multiple
threads, and summarises the outcome with CLUMPP (Jakobsson and Rosenberg
2007). Finally, it assesses the number of potential clusters using the
Puechmaille 2016 method (only suitable for certain datasets). Note: STILL
NEEDS TO BE MODIFIED FOR GENERAL USE. Input file (`.vcf`) should be sorted by
CHROM, and different params files need to be present in current path.

	usage: structure_mp.py [-h] vcf_file pop_file maxK replicates threads

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	  pop_file    population file (.txt)
	  maxK        maximum number of K (expected clusters)
	  replicates  number of replicate runs for each K
	  threads     number of parallel threads
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [pop_file.txt](input_examples/pop_file.txt).



## other

**[fasta_include.py](fasta_include.py)** - Reduces FASTA file to only those loci listed in supplied text file.

	usage: fasta_include.py [-h] fasta_file inclusion_file

	positional arguments:
	  fasta_file      FASTA input file (`.fasta`/ `.fa`)
	  inclusion_file  text file with names of loci to be included
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [fasta_file.fa](input_examples/fasta_file.fa), [inclusion_file.txt](input_examples/inclusion_file.txt).


**[README_compile.py](README_compile.py)** - Compiles README markdown file for this repository
(https://github.com/pimbongaerts/radseq). Categories are assigned based on
prefix, usage information is extracted from argparse, and example input files
are assigned based on argument names. *[File did not pass PEP8 check]*

	usage: README_compile.py [-h]

	optional arguments:
	  -h, --help  show this help message and exit
	





## fastq

**[fastq_barcodes2samplenames.py](fastq_barcodes2samplenames.py)** - Renames barcoded .fastq.gz files in a folder to sample/individual names. Takes
relative or absolute path as first argument (e.g. 'samples'; without trailing
slash) and a text file as second argument. The latter should be a tab-
separated text file, with the barcode in the first column, and the new
sample/individual name in the second column. The script conducts a trial run
first, listing the files to be renamed, and then asks for confirmation before
doing the actual renaming.

	usage: fastq_barcodes2samplenames.py [-h] path barcode_file

	positional arguments:
	  path          path (for current directory use `.`)
	  barcode_file  text file (tsv or csv) with barcodes and sample names
	
	optional arguments:
	  -h, --help    show this help message and exit
	




**[fastq_bin_paired_reads.py](fastq_bin_paired_reads.py)** - Clusters reads of paired-end RAD-seq data for downstream contig assembly. It
maps R1 reads to a reference, and then outputs those reads and the
corresponding R2 reads to a separate 'shuffled' FASTQ file per locus. Note:
when using an existing output folder, reads are being appended to existing
files (use this to append data from multiple samples). BWA needs to be
installed and accessible through `PATH` environmental variable.

	usage: fastq_bin_paired_reads.py [-h]
                                 r1_fastq_file r2_fastq_file ref_fasta_file
                                 threads output_folder

	positional arguments:
	  r1_fastq_file   file in FASTQ format with R1 reads
	  r2_fastq_file   file in FASTQ format with R2 reads
	  ref_fasta_file  file in FASTA format with reference contigs
	  threads         number of threads to be used by BWA
	  output_folder   name of output folder
	
	optional arguments:
	  -h, --help      show this help message and exit
	




**[fastq_seqcount.py](fastq_seqcount.py)** - Outputs number of reads for each fastq.gz sample to text file, and prints
mean/min/max to STDOUT. Note: Does not work with all FASTQ formats, and
correct depending on OS `zcat` or `gzcat` needs to be set in COMPRESS_UTIL
constant.

	usage: fastq_seqcount.py [-h] path output_filename

	positional arguments:
	  path             path (for current directory use `.`)
	  output_filename  name of text output file
	
	optional arguments:
	  -h, --help       show this help message and exit
	





## popfile

**[popfile_from_clusters.py](popfile_from_clusters.py)** - Outputs a popfile based on cluster assignment file (from e.g. STRUCTURE) and
outputs a popfile based on those assigments (using supplied assignment
threshold). Note: I use the formatted CLUMPP output (e.g. `clumpp_K4.out.csv`)
from the `structure_mp` wrapper as assignment file.

	usage: popfile_from_clusters.py [-h] assignment_file assign_cut_off

	positional arguments:
	  assignment_file  text file (tsv or csv) with assignment values for each
	                   individual (max. 2 clusters); e.g. a reformatted STRUCTURE
	                   output file
	  assign_cut_off   min. assignment value for an individual to be assigned to a
	                   cluster
	
	optional arguments:
	  -h, --help       show this help message and exit
	

Example input file(s):  [assignment_file.csv](input_examples/assignment_file.csv).


**[popfile_randomize.py](popfile_randomize.py)** - Pseudo-randomize the assignment of individuals to populations. Note: with
small population sizes, this can lead to very uneven simulated population
sizes. See also the alternative: `popfile_toggleassign.py`. Individuals in
original popfile should be ordered by population.

	usage: popfile_randomize.py [-h] pop_file

	positional arguments:
	  pop_file    text file (tsv or csv) with individuals and populations
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [pop_file.txt](input_examples/pop_file.txt).


**[popfile_toggleassign.py](popfile_toggleassign.py)** - Shuffle the assignment of individuals to populations by assigning the indivs
sequentially to the different pops. The assignment is not completely random,
but does generate equal population sizes (which otherwise differ substantially
when using random assignment under originally small population sizes).

	usage: popfile_toggleassign.py [-h] pop_file

	positional arguments:
	  pop_file    text file (tsv or csv) with individuals and populations
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [pop_file.txt](input_examples/pop_file.txt).



## pyrad

**[pyrad2concat_fasta.py](pyrad2concat_fasta.py)** - Concatenates PyRAD/ipyrad sequences from `.loci` file for each individual and
outputs as FASTA (order by popfile). Note: missing data are filled with gaps
(`N`)

	usage: pyrad2concat_fasta.py [-h] pyrad_file pop_file

	positional arguments:
	  pyrad_file  PyRAD file (`.loci`)
	  pop_file    text file (tsv or csv) with individuals and populations
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [pyrad_file.loci](input_examples/pyrad_file.loci), [pop_file.txt](input_examples/pop_file.txt).


**[pyrad2fasta.py](pyrad2fasta.py)** - Create FASTA file with a representative sequence (using first sample) for each
locus in pyRAD/ipyrad `.loci` or `.allele` file.

	usage: pyrad2fasta.py [-h] pyrad_file

	positional arguments:
	  pyrad_file  PyRAD allele file (`.loci` or `.allele`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [pyrad_file.loci](input_examples/pyrad_file.loci).


**[pyrad2migrate.py](pyrad2migrate.py)** - Converts PyRAD `.allele` file to migrate-n input file (population designated
indicated in supplied popfile). Note: only appropriate for PyRAD `.allele`
file (not `.loci`).

	usage: pyrad2migrate.py [-h] allele_file pop_file

	positional arguments:
	  allele_file  PyRAD allele file (.allele)
	  pop_file     text file (tsv or csv) with individuals and populations
	
	optional arguments:
	  -h, --help   show this help message and exit
	

Example input file(s):  [pop_file.txt](input_examples/pop_file.txt).


**[pyrad_filter.py](pyrad_filter.py)** - Filters PyRAD output file (`.loci`) for those loci (1) present or absent
(using --exclude flag) in supplied list, (2) genotyped for at least X number
of samples, and (3) with at least Y number of informative sites. Note: can
also be used for `.alleles` file but then 2x the number of samples should be
given (assuming diploid individual).

	usage: pyrad_filter.py [-h] [-e]
                       pyrad_file loci_file sample_threshold snp_threshold

	positional arguments:
	  pyrad_file        PyRAD file (`.loci`)
	  loci_file         text file with PyRAD loci to be included
	  sample_threshold  min. number of samples genotyped for a locus to be
	                    included
	  snp_threshold     min. number of SNPs for a locus to be included
	
	optional arguments:
	  -h, --help        show this help message and exit
	  -e, --exclude     use the loci in loci_file as exclusion list
	

Example input file(s):  [pyrad_file.loci](input_examples/pyrad_file.loci), [loci_file.txt](input_examples/loci_file.txt).


**[pyrad_trim.py](pyrad_trim.py)** - Trims sequence length in PyRAD/ipyrad `.alleles` or `.loci` file.

	usage: pyrad_trim.py [-h] pyrad_file seq_length

	positional arguments:
	  pyrad_file  PyRAD allele file (`.loci` or `.allele`)
	  seq_length  length to which sequences are trimmed
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [pyrad_file.loci](input_examples/pyrad_file.loci).


**[pyradclust2fasta.py](pyradclust2fasta.py)** - Creates one large FASTA from all PyRAD clustS files in directory. Only outputs
clusters that exceed size threshold (min. number of sequences in cluster).
First sequence of each cluster is outputted (together with size of overall
cluster - note: not of that specific sequence). Prints the outputted and total
number of clusters to STDOUT.

	usage: pyradclust2fasta.py [-h] path cluster_threshold output_file

	positional arguments:
	  path               path that contains PyRAD `.clustS` files
	  cluster_threshold  minimum size of cluster to be included
	  output_file        name of output FASTA file
	
	optional arguments:
	  -h, --help         show this help message and exit
	





## mapping

**[mapping_get_blastn_matches.py](mapping_get_blastn_matches.py)** - Extracts a list of loci that have a blastn e-value below a certain threshold,
and outputs the (first) matching reference locus, as well as the alignment
length, nident, and pident. Results are outputted to file with the chosen
e-value as post-fix, and STDOUT gives minimum alignment stats for filtered
loci. Note: fields in input file should be (in this order): query id, subject
id, alignment length, identity, perc. identity, evalue, bitscore (additional
fields beyond that are fine). This can be achieved by using `blastn` with the
following argument: `-outfmt 7 qseqid sseqid length nident pident evalue
bitscore`.

	usage: mapping_get_blastn_matches.py [-h] blastn_file evalue_cut_off

	positional arguments:
	  blastn_file     blastn output file with the following fields (in that
	                  order): query id, subject id, alignment length, identity,
	                  perc. identity, evalue, bitscore
	  evalue_cut_off  maximum e-value for match to be included
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [blastn_file.txt](input_examples/blastn_file.txt).


**[mapping_get_bwa_matches.py](mapping_get_bwa_matches.py)** - Extracts a list of succesfully mapped loci from `.sam` file (produced with
`bwa mem`). Successfully mapped loci are identified by default identified as
those with flags 0 and 16 (can be adjusted in MATCH_FLAGS constant).
Configured for use with single-end reads.

	usage: mapping_get_bwa_matches.py [-h] sam_file

	positional arguments:
	  sam_file    `bwa mem` output file (`.sam`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [sam_file.sam](input_examples/sam_file.sam).


**[mapping_identify_blast_matches.py](mapping_identify_blast_matches.py)** - Extracts a list of loci that have a blastn e-value below a certain threshold,
and outputs the (first) matching reference locus, as well as the alignment
length, nident, e-value and bitscore. It also compiles a set of all tax_ids,
which it uses to connect with the NCBI taxonomy database to get phylum ids for
each match using Entrez. Results are outputted to file with the chosen e-value
as post-fix, and STDOUT gives minimum alignment stats for filtered loci. Note:
fields in input file should be (in this order): query id, subject id,
alignment length, identity, perc. identity, evalue, bitscore, staxids, stitle.
This can be achieved by using `blastn` with the following argument: `-outfmt 7
qseqid sseqid length nident pident evalue bitscore staxids stitle`.

	usage: mapping_identify_blast_matches.py [-h] blastn_file evalue_cut_off email

	positional arguments:
	  blastn_file     blastn output file with the following fields (in that
	                  order): query id, subject id, alignment length, identity,
	                  perc. identity, evalue, bitscore
	  evalue_cut_off  maximum e-value for match to be included
	  email           email address to be used for NCBI connection
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [blastn_file.txt](input_examples/blastn_file.txt).



## vcf

**[vcf2hapmatrix.py](vcf2hapmatrix.py)** - Converts `.vcf` file to Tag Haplotype Matrix (with Chrom), with order of
individuals as indicated in optional file. Note: not yet properly tested. SNPs
of same CHROM (first column) in `.vcf` should be grouped
together/sequentially, and all individuals need to be listed in order_file.

	usage: vcf2hapmatrix.py [-h] [-o order_file] vcf_file

	positional arguments:
	  vcf_file              input file with SNP data (`.vcf`)
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -o order_file, --order_file order_file
	                        text file with preferred output order of individuals
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf2introgress.py](vcf2introgress.py)** - Converts `.vcf` file to INTROGRESS input files (4 files). Splits data into
three categories: parental1, parental2 and admixed based on cluster assignment
(provided in separate file; e.g. STRUCTURE output) and given threshold, and
outputs data for loc that exceed a certain frequency difference between the
two 'parental' categories. Note: not yet properly tested. also see similar
`vcf_ancestry_matrix.py` script. I use the formatted CLUMPP output
(`clumpp_K2.out.csv`) from the `structure_mp` wrapper as assignment file (max.
of 2 clusters).

	usage: vcf2introgress.py [-h] [--include]
                         vcf_file assignment_file assign_cut_off freq_cut_off
                         output_prefix

	positional arguments:
	  vcf_file         input file with SNP data (`.vcf`)
	  assignment_file  text file (tsv or csv) with assignment values for each
	                   individual (max. 2 clusters); e.g. a reformatted STRUCTURE
	                   output file
	  assign_cut_off   min. assignment value for an individual to be included in
	                   the allele frequency calculation (i.e. putative purebred)
	  freq_cut_off     min. allele frequency difference between the 2 clusters for
	                   a locus to be included in the output
	  output_prefix    prefix for output files
	
	optional arguments:
	  -h, --help       show this help message and exit
	  --include, -i    set this flag if parental pops need to be included in
	                   output
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [assignment_file.csv](input_examples/assignment_file.csv).


**[vcf2outflank.py](vcf2outflank.py)** - Converts `.vcf` file to OutFLANK input files (3 different files). Note: not
yet properly tested.

	usage: vcf2outflank.py [-h] vcf_file pop_file output_prefix

	positional arguments:
	  vcf_file       input file with SNP data (`.vcf`)
	  pop_file       text file (tsv or csv) with individuals and populations
	  output_prefix  prefix for output files
	
	optional arguments:
	  -h, --help     show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [pop_file.txt](input_examples/pop_file.txt).


**[vcf_ancestry_matrix.py](vcf_ancestry_matrix.py)** - Creates a genotype matrix for loci that have a large allele frequency
difference between two genetic clusters (as identified with e.g. STRUCTURE).
The script takes both a `.vcf` file and a text file with the assignment
probabilities as input. An assignment threshold (e.g. 0.98) needs to be
supplied to identify the reference individuals in the two clusters, and an
allele frequency cut-off needs to be supplied to identify divergent loci. An
optional file can be supplied with a list of loci that need to be included
regardless (e.g. previously identified outliers). Note: I use the formatted
CLUMPP output (`clumpp_K2.out.csv`) from the `structure_mp` wrapper as
assignment file (max. of 2 clusters).

	usage: vcf_ancestry_matrix.py [-h] [--include inclusion_file]
                              vcf_file assignment_file assign_cut_off
                              freq_cut_off

	positional arguments:
	  vcf_file              input file with SNP data (`.vcf`)
	  assignment_file       text file (tsv or csv) with assignment values for each
	                        individual (max. 2 clusters); e.g. a reformatted
	                        STRUCTURE output file
	  assign_cut_off        min. assignment value for an individual to be included
	                        in the allele frequency calculation (i.e. putative
	                        purebred
	  freq_cut_off          min. allele frequency difference between the 2
	                        clusters for a locus to be included in the output
	
	optional arguments:
	  -h, --help            show this help message and exit
	  --include inclusion_file, -i inclusion_file
	                        text file with loci to be included in output
	                        regardless of allele frequency differences
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [assignment_file.csv](input_examples/assignment_file.csv).


**[vcf_contrast_samples.py](vcf_contrast_samples.py)** - Contrast all samples in `.vcf` file against certain reference sample(s) (e.g.
outgroup samples), to assess for fixed / private alleles.

	usage: vcf_contrast_samples.py [-h]
                               vcf_file
                               [reference_samples [reference_samples ...]]

	positional arguments:
	  vcf_file           input file with SNP data (`.vcf`)
	  reference_samples  sample(s) against which the remainder of the dataset will
	                     be compared
	
	optional arguments:
	  -h, --help         show this help message and exit
	




**[vcf_find_clones.py](vcf_find_clones.py)** - Script compares the allelic similarity of individuals in a VCF, and outputs
all pairwise comparisons. This can be used to detect potential clones based on
percentage match. Note: highest matches can be assessed in the output file by
using `$ sort -rn --key=5 output_file.txt | head -n 50` in the terminal.

	usage: vcf_find_clones.py [-h] vcf_file

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf_gdmatrix.py](vcf_gdmatrix.py)** - Calculates Genetic Distance (Hamming / p-distance) for each pair of
individuals in a `.vcf` file and outputs as matrix. Popfile is supplied to
indicate order in matrix.

	usage: vcf_gdmatrix.py [-h] vcf_file pop_file

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	  pop_file    text file (tsv or csv) with individuals and populations
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [pop_file.txt](input_examples/pop_file.txt).


**[vcf_genotype_freqs.py](vcf_genotype_freqs.py)** - Outputs genotype frequencies for specific SNPs in each population, organised
by group.

	usage: vcf_genotype_freqs.py [-h] vcf_file factor_file SNP_file

	positional arguments:
	  vcf_file     input file with SNP data (`.vcf`)
	  factor_file  text file (tsv or csv) with individuals, their population
	               assignment and group assignment
	  SNP_file     text file (tsv or csv) with CHROM/POS of each SNP to be
	               outputted
	
	optional arguments:
	  -h, --help   show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [factor_file.txt](input_examples/factor_file.txt).


**[vcf_get_chrom_pos_from_number.py](vcf_get_chrom_pos_from_number.py)** - Translates sequential marker numbers back to CHROM/POS from original `.vcf`
file. Several programs only allow for integers to identify markers, this
script is to restore the original CHROM/POS for markers that were identified.

	usage: vcf_get_chrom_pos_from_number.py [-h] vcf_file markernumbers_file

	positional arguments:
	  vcf_file            input file with SNP data (`.vcf`)
	  markernumbers_file  text file with SNP numbers that were identified
	
	optional arguments:
	  -h, --help          show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [markernumbers_file.txt](input_examples/markernumbers_file.txt).


**[vcf_include_chrom.py](vcf_include_chrom.py)** - Retains only those loci (/CHROMs) in `.vcf` that are given in file.

	usage: vcf_include_chrom.py [-h] vcf_file inclusion_file

	positional arguments:
	  vcf_file        input file with SNP data (`.vcf`)
	  inclusion_file  text file with loci to be retained
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [inclusion_file.txt](input_examples/inclusion_file.txt).


**[vcf_minrep_filter.py](vcf_minrep_filter.py)** - Filters `.vcf` file for SNPs that are genotyped for a minimum proportion of
individuals in each of the populations (rather than overall proportion of
individuals). This can help to guarantee a minimum number of individuals to
calculate population-based statistics, and eliminate loci that might be
suffering from locus drop-out in particular populations. Note: only
individuals that are listed in popfile are taken into account to determine
proportion of individuals genotyped (but all indivs are outputted).

	usage: vcf_minrep_filter.py [-h]
                            vcf_file pop_file min_proportion output_filename

	positional arguments:
	  vcf_file         input file with SNP data (`.vcf`)
	  pop_file         text file (tsv or csv) with individuals and populations
	  min_proportion   proportion of individuals required to be genotyped in each
	                   population for a SNP to be included (e.g `0.8` for 80
	                   percent of individuals)
	  output_filename  name of output file (`.vcf`)
	
	optional arguments:
	  -h, --help       show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [pop_file.txt](input_examples/pop_file.txt).


**[vcf_missing_data.py](vcf_missing_data.py)** - Outputs list of missing data (# and % of SNPs) for each sample in VCF, to
identify poor-performing samples to eliminate prior to SNP filtering. Takes
vcf_filename as argument. Outputs to STDOUT (no output file).

	usage: vcf_missing_data.py [-h] vcf_file

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf_pos_count.py](vcf_pos_count.py)** - Counts SNP occurrence frequency for each POS in `.vcf` file.

	usage: vcf_pos_count.py [-h] vcf_file

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf_read_trim.py](vcf_read_trim.py)** - Removes SNPs from `.vcf` that are above a certain POS value.

	usage: vcf_read_trim.py [-h] vcf_file highest_pos

	positional arguments:
	  vcf_file     input file with SNP data (`.vcf`)
	  highest_pos  max. POS value allowed in `.vcf`
	
	optional arguments:
	  -h, --help   show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf_reference_loci.py](vcf_reference_loci.py)** - Lists all loci (using CHROM column) in `.vcf` that are genotyped for at least
one of the indicated samples/individuals. This can be used to reduce the
dataset to loci matching an included reference (e.g. aposymbiotic) samples.
Note: `vcf` can subsequently be filtered by using the output as inclusion_file
for `vcf_include_chrom.py`.

	usage: vcf_reference_loci.py [-h]
                             vcf_file
                             [reference_samples [reference_samples ...]]

	positional arguments:
	  vcf_file           input file with SNP data (`.vcf`)
	  reference_samples  sample(s) against which the remainder of the dataset will
	                     be compared
	
	optional arguments:
	  -h, --help         show this help message and exit
	




**[vcf_remove_chrom.py](vcf_remove_chrom.py)** - Excludes those loci (/CHROMs) in `.vcf` that are listed in exclusion list.
Also outputs a logfile with loci that were listed but not present in `.vcf`.

	usage: vcf_remove_chrom.py [-h] vcf_file exclusion_file

	positional arguments:
	  vcf_file        input file with SNP data (`.vcf`)
	  exclusion_file  text file loci to be excluded
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [exclusion_file.txt](input_examples/exclusion_file.txt).


**[vcf_rename_loci.py](vcf_rename_loci.py)** - Renames CHROMS in `.vcf` file according to list with old/new names, and only
outputs those loci that are listed.

	usage: vcf_rename_loci.py [-h] vcf_file locusnames_file

	positional arguments:
	  vcf_file         input file with SNP data (`.vcf`)
	  locusnames_file  text file (tsv or csv) with old and new name for each locus
	                   (/CHROM)
	
	optional arguments:
	  -h, --help       show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [locusnames_file.txt](input_examples/locusnames_file.txt).


**[vcf_rename_samples.py](vcf_rename_samples.py)** - Renames sample in `.vcf` file according to list with old/new names; also
outputs samples that are not listed in name_change file.

	usage: vcf_rename_samples.py [-h] vcf_file samplenames_file

	positional arguments:
	  vcf_file          input file with SNP data (`.vcf`)
	  samplenames_file  text file (tsv or csv) with old and new name for each
	                    sample (not all samples have to be listed)
	
	optional arguments:
	  -h, --help        show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [samplenames_file.txt](input_examples/samplenames_file.txt).


**[vcf_single_snp.py](vcf_single_snp.py)** - Reduces `.vcf` file to a single 'random' SNP per locus/chrom. Use for analyses
that require SNPs that are not physically linked. (although note that they of
course still may be - particularly so when dealing with short loci)

	usage: vcf_single_snp.py [-h] vcf_file

	positional arguments:
	  vcf_file    input file with SNP data (`.vcf`)
	
	optional arguments:
	  -h, --help  show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf).


**[vcf_spider.py](vcf_spider.py)** - Wrapper for PGDspider on Mac OS to convert `.vcf` files to various formats.
Note : set PGDSPIDER_PATH constant before use, and make script executable in
terminal with `$ chmod +x vcf_spider.py`.

	usage: vcf_spider.py [-h] vcf_filename pop_filename output_filename

	positional arguments:
	  vcf_filename     original vcf file
	  pop_filename     pop filename (.txt)
	  output_filename  output filename (extension used to determine file format
	                   (.genepop, .bayescan, .structure or .arlequin)
	
	optional arguments:
	  -h, --help       show this help message and exit
	




**[vcf_splitfst.py](vcf_splitfst.py)** - Filter original SNP dataset (`.vcf`) for a particular Fst percentile bin.
Note: order in Fst file needs to correspond with (`.vcf`) file, currently set
(see script CONSTANTS) to work with LOSITAN output file, and output filename
automatically generated from percentile bins.

	usage: vcf_splitfst.py [-h] vcf_file fst_file min_percentile max_percentile

	positional arguments:
	  vcf_file        input file with SNP data (`.vcf`)
	  fst_file        text file (tsv or csv) with Fst values for each SNP (same
	                  order as in vcf) - currently set to work with LOSITAN output
	                  file
	  min_percentile  min. Fst value for a SNP to be included
	  max_percentile  max. Fst value for a SNP to be included
	
	optional arguments:
	  -h, --help      show this help message and exit
	

Example input file(s):  [vcf_file.vcf](input_examples/vcf_file.vcf), [fst_file.txt](input_examples/fst_file.txt).


