# RAD-seq script library
Collection of Python scripts for parsing/analysis of reduced representation sequencing data (e.g. RAD-seq, nextRAD). While many of the scripts are functional, some still need considerable cleaning up and more thorough testing - and this repository therefore very much represents a *work in progress*.

These scripts all require [Python 3](https://www.python.org/download/releases/3.0/), with some requiring additional packages ([BioPython](https://github.com/biopython/biopython.github.io/) and [NumPy](http://www.numpy.org/) - both of which can be easily installed using the [Miniconda](http://conda.pydata.org/miniconda.html) or [Anaconda](https://www.continuum.io/downloads) installers, or [PyVCF](https://github.com/jamescasbon/PyVCF) - which can be installed using e.g. `pip install PyVCF`). Usage information for each script can be obtained using the `-h` or `--help` flag (e.g. `python3 name_of_script.py -h`, or is also listed in this README.

This documentation is dynamically generated using the listed [README_compile.py](README_compile.py) script, extracting purpose, usage and links to example files from the [argparse](https://docs.python.org/3/library/argparse.html) information of each script.

## Recently added
**[vcf_remap2genome.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_remap2genome.py)** - script to remap VCF from de novo RAD assembly back to a reference genome

**[pyrad_find_caps_markers.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad_find_caps_markers.py)** - search PyRAD output file for diagnostic CAPS loci that can distinguish two groups (or one group and all other samples)

**[vcf_clone_detect.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_clone_detect.py)** - script to facilitate identification of clones in dataset

