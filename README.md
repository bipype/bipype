# bipype
bipype stands for BioInformatics-PYthon-PipE.


## Available tools
* fastx-toolkit v0.0.13
* usearch v7.0.959_i86linux32
* FLASH v1.2.7
* bowtie2 v2.2.4
* samtools v0.1.18
* RAPsearch 2.12_64bits
* MetaVelvet v1.2.01
* HUMAnN v0.99
* SARTools v1.2.0
* KRONA 2.0


## Available workflows
Bipype accepts three types of input: amplicons, whole genome sequences (WGS)
and metatranscriptomic data. Bipype works with both paired-end and single read
sequences. For amplicon (prokaryotic or fungal) data, if needed, paired-end
reads are merged and sequence reconstruction is performed. Reconstructed 16S or
ITS sequences are searched in proper reference databases and hits to related
taxonomic units are counted.

For WGS reads three paths are available. In one they will be used for
reconstruction of contigs, which will be further used for reference database
search (outside of the pipeline). In the second they will be used directly in
taxonomic diversity search, in third they will be compared to sequences related
to metabolic pathways.

The biodiversity data at taxonomical and functional levels are similar in their
tree structure - we provide display of results in a common form of an
interactive HTML tree.


## Installation
bipype is written in Python 2.7.
To install it, clone the repository with command
`git clone https://github.com/bipype/bipype.git`

## Usage
Chack "bipype -h" for details about running bipype
