==============================
DESCRIPTION OF PROGRAM OPTIONS
==============================

Bipype is a very useful program, which prepare a lot of types of bioinformatics analyses.
There are three input options: amplicons, WGS (whole genome sequences) and metatranscriptomic data. If amplicons are input data, then bipype does reconstruction and pairs merging. After that biodiversity is searching. There are two types of searching depending on the amplicons types (ITS or 16S).
If WGS are chosen, then bipype finds the SA coordinates of the input reads and generates alignments in the SAM format given single-end reads, aligns reads to reference sequence(s).
All of these analyses will be shown with Krona program, which allows to show hierarchical data with pie charts.


bipype stands for bioinformatics-python-pipe.


------------------------
Proposed analysis
------------------------

Below are listed a few analysis which could be performed with bipype program.
As input for all following analysis FASTA files should be given.

^^^^^^^^^^^^^
Metagenomics
^^^^^^^^^^^^^

output: krona.html

* **Analysis 1**

.. code-block:: none
	
	bipype -m run -rapsearch rap_KEGG -humann
	bipype -m run --order prepare_taxonomy_stats -ot txt

* **Analysis 2**

This analysis can be done in the same way as presented in case of analysis of amplicons.

^^^^^^^^^
Amplicons
^^^^^^^^^

* **Both (16S and ITS)**

output: 16S_ITS.krona

.. code-block:: none

	bipype -m run --cutadapt use_filenames ITS 16S
	bipype -m run --order taxonomy_stats -ot ITS 16S

* **Only ITS**

output: ITS.krona 

.. code-block:: none

	bipype -m run --cutadapt use_filenames ITS
	bipype -m run --order taxonomy_stats -ot ITS

* **Only 16S**

output: 16S.krona

.. code-block:: none

	bipype -m run --cutadapt use_filenames 16S
	bipype -m run --order taxonomy_stats -ot  16S

^^^^^^^^^^^^^^^^^^^
Metatranscriptomics
^^^^^^^^^^^^^^^^^^^

---------------------------------------------
More detailed description of bipype workflow
---------------------------------------------

Bipype consists of three major, parts:

* **sample**, which performs most analysis and runs other programs

* **taxonomy_stats**, which generates taxonomy stats in Krona format

* **metatranscriptomics** TODO

^^^^^^^^^^
Sample
^^^^^^^^^^

**General schema of workflow**

1. Samples mapping
2. Alignment reads to reference sequence(s)
3. Determining the presence/absence and abundance of microbial pathways

**Samples mapping**

This step checks to which species and taxonomic group given sample belongs.
It searches databases for fungi, plants or both, accordingly to value of <to_calculate> option:

* p - plant,

* f - fungi,

* b - both

**Taxonomy database loading**

By default the taxonomy database is loaded from path defined by <db_NCBI_taxonomy> option. This should by a “pickle” file and by default it is:


${PATH_NCBI_TAXA_DB}


Databases will be loaded only if mode is “run”.

**What if mentioned, default taxonomy database doesn’t exists?**

Then, there are loaded other, also default databases:

+-------------------------+------------------------+
| Database                | default path           |
+=========================+========================+
| GI → TaxID              | ${PATH_GI_TAX}         |
+-------------------------+------------------------+
| TaxID → scientific_name | ${PATH_TAX_NAME}       |
+-------------------------+------------------------+

Formatting of these databases are described in appendix “databases formatting”.

**How to replace default taxonomy database?**

To load your own database, you need to supersede the <db_NCBI_taxonomy> option with the path to file, where pickled dict with database is saved.

One of possible ways to create this database:

1. Modify databases: “GI → TaxID” and “TaxID → scientific_name”, (discussed below), and then
2. Load modified databases and pickle them to new file (here: ‘new_db’), for example with use of the following python script:

.. code-block:: python
	
	#!/usr/bin/python
	from refseq_bipype import taxa_read
	import cPickle as pickle
	with open('new_db', 'wb') as output_file:
		list_with_dicts = taxa_read("manual")
		pickle.dump(list_with_dicts, output_file)

**Getting input files**

Finds in current working directory and subdirectories fastq files, which are “paired_end”.
If mode is run, also unpacks compressed archives in search of input files.

**Mapping**

A list of path to refseq databases is given in db_refseq_fungi and db_refseq_plant parameters for plant and fungi analysis. Up to two paths are allowed, so for one path program will add “False”, as a second element.

**Reconstruct option:**

If that option was chosen, bwa (Burrows-Wheeler Alignment Tool) program is run. The program finds the SA coordinates of the input reads and generates alignments in the SAM format given single-end reads. Repetitive hits will be randomly chosen. If mode “run” was chosen, commands "bwa aln", "bwa samse" and “samtools mpileup” were run.
Bwa aln command:

+---------+-------------------------------------------------------------------------------------+
| **aln** | bwa aln [-n maxDiff] [-o maxGapO] [-e maxGapE] [-d nDelTail] [-i nIndelEnd] [-k     |
|         | maxSeedDiff] [-l seedLen] [-t nThrds] [-cRN] [-M misMsc] [-O gapOsc] [-E gapEsc] [-q|
|         | trimQual] <in.db.fasta> <in.query.fq> > <out.sai>                                   |
+---------+-------------------------------------------------------------------------------------+

Find the SA coordinates of the input reads. Maximum *maxSeedDiff* differences are allowed in the first *seedLen* subsequence and maximum *maxDiff* differences are allowed in the whole sequence.
Parameter t (number of threads) is taken from reconstruct function (thr parameter), all others will be default.

Bwa samse command:

+-----------+----------------------------------------------------------------------------------+
| **samse** | bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>                 |
+-----------+----------------------------------------------------------------------------------+

Generate alignments in the SAM format given single-end reads. Repetitive hits will be randomly chosen.
All parameters will be default.
Also SAMtools program is run. For example, we have reference sequences in ref.fa, indexed by samtools faidx, and position sorted alignment file aln.bam, the following command lines call SNPs and short INDELs:

**samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq**

**Alignment reads to reference sequence(s)**

After choosing fungi or plant, program runs refseq_mapping function with appropriate parameters for this type of input. Biopype uses Bowtie 2 program to align reads to reference sequence(s). Firstly, Bowtie 2 aligns reads to reference sequence (or sequences, if refseq_2 is not False and merge output SAM files). Secondly, BAM files are made, sorted and indexed. Thirdly, Function launches 'samtools idxstats'. Finally, idx_map() function is called: parse data from samtools idxstats output file and writes data to outfile in key; value format, where:
GI/TaxID/scientific name → number of mapped reads.
For fungies and plants, different basenames of the indexes for the reference genome are used as refseq parameter for refseq_mapping function.

There are values for adapters` list (cutadapt parameter): ‘both’, ‘ITS’, ‘16S’. If list of adapters types isn`t empty, then program searches for the adapters in reads from input files, removes them, when it finds them and writes to output files which have .cutadapt.fastq extension. Then function runs FLASH (software tool to merge paired-end reads) with cutadapt.fastq files as input and .amplicons.cutadapt.flash.merged.fastq files as output. These results are input for fastq_to_fasta which converts file format from fastq to fasta.
If Velvet wasn`t choosen, then rapsearch option is available. RAPSearch will be run with default KEGG=None, if ‘rap_prot’ option is in to_calculate option. RAPSearch will be runned with options:

-z 10 -e 0.001 -b 100 -v 100 -g T -a T 1, if mode is ‘run’.
If ‘rap_KEGG’ or ‘rap_KO’ are in to_calculate list, then all files .fastq will be changed for .fasta.
If rap_KEGG is in to_calculate list, then rapsearch will be runned with KEGG=’masl’.
In other case, rapsearch will be runned with KEGG=’KO’. With this parameter RAPSearch will be runned with options:
-z 12 -v 20 -b 1 -t n -a t 1, if mode is ‘run’.

**Determining the presence/absence and abundance of microbial pathways**

HUMAnN
If ‘humann’ is in to_calculate list, then program checks if m8 files exist. In that case humann function will be runned and new catalog humann-0.99 will be created in rapsearch result folder. This function copies HUMAnN program to the current directory, moves input (\*.m8) files to the input directory, copies hmp_metadata.dat file to the input directory and runs HUMAnN. HUMAnN is a pipeline for efficiently and accurately determining the presence/absence and abundance of microbial pathways in a community from metagenomic data.


In that case, human function will be runned and analysis results will be added in rapsearch result folder.
If ‘16S’ or ‘ITS’ are in to_calculate list (or both of them), then usearch function will be runned for these types of sequences.
Usearch function runned USEARCH command, if mode=’run’:

-usearch_local [katalog z USEARCH] -db [input file] -evalue 0.01 -id 0.9 -blast6out [output file] -strand both -threads [threads (integer)]

^^^^^^^^^^^^^^^^^^^^^^^^
Preparing taxonomy stats
^^^^^^^^^^^^^^^^^^^^^^^^

Performs statistical analysis of taxonomy from appropriate files from current working directory: counts occurrences of different taxa and prepares the results to be presented in HTML format.


Results will be converted to HTML (with the Krona program), but only when <mode> is set to 'run'.


Takes following options: output_type, mode, out_dir.

**Input**

Analysis will be performed on files, which simultaneously:

* are located in current working directory or in subdirectories,

* have suffixes of filenames equal to value of <output_type> option,

* have suffixes of filenames equal to 'usearch_' + value of <output_type> option, but only if <output_type> is 'ITS' or '16S'.

**Output**

Output files will be placed in <out_dir> directory, with basenames depending on <output_type>.
        
Examples:

+-------------------------+------------------------+-----------------------+
| Given <output_type>     | Filename of krona file | Filename of html file |
+=========================+========================+=======================+
| ITS                     | ITS.krona              |  ITS.html             |
+-------------------------+------------------------+-----------------------+
| ITS, 16S                | ITS_16S.krona          |  ITS_16S.html         |
+-------------------------+------------------------+-----------------------+




