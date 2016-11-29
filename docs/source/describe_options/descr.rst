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


