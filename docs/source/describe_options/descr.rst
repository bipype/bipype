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



