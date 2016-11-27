====================
Command line options
====================

General
-------
| Performance and I/O related options.
|
|

Help
~~~~

.. option:: -h, --help

  Show help message and exit.

Adjusting the number of used threads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: -t <threads>, --threads <threads>

  Number of threads to be used. 
  Default: 8

Mode of running
~~~~~~~~~~~~~~~

.. option:: -m <mode>, --mode <mode>

   A mode in which the program will be run.
   Available modes: test, run
   Default: test

Directory with input files
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option::  -r <root_directory>, --root_dir <root_directory>

   Root of the directory tree to be searched.
   Default: ~/bipype

Directory for output files
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption:: -o <output_directory>, --out_dir <output_directory>

   Indicates, where the output files should be located.
   To specify current working directory, use 'in_situ'.
   Otherwise type path to chosen directory.
   Default: in_situ

Insert length
~~~~~~~~~~~~~

.. cmdoption:: --ins_len <insert_length>
  
   Lenght of instert.
   Be advised - better use it for single run.
   Default: 9999

Processed file postfix
~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption::  -postfix <postfix>
   
   Alphanumerical postfix of processed file
   Default: (empty string) 


Continuation of work, after aborting or using manually improved partial results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption:: -e

   Tells bipype, to use existing files with partial results, so parts which are already done (for example by previous, killed instance of program) will be incorporated into pipe. It also creates possibility of manually improving existing files.
   Default: False

Taxonomy stats
--------------
| Options related to prepare_taxonomy_stats results.
|
|

Output type
~~~~~~~~~~~

.. cmdoption:: -ot <output_type_list>, --output_type <output_type_list>

   Simultaneously defines input type!
   Allows to choice on which files (for example ITS, 16S, map_count) the analysis will be performed and also determines basenames of output files. Names of input files should end with (respectively) .usearch_ITS or .usearch_16S, .map_count.
   Default: ['ITS, 16S']

Input cleaning
--------------
| Methods for cleaning input from noise.
|
|

Initial cleaning
~~~~~~~~~~~~~~~~

.. cmdoption:: -ic <method>, --initial_cleaning <method>

   Initial cleaning method
   Available methods: usearch, fastx
   Default: (empty string)

Adapters cutting
~~~~~~~~~~~~~~~~

.. cmdoption:: --cutadapt <adapter_file> <search_options>

   Location of file with adapters to be used by cutadapt (possible use of "use_filenames" to determine adapters from hardcode), and list of usearches to be run on created files - possible options are 16S, ITS, both. Please note, that mapping options -16S, -ITS are completely irrelevant if you use cutadapt. Other note - this is !!!IMPORTANT!!! to present location of file with adapters as first option of this argument
   Default: (empty string)

Mappings
--------
| Mappings to be done during the run.
|
|

ITS usearch
~~~~~~~~~~~

.. cmdoption:: -ITS

   Perform usearch on ITS database

16S usearch
~~~~~~~~~~~

.. cmdoption:: -16S

   Perform usearch on 16S database

RefSeq
~~~~~~

.. cmdoption::  -refseq <kingdom>

   Map samples on refseq.
   Available kingdoms: p, f, b
   p states for plantae, f for fungi, b - both
   Default: f


Contig reconstruction
---------------------
| Methods of contig reconstruction to be used during the run.
|
|

MetaVelvet
~~~~~~~~~~

.. cmdoption::  -MV <parameters> 	

   Parameters for MetaVelvet.
   Parameters format:
   [initial_k_mer_size, final_k-mer_size, step]

Reconstruct
~~~~~~~~~~~

.. cmdoption:: -reconstruct <option>

   Reconstruct relating to database. 
   Available options: database_loc, prefix

Humann
~~~~~~

.. cmdoption:: -humann

   Mapping rapsearch using humann (default: None)

Rapsearch
~~~~~~~~~

.. cmdoption:: -rapsearch <database>

   Perform RAPsearch on selected protein database
   Available databases: masl28,rap_prot,rap_KO
   Default: rap_prot

Databases
---------
| Locations of databases.
|
|

16S for usearch
~~~~~~~~~~~~~~~

.. cmdoption:: --db_16S <database>

   16S database to use in usearch (bowtie indexed)
   Default: ${PATH_X16S_DB}


ITS for usearch
~~~~~~~~~~~~~~~

.. cmdoption:: --db_ITS <database>

   ITS database to use in usearch (bowtie indexed)
   Default: ${PATH_ITS_DB}

Fungi for refseq
~~~~~~~~~~~~~~~~

.. cmdoption:: --db_refseq_fungi <databases>

   Refseq database to use for fungi analysis. All files found under your_database_path/*.some_suffix path will be loaded and treated as subparts of your database. 
   Up to two paths are allowed (paths should be separated with space).
   Default: ${PREF_PATH_REF_FUNGI}

Plant for refseq
~~~~~~~~~~~~~~~~

.. cmdoption:: --db_refseq_plant

   Refseq database to use for plants analysis. All files found under your_database_path/*.some_suffix path will be loaded and treated as subparts of your database. 
   Up to two paths are allowed.
   Default: ${PREF_PATH_REF_PLANT_1} ${PREF_PATH_REF_PLANT_2}

Taxonomy database
~~~~~~~~~~~~~~~~~

.. cmdoption:: --db_NCBI_taxonomy

   Location of cPickled database with mappings of NCBI_tax_id, NCBI tax names and NCBI tax ids.
   Default: ${PATH_NCBI_TAXA_DB}


Database for reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption:: --db_reconstruct

   This database will be passed as parameter to bwa program.
   Format: Fasta file.
   Default: ${PATH_RECONSTRUCT_DB}

Database for ITS taxonomy statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption:: --db_taxonomy_16S

   This database is used to create a taxonomy tree for results visualisations.
   By default it is specially formatted FASTA file.
   For more informations check “databases formatting” chapter.
   Default: ${PATH_16S_DATABASE}

Database for ITS taxonomy statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmdoption:: --db_taxonomy_ITS

   This database is used to create a taxonomy tree for results visualisations.
   By default it is specially formatted FASTA file (with headers like in UNITE database).
   For more informations check “databases formatting” chapter.
   Default: ${PATH_ITS_DATABASE}
