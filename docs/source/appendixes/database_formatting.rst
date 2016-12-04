===================
Database formatting
===================


16S taxonomic database
----------------------
| This database should be given as FASTA file.
| Format of headers for 16S should be like in following example: 
|


Example::
   
   >AF093247.1.2007 Eukaryota;Amoebozoa;Mycetozoa;Myxogastria;;Hyperamoeba_sp._ATCC50750


ITS taxonomic database
----------------------
| This database should be given as FASTA file.
| Format of headers for ITS is the same as one used in UNITE database: 
|

Example::

   >DQ233785|uncultured ectomycorrhizal fungus|Fungi|Thelephora terrestris|Fungi; Basidiomycota; Agaricomycotina; Agaricomycetes; Incertae sedis; Thelephorales; Thelephoraceae; Thelephora; Thelephora terrestris


Database “GI → TaxID”
---------------------
| This database should be prepared in tsv (tab-separed values) format.
| First column is a GI, second is a TaxID.
| 

Example::

   13	9913
   15	9915
   16	9771
   17	9771


Database “TaxID → scientific_name”
----------------------------------

From file formated as in example below::

   2	|	prokaryotes	|	prokaryotes <Bacteria>	|	in-part	|
   6	|	Azorhizobium	|		|	scientific name	|
   6	|	Azorhizobium Dreyfus et al. 1988	|		|	synonym	|
   6	|	Azotirhizobium	|		|	equivalent name	|
   7	|	ATCC 43989	|		|	type material	|
   7	|	Azorhizobium caulinodans	|		|	scientific name	|


| Only following data will be extracted:
|

===== ========================
TaxID scientific_name
===== ========================
6     Azotirhizobium
7     Azorhizobium caulinodans
===== ========================
