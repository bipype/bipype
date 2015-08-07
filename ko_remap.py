#!/usr/bin/python
import sqlite3
from pprint import pprint
from glob import glob
from multiprocessing import Process
from os.path import exists as pexists
from collections import Counter

db_path = 'ko.db'
edger_paths = 'tables_edgeR/*txt'
deseq_paths = 'tables_DESeq2/*txt'

def get_tables(database):
    """Prints all tables included in SQLite database.
                    
    Arg:
        database: Cursor object to SQLite database.    
    """
    database.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = database.fetchall()
    for row in tables:
        print(row[0])

def get_pathways(database):
    """Make dictionary from pathways table from SQLite database.
    
    Arg:
        database: Cursor object to SQLite database. 

    Returns:
        Dictionary in following format:
                {PID:Name}
            For example:
                {   'ko04060':'Cytokine-cytokine receptor interaction',
                    'ko00910':'Nitrogen metabolism'    }
    """
    database.execute("select * from Pathways")
    paths = c.fetchall()
    pathways = {}
    for path in paths:
        pathways[path[0]] = path[1]        
    return pathways

def get_kopathways(database):
    """Makes dictionaries from kopathways table from SQLite database.
    
    Arg:
        database: Cursor object to SQLite database. 

    Returns:
        Two dictionaries:
            {KO identifier:PID}
            For example:
                {   'K01194':'ko00500',
                    'K04501':'ko04390'    }
            &
            {PID:set[KO identifiers]}
            For example:
                {ko12345:set([K12345, K12346,...]),...}
    """
    database.execute("select * from KoPathways")
    kopaths = c.fetchall()
    kopathways = {}
    kopath_path = {}
    for kopath in kopaths:
        kopathways[kopath[0]] = kopath[1]
        try:
            kopath_path[kopath[1]].add(kopath[0])
        except KeyError:
            kopath_path[kopath[1]] = set([kopath[0]])
    return kopathways, kopath_path


#---------------SQLite database connection----------------------------------
conn = sqlite3.connect(db_path)
conn.row_factory = sqlite3.Row
conn.text_factory = str
c = conn.cursor()


#---------------Database content acquisition--------------------------------
###get_tables(c)
path_names = get_pathways(c)
kopath_keys, kopath_count = get_kopathways(c)


def out_content(filelist, kopath_count, path_names, method='DESeq2'):
    """For every item in 'kopath_count' dictionary and for every file 
    in 'filelist', writes to output file line with KOs, which are common
    for item.value and the set of KOs obtained from file.   
    
    Args:
        filelist:      List of paths to tab-delimited .txt files, where
                       first column is a KO identifier.
                       
        kopath_count:  Dictionary in {PID:set[KO identifiers]} format.
                       For example:
                            {ko12345:set([K12345, K12346,...]),...}
                       
        path_names:    Dictionary in {PID:Name} format.     
                       For example:
                   {'ko04060':'Cytokine-cytokine receptor interaction',
                    'ko00910':'Nitrogen metabolism'}
                    
        method:        Argument used only as a part of output file name
                       (default: 'DESeq2')
   
    Output file has following name:
            (method+'_'+filename.replace('txt', 'path_counts.csv'))    
        where:
            filename = filepath.split('\\')[-1], if '\\' in filepath.
            filename = filepath.split('/')[-1],  if '/' in filepath.
            filename = filepath,                 in other cases.
        
    & following headline (format):
        ko_path_id;ko_path_name;percent common;common KOs
        
    Writes only lines with non-zero common KOs.
    """
    for filepath in filelist:
        if '\\' in filepath:
            filename = filepath.split('\\')[-1]
        elif '/' in filepath:
            filename = filepath.split('/')[-1]
        else:
            filename = filepath
        outname = method+'_'+filename.replace('txt', 'path_counts.csv')
        Kids = set()
        with open(filepath, 'r') as file_:
            filecontent = file_.readlines()[1:]
            for line in filecontent:
                Kid = line.rstrip().split('\t')[0]
                Kids.add(Kid)
        with open(outname, 'w') as outfile:
            outfile.write('ko_path_id;ko_path_name;percent common;common KOs\n')
            for path, Kset in kopath_count.items(): # {PID:set[KO identifiers]}
                common = Kids&Kset
                if len(common) > 0:
                    percent_ko = str(int(len(common)*100/len(Kset)))
                    print_ko = ' '.join(common)
                    path_name_comma = path_names[path]
                    path_name = path_name_comma.replace(',', ' _')
                    outline = ';'.join([path, path_name, percent_ko, print_ko])+'\n'
                    outfile.write(outline)
        
 
edger_files = glob(edger_paths)
deseq_files = glob(deseq_paths)
out_content(deseq_files, kopath_count, path_names)
out_content(edger_files, kopath_count, path_names, 'edgeR')
