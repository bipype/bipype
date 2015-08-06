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

def get_tables(baza):
    #-------GET LIST OF TABLES--------------------------------------
    baza.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = baza.fetchall()
    for row in tables:
        print(row[0])

def get_pathways(baza):
    #------Get pathways-------------------------------------------
    #columns (pid, name)
    #EXAMPLE: (ko04060, Cytokine-cytokine receptor interaction)

    baza.execute("select * from Pathways")
    paths = c.fetchall()
    pathways = {}
    for path in paths:
        pathways[path[0]] = path[1]        
    return pathways

def get_kopathways(baza):
    #------Get pathways-------------------------------------------
    #columns (koid, pid)
    #EXAMPLE: (K01194 ko00500)

    baza.execute("select * from KoPathways")
    kopaths = c.fetchall()
    kopathways = {}
    kopath_path = {}
    for kopath in kopaths:
        kopathways[kopath[0]] = kopath[1]
        #kopath = {K12345:ko12345, ...}
        try:
            kopath_path[kopath[1]].add(kopath[0])
        except KeyError:
            kopath_path[kopath[1]] = set([kopath[0]])
        #kopath_path = {ko12345:set([K12345, K12346,...]), ...}
    return kopathways, kopath_path


#---------------sqlite db connection----------------------
conn = sqlite3.connect(db_path)
conn.row_factory = sqlite3.Row
conn.text_factory = str
c = conn.cursor()

#---------------db content acquisition-----------------------
###get_tables(c)
path_names = get_pathways(c)
kopath_keys, kopath_count = get_kopathways(c)

def out_content(filelist, kopath_count, path_names, method='DESeq2'):
    for filepath in filelist:
        if '\\' in filepath:
            filename = filepath.split('\\')[-1]
        elif '/' in filepath:
            filename = filepath.split('/')[-1]
        else:
            filename = filepath
        outname = method+'_'+filename.replace('txt', 'path_counts.csv')
        Kids = set()
        with open(filepath, 'r') as plik:
            filecontent = plik.readlines()[1:]
            for linia in filecontent:
                Kid = linia.rstrip().split('\t')[0]
                Kids.add(Kid)
        with open(outname, 'w') as outplik:
            outplik.write('ko_path_id;ko_path_name;percent common;common KOs\n')
            for path, Kset in kopath_count.items():
                kommon = Kids&Kset
                if len(kommon) > 0:
                    percent_ko = str(int(len(kommon)*100/len(Kset)))
                    print_ko = ' '.join(kommon)
                    path_name_comma = path_names[path]
                    path_name = path_name_comma.replace(',', ' _')
                    outline = ';'.join([path, path_name, percent_ko, print_ko])+'\n'
                    outplik.write(outline)
        
 
edger_files = glob(edger_paths)
deseq_files = glob(deseq_paths)
out_content(deseq_files, kopath_count, path_names)
out_content(edger_files, kopath_count, path_names, 'edgeR')
