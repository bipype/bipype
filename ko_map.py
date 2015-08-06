#!/usr/bin/python
import sqlite3
from pprint import pprint
from glob import glob
from multiprocessing import Process
import cPickle as pickle 
from os.path import exists as pexists
from collections import Counter
from time import time

start_time = time()

conn = sqlite3.connect('ko.db')
conn.row_factory = sqlite3.Row
conn.text_factory = str
c = conn.cursor()
#ko_naming = {}
#pid_d = {}
#Path_db = c.execute("select * from Pathways")
#Path_db_all = c.fetchall()
#for pid, Path_name in Path_db_all:
#    pid_d[pid] = Path_name
#KoPath_db = c.execute("select * from KoPathways")
#KoPath_db_all = c.fetchall()
#for koid, pid in KoPath_db_all:
#    if koid not in ko_naming:
#        ko_naming[koid] = {pid:pid_d[pid]}
#    else:
#        ko_naming[koid][pid] = pid_d[pid]

def auto_tax_read(db_loc):
    """Reads pickled {KEGG GENES number : set[KO identifiers]} dictionary."""    
    with open(db_loc, 'rb') as fp:
        imported = pickle.load(fp)
        fp.close()
    return imported
    
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# This bunch of code cheks if appropriate pickle is available 
# (HARDODED: 'kogenes.pckl') and loads it to multi_id dictionary. If there
# is no such file, 'kogenes.pckl' is made from SQL database (Variable 'c').
# Pickled: {GI identifier : set[KO identifiers]} dictionary.
multi_id = {}
if pexists('kogenes.pckl'):
    multi_id = auto_tax_read('kogenes.pckl')
    kogenes_time = time()
    print('kogenes reading time', start_time-kogenes_time)
else:
    KoPath_gid = c.execute("select * from KoGenes")
    KoPath_gid_all = c.fetchall()
    for koid, gid in KoPath_gid_all:
        if gid not in multi_id:
            multi_id[gid] = set([koid])
        else:
            multi_id[gid].add(koid)
    output = open('kogenes.pckl', 'w')
    pickle.dump(multi_id, output)
    output.close()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

def dicto_reduce(present, oversized):
    """Removes all elements from dictionaries, which keys aren't present
    in both.
    
    Args:
        present:                Dictionary
        oversized:              Dictionary
    
    Results:
        oversized, present      Dictionaries      
        
    ATTENTION: Order of parametres is opposite to results. 
    
    Example:
        >>> dict_1={'a':1,'c':3,'d':4}
        >>> dict_2={'a':3,'b':4,'c':4}
        >>> dicto_reduce(dict_1, dict_2)
        ({'a': 3, 'c': 4}, {'a': 1, 'c': 3})
        >>>
    """
    surplus = set(oversized.keys()) - set(present.keys())
    for gid in surplus:
        del oversized[gid]
    unaccounted = set(present.keys()) - set(oversized.keys()) 
    for gid in unaccounted:
        del present[gid]
    return oversized, present

def m8_to_ko(file_): 
    # multi_id should be a function argument
    # Computational complexity is really bad
    """Assigns and counts GI identifiers from BLAST Tabular (flag: -m 8)
    output format file, for every KO from multi_id variable.

    After mapping, writes data to output file.
    
    Args:
        file_: Path to BLAST Tabular (flag: -m 8) output format file.
        
        VARIABLE !!! 
        multi_id format: {GI identifier : set[KO identifiers]}
        VARIABLE !!!   
    
    Output file (outname) has following path:
        outname = file_.replace('txt.m8', 'out')
        outname = outname.replace('Sample_GB_RNA_stress_', '')    
    & following format:
        K00161  2
        K00627  0
        K00382  11
    """
    #print('working on %s'%(file_))
    tmp_ko_dict = {}
    outname = file_.replace('txt.m8', 'out')
    outname = outname.replace('Sample_GB_RNA_stress_', '')
    content = open(file_, 'r')
    hit_gid = [] # List of GI identifiers from file_
    for linia in content:
        gid = linia.split('\t')[1]
        hit_gid.append(gid)
    file_reading_time = time()
    print(file_, 'file_reading seconds', kogenes_time-file_reading_time, 'total time', start_time-file_reading_time)
    gid_count = Counter(hit_gid)
    multi_clean, gid_clean = dicto_reduce(gid_count, multi_id)
    #multi_id format: {GI identifier : set[KO identifiers]}
    cleaning_time = time()
    print(file_, 'cleaning time seconds', file_reading_time-cleaning_time, 'total time', start_time-cleaning_time)
    todo = len(gid_clean.keys()) # not used
    for gid in gid_clean:
        for ko in multi_clean[gid]:
            try:
                tmp_ko_dict[ko] += gid_clean[gid]
            except KeyError:
                tmp_ko_dict[ko] = gid_clean[gid]

    comparison_time = time()
    print(file_, 'comparing time seconds', cleaning_time-comparison_time, 'total time', start_time-comparison_time)
    out_file = open(outname, 'w')
    for ko in tmp_ko_dict:
        to_print = '%s\t%i\n'%(ko, tmp_ko_dict[ko])
        out_file.write(to_print)
    out_file.close()
    writing_time = time()
    print(file_, 'comparing time seconds', comparison_time-writing_time, 'total time', start_time-writing_time)


#m8_list = glob('*m8')
#for plik in m8_list:
#    p = Process(target=m8_to_ko, args=(plik,))
#    p.start()
