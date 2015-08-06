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
multi_id = {}
#pid_d = {}
ko_naming = {}
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
    with open(db_loc, 'rb') as fp:
        imported = pickle.load(fp)
        fp.close()
    return imported

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

def dicto_reduce(present, oversized):
    surplus = set(oversized.keys()) - set(present.keys())
    for gid in surplus:
        del oversized[gid]
    unaccounted = set(present.keys()) - set(oversized.keys()) 
    for gid in unaccounted:
        del present[gid]
    return oversized, present

def m8_to_ko(plik):
    print('working on %s'%(plik))
    tmp_ko_dict = {}
    outname = plik.replace('txt.m8', 'out')
    outname = outname.replace('Sample_GB_RNA_stress_', '')
    content = open(plik, 'r')
    hit_gid = []
    for linia in content:
        gid = linia.split('\t')[1]
        hit_gid.append(gid)
    file_reading_time = time()
    print(plik, 'file_reading seconds', kogenes_time-file_reading_time, 'total time', start_time-file_reading_time)
    gid_count = Counter(hit_gid)
    multi_clean, gid_clean = dicto_reduce(gid_count, multi_id)
    cleaning_time = time()
    print(plik, 'cleaning time seconds', file_reading_time-cleaning_time, 'total time', start_time-cleaning_time)
    todo = len(gid_clean.keys())
    for gid in gid_clean:
        for ko in multi_clean[gid]:
            try:
                tmp_ko_dict[ko] += gid_clean[gid]
            except KeyError:
                tmp_ko_dict[ko] = gid_clean[gid]

    comparison_time = time()
    print(plik, 'comparing time seconds', cleaning_time-comparison_time, 'total time', start_time-comparison_time)
    out_file = open(outname, 'w')
    for ko in tmp_ko_dict:
        to_print = '%s\t%i\n'%(ko, tmp_ko_dict[ko])
        out_file.write(to_print)
    out_file.close()
    writing_time = time()
    print(plik, 'comparing time seconds', comparison_time-writing_time, 'total time', start_time-writing_time)


#m8_list = glob('*m8')
#for plik in m8_list:
#    p = Process(target=m8_to_ko, args=(plik,))
#    p.start()
