import sqlite3
import cPickle
import subprocess
from bs4 import BeautifulSoup
from glob import glob
from urllib import urlopen
from multiprocessing import Process
from collections import Counter
from time import time
from datetime import datetime
from os.path import exists as pexists
from os.path import join as pjoin
from os.path import dirname, realpath
from os import system, chdir, getcwd
from settings_bipype import *


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
    return (oversized, present)


def connect_db(db):
    """Connects database

        Arg:
            db:     Path to SQL database

        Returns:
            Cursor object to database
    """
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    return conn.cursor()


def get_tables(database):
    """Prints all tables included in SQLite3 database.

        Arg:
            database: Cursor object to SQLite3 database.
    """
    database.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = database.fetchall()
    for row in tables:
        print row[0]


def auto_tax_read(db_loc):
    """Reads pickled {KEGG GENES number : set[KO identifiers]} dict."""
    with open(db_loc, 'rb') as file_:
        dictionary = cPickle.load(file_)
        file_.close()
    return dictionary


def pickle_or_db(pickle, db):
    """Reads pickle or SQL database, than makes a dict.

        If appropriate pickle (a dict) is available, it is read.
        In the other case function reads 'kogenes' table from
        SQL database and makes missing pickle. Eventually returns dict.


        Args:
            pickle: Path to pickled dict in following format:
                    {KEGG GENES identifier : set[KO identifiers]}

            db:     Cursor object to SQL database with 'kogenes' table
                    (KO identifier          KEGG GENES identifier)
        Returns:
            Dict in {KEGG GENES identifier : set[KO identifiers]} format.


        Some information for Bipype's developers
        (delete this before final version):
        Code from this fuction was not a fuction in previous version and
        'args' was hardcoded to:
                        'kogenes.pckl' & c (variable with db's cursor)
    """
    start_time = time()
    multi_id = {}
    if pexists(pickle):
        multi_id = auto_tax_read(pickle)
        kogenes_time = time()
        print ('kogenes reading time', kogenes_time - start_time)
    else:
        db.execute('select * from KoGenes')
        KoPath_gid_all = db.fetchall()
        for (koid, gid,) in KoPath_gid_all:
            if gid not in multi_id:
                multi_id[gid] = set([koid])
            else:
                multi_id[gid].add(koid)
        with open(pickle, 'w') as output:
            cPickle.dump(multi_id, output)
    return multi_id


def get_pathways(database):
    """Make dictionary from pathways table from SQLite3 database.

        Arg:
            database: Cursor object to SQLite3 database.

        Returns:
            Dictionary in following format:
                    {KEGG_Pathway_id:Name}
                For example:
                    {   'ko04060':'Cytokine-cytokine receptor interaction',
                        'ko00910':'Nitrogen metabolism'    }
    """
    database.execute('select * from Pathways')
    paths = database.fetchall()
    pathways = {}
    for path in paths:
        pathways[path[0]] = path[1]
    return pathways


def get_kopathways(database):
    """Makes dictionaries from kopathways table from SQLite3 database.

        Arg:
            database: Cursor object to SQLite3 database.

        Returns:
            Two dictionaries:
                {KO identifier:set[KEGG_Pathway_ids]}
                For example:
                    {   'K01194':set(['ko00500','ko00600',...]),
                        'K04501':set(['ko04390',...])   }
                &
                {KEGG_Pathway_id:set[KO identifiers]}
                For example:
                    {ko12345:set([K12345, K12346,...]),...}
    """
    database.execute('select * from KoPathways')
    kopaths = database.fetchall()
    kopathways = {}
    kopath_path = {}
    for kopath in kopaths:
        if kopath[0] not in kopathways:
            kopathways[kopath[0]] = set([kopath[1]])
        else:
            kopathways[kopath[0]].add(kopath[1])
        try:
            kopath_path[kopath[1]].add(kopath[0])
        except KeyError:
            kopath_path[kopath[1]] = set([kopath[0]])
    return (kopathways, kopath_path)


def m8_to_ko(file_, multi_id):
    """Assigns and counts KEGG GENES identifiers from BLAST Tabular
        (flag: -m 8) output format file, for every KO from multi_id.

        After mapping, writes data to output file.

        Args:
            file_:    Path to BLAST Tabular (flag: -m 8) format file
            multi_id: Dict {KEGG GENES identifier : set[KO identifiers]}

        Output file (outname) has following name:
            outname = file_.replace('txt.m8', 'count')
        & following format:
            K00161  2
            K00627  0
            K00382  11
    """
    start_time = time()
    tmp_ko_dict = {}
    outname = file_.replace('txt.m8', 'count')
    content = open(file_, 'r')
    hit_gid = []
    for line in content:
        if line[0] != '#':
            gid = line.split('\t')[1]
            hit_gid.append(gid)
    file_reading_time = time()
    print (file_, 'file_reading seconds', file_reading_time - start_time)
    gid_count = Counter(hit_gid)
    (multi_clean, gid_clean,) = dicto_reduce(gid_count, multi_id)
    cleaning_time = time()
    print (file_, 'cleaning time seconds', cleaning_time - file_reading_time)
    for gid in gid_clean:
        for ko in multi_clean[gid]:
            try:
                tmp_ko_dict[ko] += gid_clean[gid]
            except KeyError:
                tmp_ko_dict[ko] = gid_clean[gid]
    comparison_time = time()
    print (file_, 'comparing time seconds', comparison_time-cleaning_time)
    with open(outname.replace('m8', 'counts'), 'w') as out_file:
        for ko in tmp_ko_dict:
            to_print = '%s\t%i\n' % (ko, tmp_ko_dict[ko])
            out_file.write(to_print)
    writing_time = time()
    print (file_, 'comparing time seconds', writing_time - comparison_time,
        'total time', writing_time - start_time)


def out_content(filelist, kopath_values, path_names, method = 'DESeq2'):
    """For every item in 'kopath_values' dictionary and for every file
        in 'filelist', writes to output file line with KOs, which are common
        for item.value and the set of KOs obtained from file.

        Args:
            filelist:      List of paths to tab-delimited .txt files, where
                           first column is a KO identifier.

            kopath_values: {KEGG_Pathway_id:set[KO identifiers]} dict.
                           For example:
                                {ko12345:set([K12345, K12346,...]),...}

            path_names:    Dictionary in {KEGG_Pathway_id:Name} format.
                           For example:
                       {'ko04060':'Cytokine-cytokine receptor interaction',
                        'ko00910':'Nitrogen metabolism'}

            method:        Argument used only as a part of output file name
                           (default: 'DESeq2')

        Output file has following name:
                (method+'_'+filename.replace('txt', 'path_counts.csv'))
            where:
                filename = filepath.split('\')[-1], if '\' in filepath.
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
        outname = method + '_' + filename.replace('txt', 'path_counts.csv')
        Kids = set()
        with open(filepath, 'r') as file_:
            filecontent = file_.readlines()[1:]
            for line in filecontent:
                Kid = line.rstrip().split('\t')[0]
                Kids.add(Kid)
        with open('ko_remap/' + outname, 'w') as outfile:
            outfile.write('ko_path_id;ko_path_name;percent common;common KOs\n')
            for (path, Kset,) in kopath_values.items():
                common = Kids & Kset
                if len(common) > 0:
                    percent_ko = str(len(common) * 100.0 / len(Kset))
                    print_ko = ' '.join(common)
                    path_name_comma = path_names[path]
                    path_name = path_name_comma.replace(',', ' _')
                    outline = ';'.join([path,
                     path_name,
                     percent_ko,
                     print_ko]) + '\n'
                    outfile.write(outline)


def fastq_to_fasta(fastq):
    """Runs fastq_to_fasta for _file and write output as out_file.

        Writes output in "fasta/" directory.

        GLOBAL:
            - path to fastq_to_fasta program:               PATH_FQ2FA
    """
    out_file = fastq.rsplit('.', 1)[0] + '.fasta'
    out_file = out_file.rsplit('/', 1)[-1]
    out_file = 'fasta/' + out_file
    subprocess.check_call([PATH_FQ2FA,
     '-i', fastq, '-o', out_file, '-Q33'])


def rapsearch2(input_file, threads):
    """Runs rapsearch2() for input_file in fasta format.

        Writes outputs in "m8/" directory.

        GLOBALS:
            - path to RAPSearch2 program:                   PATH_RAPSEARCH
            - path to similarity search database:           PATH_REF_PROT_KO
    """
    out_name = input_file.replace('tmp.fasta', 'txt')
    out_name = out_name.replace('fasta/', 'm8/')
    subprocess.check_call([PATH_RAPSEARCH,'-q', input_file, '-d',
     PATH_REF_PROT_KO, '-o', out_name, '-z', str(threads),
     '-v', '20', '-b', '1', '-t', 'n', '-a', 't'])


def get_ko_fc(ko_dict, ref_cond, filepath):
    """From given table file (SARTool), adds found fold changes to ko_dict.

        Args:
            ko_dict:    {KO_id:{cond1:value1, cond2:value2...}...} dict
            ref_cond:   reference condition (string)
            filepath:   filepath to output table file from edgeR or DESeq2

        Returns:
            ko_dict with added fold changes from table file
        """
    with open(filepath) as _file:
        l_cond = filepath.split('/')[-1].split('.')[0].split('vs', 2)
        if l_cond[0] != ref_cond:
            cond = l_cond[0]
        else:
            cond = l_cond[1]
        file_content = _file.readlines()
        fc_index = file_content[0].rstrip().split('\t').index('FC')
        for i in range(1, len(file_content)):
            line = file_content[i].rstrip().split('\t')
            KO = line[0]
            base_fc = line[fc_index]
            fc_num = True
            try:
                fc = float(base_fc)
            except ValueError:
                fc = base_fc
                fc_num = False
            if '.down.' in filepath and fc_num:
                fc = -1.0 / linia[fc_index]
            if KO in ko_dict.keys():
                ko_dict[KO][cond] = fc
            else:
                ko_dict[KO] = {cond: fc}
    return ko_dict


def low_change(ko_dict, all_conds):
    """For every KO adds condition:0, where condition is missing.

        Args:
            ko_dict: {KO_id:{cond1:value1, cond2:value2...}...} dict
            all_conds: list of conditions (list of strings)

        Returns:
            suplemented ko_dict

        For example:
            low_change( { 'K12345':{'pH5':1.41, 'pH6':1.73},
                          'K23456':{'pH6':2.0, 'pH8':2.24}   },
                        ['pH5','pH6','pH8'] )
            give        { 'K12345':{'pH5':1.41, 'pH6':1.73, 'pH8':0.0},
                          'K23456':{'pH5':0.0, 'pH6':2.0, 'pH8':2.24}   }

    """
    for (KO, conds,) in ko_dict.items():
        for cond in all_conds:
            if cond not in conds:
                ko_dict[KO][cond] = 0.0
    return ko_dict


def get_kegg_name(ko):
    """Returns name assigned to given KO identifier (from kegg.jp)

        Arg:
            ko: KO identifier (string)

        Returns:
            name assigned to ko (string)
    """
    url = 'http://www.kegg.jp/dbget-bin/www_bget?ko:' + ko
    name = 'NA'
    ko_kegg = urlopen(url).read().decode('utf-8')
    soup = BeautifulSoup(ko_kegg, 'html.parser')
    if soup:
        table = soup.find('td', attrs={'class': 'fr4'})
        if table:
            rows = table.findAll('tr')
            if rows:
                for tr in rows:
                    th = tr.find('th')
                    if th:
                        colname = th.find('nobr')
                        if colname:
                            if 'Name' in colname:
                                content = tr.find('td')
                                if content:
                                    name = content.find(text=True)
    return name


def mapper(ko_dict, ko_set):
    """Assings every KO_id from ko_dict to KEGG_Pathway_id from ko_set

        Args:
            ko_dict:    {KO_id:{cond1:value1, cond2:value2...}...} dict
            ko_set:     {KEGG_Pathway_id:set[KO identifiers]} dict

        Returns:
            mapper_d:
              {KEGG_Pathway_id:{KO_id:{cond1:value1, cond2:value2...}...}...}
    """
    all_ko = set(ko_dict.keys())
    mapper_d = {}
    for (path_name, path_Kids,) in ko_set.items():
        kommon = all_ko & path_Kids
        if len(kommon) > 0:
            mapper_d[path_name] = {}
            for Kid in kommon:
                fc = ko_dict[Kid]
                mapper_d[path_name][Kid] = fc
    return mapper_d


def mapper_write(ko_path_dict, all_conds, out_dir):
    """Writes file with KO and corresponding fold change,
        for every combination of condition & KEGG_Pathway_id .

        Args:
            ko_path_dict:
              {KEGG_Pathway_id:{KO_id:{cond1:value1, cond2:value2...}...}...}
            all_conds:  list of conditions (list of strings)
            out_dir:    relative output directory path

        Output file has following path:
            out_dir/condX/
                      , following name:
            KEGG_Pathway_id.txt
                      , following header:
            # KO KEGG_Pathway_id
                      & following format:
            KO_id corresponding_fold_change
    """
    for cond in all_conds:
        path = pjoin(out_dir, cond)
        system('mkdir ' + path)
        for pathway in ko_path_dict.keys():
            header = '#KO ' + path + '\n'
            filename = pathway + '.txt'
            output_path = pjoin(path, filename)
            with open(output_path, 'w') as _file:
                _file.write(header)
                for Kid in ko_path_dict[pathway]:
                    line = '\t'.join([Kid, str(ko_path_dict[pathway][Kid][cond]), '\n'])
                    _file.write(line)


def config_from_file(work_dir, _file):
    """Reads paramaters from configuration _file. Prepares target.txt

        Args:
            work_dir:   current working directory
            _file       configuration file for metatranscriptomic pipeline

        Returns:
            ref_cond:   reference condition defined by user
            all_conds:  set of conditions (groups) from target.txt
            fastqs:     list of fastq files on which analysis will be done
    """
    all_conds = []
    fastqs = []
    with open(_file) as f:
        lines = f.readlines()
        ref_cond = lines[0].split()[0]
    with open('target.txt', 'w') as f:
        f.write('label\tfiles\tgroup')
        idents = []
        for line in lines[1:]:
            id_num = 1
            line = line.split()
            all_conds.append(line[2])
            while True:
                hyp_ident = line[2] + '_' + str(id_num)
                if hyp_ident in idents:
                    id_num+=1
                else:
                    idents.append(hyp_ident)
                    break
            fastqs.append(line[0])
            fastqs.append(line[1])
            target_name = line[0].rsplit('/', 1)[-1]
            target_name = target_name.replace('R1_', '')
            target_name = target_name.replace('.fastq', '.count')
            f.write('\n' + hyp_ident '_' + '\t' + target_name + '\t' + line[2])
    with open('template_script_DESeq2.r') as f:
        lines = f.readlines()
        lines[24] = 'condRef <- "' + ref_cond + '"' + '\n'
    with open('template_script_DESeq2.r', 'w') as f:
        for line in lines:
            f.write(line)
    with open('template_script_edgeR.r') as f:
        lines = f.readlines()
        lines[24] = 'condRef <- "' + ref_cond + '"' + '\n'
    with open('template_script_edgeR.r', 'w') as f:
        for line in lines:
            f.write(line)
    return (ref_cond, set(all_conds), fastqs)


def run_fastq_to_fasta(fastqs):
    """Runs fastq_to_fasta() for every .fastq in fastqs"""
    for fastq in fastqs:
        fastq_to_fasta(fastq)


def run_cat_pairing():
    """Merges fasta files with paired-end reads in fasta/"""
    for file_R1 in glob('fasta/*R1*fasta'):
        for file_R2 in glob('fasta/*R2*fasta'):
            if file_R1.split('R1') == file_R2.split('R2'):
                outname = file_R1.replace('R1_', '')
                outname = outname.replace('.fasta', '.tmp.fasta')
                subprocess.check_call(['touch', outname])
                with open(outname, 'w') as out:
                    with open(file_R1) as f1:
                        with open(file_R2) as f2:
                            out.write(f1.read())
                            out.write(f2.read())


def run_rapsearch(threads):
    """Runs rapsearch2() for every .tmp.fasta in fasta/"""
    for _file in glob('fasta/*tmp.fasta'):
        rapsearch2(_file, threads)


def run_ko_map():
    """Runs m8_to_ko() for every .m8 file in raw directory.

        GLOBALS:
            - path to KO database:                                  PATH_KO_DB
            - pickle to dict from KO GENES table from KO database:  PATH_KO_PCKL
    """
    data = pickle_or_db(PATH_KO_PCKL, connect_db(PATH_KO_DB))
    p_list = []
    for file_ in glob('m8/*m8'):
        p = Process(target=m8_to_ko, args=(file_, data))
        p.start()
        p_list.append(p)
    for p in p_list:
        p.join()


def run_SARTools():
    """Runs SARTools in R.

        HARDCODED: R templates:
                        edger: template_script_DESeq2.r'
                        deseq: template_script_edgeR.r'
    """
    system('Rscript template_script_DESeq2.r')
    system('Rscript template_script_edgeR.r')


def run_pre_ko_remap(ref_cond):
    """Prepares args for run_(pre_/new_)ko_remap()

        Arg:
            ref_cond:       Reference condition (group) - string

        Returns:
            path_names:     {KEGG_Pathway_id:Name} dict
            kopath_keys:    {KO identifier:set[KEGG_Pathway_ids]} dict
            kopath_values:  {KEGG_Pathway_id:set[KO identifiers]} dict
            edger_files:    list of edgeR outputs paths
            deseq_diles:    list of DESeq outputs paths

        HARDCODED: Paths to files from SARTools:
                        edger: 'edger/*[pn].txt'
                        deseq: 'deseq/*[pn].txt'
        GLOBALS:
            - path to KO database:  PATH_KO_DB
    """
    cursor = connect_db(PATH_KO_DB)
    path_names = get_pathways(cursor)
    (kopath_keys, kopath_values,) = get_kopathways(cursor)
    edger_files = glob('edger/tables/*.txt')
    #edger_files = glob('edger/tables/*[pn].txt') change it after tests!!!!!
    deseq_files = glob('deseq/tables/*.txt')
    #deseq_files = glob('deseq/tables/*[pn].txt') change it after tests!!!!!
    return (path_names, kopath_keys, kopath_values, edger_files, deseq_files)


def run_ko_remap(deseq_files, edger_files, kopath_values, path_names):
    """Runs out_content(files, kopath_values, path_names (,'edgeR'))
        for files from 'edger_paths' & 'deseq_paths'.

        Args:
            deseq_diles:    list of DESeq outputs paths
            edger_files:    list of edgeR outputs paths
            kopath_values:  {KEGG_Pathway_id:set[KO identifiers]} dict
            path_names:     {KEGG_Pathway_id:Name} dict
    """
    out_content(deseq_files, kopath_values, path_names)
    out_content(edger_files, kopath_values, path_names, 'edgeR')


def run_new_ko_remap(deseq_files, edger_files, kopath_values, all_conds, ref_cond):
    """Runs get_ko_fc(), low_change(), mapper() and mapper_write()
        in appropriate way for files from deseq_files and edger_files.

        Args:
            deseq_diles:    list of DESeq outputs paths
            edger_files:    list of edgeR outputs paths
            ref_cond:       Reference condition (group) - string
            kopath_values:  {KEGG_Pathway_id:set[KO identifiers]} dict
            all_conds:      list of conditions (list of strings)

        Returns:
            ko_dict_deseq:  {KO_id:{cond1:value1, cond2:value2...}...} dict
            ko_dict_edger:  {KO_id:{cond1:value1, cond2:value2...}...} dict

        HARDCODED: Output directories paths:
                    deseq: 'new_ko_remap/deseq/'
                    edger: 'new_ko_remap/edger/'
    """
    ko_dict_deseq = {}
    ko_dict_edger = {}
    for _file in deseq_files:
        if ref_cond in _file:
            ko_dict_deseq = get_ko_fc(ko_dict_deseq, ref_cond, _file)
    for _file in edger_files:
        if ref_cond in _file:
            ko_dict_edger = get_ko_fc(ko_dict_edger, ref_cond, _file)
    ko_dict_deseq = low_change(ko_dict_deseq, all_conds)
    ko_dict_edger = low_change(ko_dict_edger, all_conds)
    mapper_deseq = mapper(ko_dict_deseq, kopath_values)
    mapper_edger = mapper(ko_dict_edger, kopath_values)
    mapper_write(mapper_deseq, all_conds, 'new_ko_remap/deseq/')
    mapper_write(mapper_edger, all_conds, 'new_ko_remap/edger/')
    return (ko_dict_deseq, ko_dict_edger)


def run_ko_csv(ko_dict_deseq, ko_dict_edger, all_conds, kopath_keys, path_names, ref_cond):
    """For given ko_dicts writes CSV files with pathways and foldchanges

        Args:
            ko_dict:        {KO_id:{cond1:value1, cond2:value2...}...} dict
            all_conds:      list of conditions (list of strings)
            kopath_keys:    {KO identifier:set[KEGG_Pathway_ids]} dict
            path_names:     {KEGG_Pathway_id:Name} dict
            filepath:       output filepath

        Output files have following format (and header):
            KO_id;Gene_name;paths ids;paths names;FC vs cond1;FC vs cond2;...;

        HARDCODED: Output files paths:
                deseq: 'csv/deseq.csv'
                edger: 'csv/edger.csv'
    """
    for touple in [(ko_dict_deseq, 'csv/deseq.csv'),
                   (ko_dict_edger, 'csv/edger.csv')]:
        filepath = touple[1]
        ko_dict = touple[0]
        with open(filepath, 'wb') as _file:
            l_header = ['KO_id', 'Gene_name', 'paths ids', 'paths names']
            for cond in all_conds:
                l_header.append('FC ' + ref_cond + ' vs ' + cond)
            header = ';'.join(l_header + ['\n'])
            _file.write(header)
            for Kid in ko_dict.keys():
                to_write = [Kid]                    # ['K01369']
                to_write.append(get_kegg_name(Kid)) # ['K01369','LGMN']
                if Kid in kopath_keys.keys():
                    koids = []
                    konames = []
                    for path in kopath_keys[Kid]:
                        koids.append(path)
                        konames.append(path_names[path])
                    to_write.append(','.join(koids))    # ['K01369','LGMN', 'ko04612, ko04142']
                    to_write.append(','.join(konames))  # ['K01369','LGMN', 'ko04612, ko04142','Antigen processing and presentation,Lysosome']
                else:
                    to_write.append('NA')
                    to_write.append('NA')
                for cond in all_conds:
                # ['K01369','LGMN', 'ko04612, ko04142','Antigen processing and presentation,Lysosome','0','1.12','2.32']
                    to_write.append(str(ko_dict[Kid][cond]))
                _file.write(';'.join(to_write) + '\n')


def metatranscriptomics(opts):
    """Performs analyse of metagenomic data.

        For more information please refer to
        run_fastq_to_fasta(), run_rapsearch() run_ko_map(), run_SARTools(),
        run_pre_ko_remap(), run_ko_remap(), run_new_ko_remap(), run_ko_csv().
    """
    before_cwd = getcwd()
    chdir(dirname(realpath(__file__)))
    timestamp = '_'.join(str(datetime.now()).split())
    work_dir = 'metatr_results/metatr_results_' + timestamp + '/'
    system('mkdir ' + work_dir)
    system('cp -r metatr_results/metatr_pattern/* ' + work_dir)
    chdir(work_dir)
    ref_cond, all_conds, fastqs = config_from_file(work_dir, opts.metatr_config)
    print '\nconfiguration file reading: DONE\n\n'
    run_fastq_to_fasta(fastqs)
    print '\nfastq_to_fasta: DONE\n\n'
    run_cat_pairing()
    print '\ncat: DONE\n\n'
    run_rapsearch(opts.threads)
    print '\nrapsearch: DONE\n\n'
    run_ko_map()
    print '\nKO mapping: DONE\n\n'
    run_SARTools()
    print '\nSARTools: DONE\n\n'
    path_names, kopath_keys, kopath_values, edger_files, deseq_files = \
     run_pre_ko_remap(ref_cond)
    if opts.metatr_output_type != 'new':
        run_ko_remap(deseq_files, edger_files, kopath_values, path_names)
    if opts.metatr_output_type != 'old':
        ko_dict_deseq, ko_dict_edger = run_new_ko_remap(
         deseq_files, edger_files, kopath_values, all_conds, ref_cond)
    print '\npathway mapping: DONE\n\n'
    if opts.metatr_output_type != 'old':
        run_ko_csv(ko_dict_deseq, ko_dict_edger, all_conds, kopath_keys,
         path_names, ref_cond)
    print '\ngenerating summative CSV: DONE\n\n'
    if (opts.out_dir=='in_situ'):
        out_dir=before_cwd
    else:
        out_dir=opts.out_dir
    system('mkdir out_dir')
    if opts.metatr_output_type != 'new':
        old_path = out + '/old'
        system('mkdir '+out_path)
        system('cp ko_remap/* '+old_path)
    if opts.metatr_output_type != 'old':
        system('cp csv/* '+out_dir)
    system('cp edger/_report.html ' + out_dir)
    system('cp deseq/_report.html ' + out_dir)
    system('rm -rf ../metatr_results_' + timestamp + '/'
    chdir(before_cwd)
    print '\nMETATRANSCRIPTOMIC WORKFLOW DONE\n\n'
