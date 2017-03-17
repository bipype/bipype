import sqlite3
import cPickle
import subprocess
from bs4 import BeautifulSoup
from glob import glob
from urllib import urlopen
from multiprocessing import Process
from collections import Counter
from time import time
from os.path import exists as pexists
from os.path import join as pjoin
from os.path import dirname, realpath
from os import system, chdir, getcwd
from settings_bipype import *


def dicto_reduce(present, oversized):
    """Removes all elements from dictionaries,
    which keys aren't present in both.

    Args:
        present: dict
        oversized: dict

    Returns:
        (oversized, present):
            tuple of dicts

    Warning:
        Order of parametres is opposite to results.

    Example:
        >>> dict_1={'a':1,'c':3,'d':4}
        >>> dict_2={'a':3,'b':4,'c':4}
        >>> dicto_reduce(dict_1, dict_2)
        ({'a': 3, 'c': 4}, {'a': 1, 'c': 3})
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

    Args:
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

    Args:
        database: Cursor object to SQLite3 database.
    """
    database.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = database.fetchall()
    for row in tables:
        print row[0]


def auto_tax_read(db_loc):
    """Reads pickled ``{KEGG GENES number: set[KO identifiers]}`` dict."""
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
                ``{KEGG GENES identifier : set[KO identifiers]}``

        db:     Cursor object to SQL database with 'kogenes' table
                ``(KO identifier          KEGG GENES identifier)``
    Returns:
        Dict in ``{KEGG GENES identifier: set[KO identifiers]}`` format.


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

    Args:
        database: Cursor object to SQLite3 database.

    Returns:

        dict: dictionary in following format::

            {KEGG_Pathway_id:Name}

        For example::

            {
                'ko04060': 'Cytokine-cytokine receptor interaction',
                'ko00910': 'Nitrogen metabolism'
            }
    """
    database.execute('select * from Pathways')
    paths = database.fetchall()
    pathways = {}
    for path in paths:
        pathways[path[0]] = path[1]
    return pathways


def get_kopathways(database):
    """Makes dictionaries from kopathways table from SQLite3 database.

    Args:
        database: Cursor object to SQLite3 database.

    Returns:
        Two dictionaries:

            id -> pathways::

                {KO identifier: set[KEGG_Pathway_ids]}

            For example::

                {
                    'K01194': set(['ko00500','ko00600',...]),
                    'K04501': set(['ko04390',...])
                }

            pathway -> ids mappings::

                {KEGG_Pathway_id: set[KO identifiers]}

            For example::

                {ko12345: set([K12345, K12346,...]),...}
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
        `file_`:    Path to BLAST Tabular (flag: -m 8) format file
        multi_id: Dict ``{KEGG GENES identifier : set[KO identifiers]}``

    Output file (outname) has following name::

        outname = file_.replace('txt.m8', 'count')

    and following format::

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
    with open(outname, 'w') as out_file:
        for ko in tmp_ko_dict:
            to_print = '%s\t%i\n' % (ko, tmp_ko_dict[ko])
            out_file.write(to_print)
    writing_time = time()
    print (
        file_, 'comparing time seconds', writing_time - comparison_time,
        'total time', writing_time - start_time
    )


def out_content(filelist, kopath_values, path_names, method='DESeq2'):
    """For every item in 'kopath_values' dictionary and for every file
    in 'filelist', writes to output file line with KOs, which are common
    for item.value and the set of KOs obtained from file.

    Args:
        filelist:
            List of paths to tab-delimited .txt files, where
            first column is a KO identifier.

        kopath_values:
            ``{KEGG_Pathway_id:set[KO identifiers]}`` dict.

            For example::

                {ko12345:set([K12345, K12346,...]),...}

        path_names:

            Dictionary in ``{KEGG_Pathway_id:Name}`` format.

            For example::

                {
                    'ko04060': 'Cytokine-cytokine receptor interaction',
                    'ko00910': 'Nitrogen metabolism'
                }

        method:
            Argument used only as a part of output file name

    Output file has following name::

            (method+'_'+filename.replace('txt', 'path_counts.csv'))
        where:
            filename = filepath.split('\')[-1], if '\' in filepath.
            filename = filepath.split('/')[-1],  if '/' in filepath.
            filename = filepath,                 in other cases.

    anf following headline::

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
        with open(outname, 'w') as outfile:
            outfile.write('ko_path_id;ko_path_name;percent common;common KOs\n')
            for (path, Kset,) in kopath_values.items():
                common = Kids & Kset
                if len(common) > 0:
                    percent_ko = str(len(common) * 100.0 / len(Kset))
                    print_ko = ' '.join(common)
                    path_name_comma = path_names[path]
                    path_name = path_name_comma.replace(',', ' _')
                    outline = ';'.join([
                        path,
                        path_name,
                        percent_ko,
                        print_ko
                    ]) + '\n'
                    outfile.write(outline)


def fastq_to_fasta(fastq):
    """Runs fastq_to_fasta on fastq.

    GLOBAL:
        - path to fastq_to_fasta program:               PATH_FQ2FA
    """
    out_file = fastq.rsplit('.', 1)[0] + '.fasta'
    out_file = out_file.rsplit('/', 1)[-1]
    subprocess.check_call([PATH_FQ2FA, '-i', fastq, '-o', out_file, '-Q33'])


def rapsearch2(input_file, threads):
    """Runs ``rapsearch2`` for `input_file` in fasta format.

    Writes outputs in "m8/" directory.

    GLOBALS:
        - path to RAPSearch2 program:                   PATH_RAPSEARCH
        - path to similarity search database:           PATH_REF_PROT_KO
    """
    out_name = input_file.replace('tmp.fasta', 'txt')
    subprocess.check_call(
        [
            PATH_RAPSEARCH, '-q', input_file, '-d',
            PATH_REF_PROT_KO, '-o', out_name, '-z', str(threads),
            '-v', '20', '-b', '1', '-t', 'n', '-a', 't'
        ]
    )


def get_ko_fc(ko_dict, ref_cond, filepath, deseq=False):
    """From given table file (SARTool), adds found fold changes to ko_dict.

        Args:
            ko_dict:    ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict
            ref_cond:   reference condition (string)
            filepath:   filepath to output table file from edgeR or DESeq2
            deseq:      True, if filepath points to DESeq2 table file
                        False, if filepath points to edgeR table file

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
        if deseq:
            fc_index = file_content[0].rstrip().split('\t').index('FoldChange')
        else:
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
                fc = -1.0 / float(line[fc_index])
            if KO in ko_dict.keys():
                ko_dict[KO][cond] = fc
            else:
                ko_dict[KO] = {cond: fc}
    return ko_dict


def low_change(ko_dict, all_conds):
    """For every KO adds condition: 0, if condition is missing.

    Args:
        ko_dict: ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict
        all_conds: list of conditions (list of strings)

    Returns:
        suplemented ko_dict

        For example::

            low_change(
                {
                    'K12345': {'pH5': 1.41, 'pH6': 1.73},
                    'K23456': {'pH6': 2.0, 'pH8': 2.24}
                },
                ['pH5', 'pH6', 'pH8']
            )

        gives::

            {
                'K12345': {'pH5': 1.41, 'pH6': 1.73, 'pH8': 0.0},
                'K23456': {'pH5': 0.0, 'pH6': 2.0, 'pH8': 2.24}
            }

    """
    for (KO, conds,) in ko_dict.items():
        for cond in all_conds:
            if cond not in conds:
                ko_dict[KO][cond] = 0.0
    return ko_dict


def get_kegg_name(ko):
    """Returns name assigned to given KO identifier (from kegg.jp)

    Args:
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
    """Assings every KO_id from ko_dict to KEGG_Pathway_id from `ko_set`

    Args:
        ko_dict:    ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict
        ko_set:     ``{KEGG_Pathway_id:set[KO identifiers]}`` dict

    Returns:
        dict: Dict with structure::

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
    for every combination of condition & KEGG_Pathway_id.

    Args:
        ko_path_dict:
            ``{KEGG_Pathway_id:{KO_id:{cond1:value1, cond2:value2...}...}...}``
        all_conds:  list of conditions (list of strings)
        out_dir:    relative output directory path

    Output file has following path::

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
                    line = '\t'.join([
                        Kid,
                        str(ko_path_dict[pathway][Kid][cond]),
                        '\n'
                    ])
                    _file.write(line)


def config_from_file(_file):
    """Reads parameters from configuration `_file`.
    Prepares target.txt and templates for SARTools.

    Args:
        `_file`:       configuration file for metatranscriptomic pipeline

    Returns:
        (ref_cond, all_conds, fastqs):

        - ref_cond:   reference condition defined by user
        - all_conds:  set of conditions (groups) from target.txt
        - fastqs:     list of fastq files on which analysis will be done
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
                    id_num += 1
                else:
                    idents.append(hyp_ident)
                    break
            fastqs.append(line[0])
            fastqs.append(line[1])
            target_name = line[0].rsplit('/', 1)[-1]
            target_name = target_name.replace('R1_', '')
            target_name = target_name.replace('.fastq', '.count')
            f.write('\n' + hyp_ident + '_\t' + target_name + '\t' + line[2])
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
    """Runs :func:`fastq_to_fasta` for every .fastq in fastqs."""
    for fastq in fastqs:
        fastq_to_fasta(fastq)


def run_cat_pairing():
    """Merges fasta files with paired-end reads in cwd."""
    for file_R1 in glob('*R1*fasta'):
        for file_R2 in glob('*R2*fasta'):
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
    """Runs :func:`rapsearch2` for every .tmp.fasta in cwd."""
    for _file in glob('*tmp.fasta'):
        rapsearch2(_file, threads)


def run_ko_map():
    """Runs :func:`m8_to_ko` for every .m8 file in cwd.

    GLOBALS:
        - path to KO database:                                  PATH_KO_DB
        - pickle to dict from KO GENES table from KO database:  PATH_KO_PCKL
    """
    data = pickle_or_db(PATH_KO_PCKL, connect_db(PATH_KO_DB))
    p_list = []
    for file_ in glob('*m8'):
        p = Process(target=m8_to_ko, args=(file_, data))
        p.start()
        p_list.append(p)
    for p in p_list:
        p.join()


def run_SARTools():
    """Runs SARTools in R.

    HARDCODED:
        R templates:
            - edger: template_script_DESeq2.r
            - deseq: template_script_edgeR.r
    """
    system('Rscript template_script_DESeq2.r')
    system('Rscript template_script_edgeR.r')


def run_pre_ko_remap():
    """Prepares args for :func:`run_ko_remap` or :func:`run_new_ko_remap`

    Returns:
        path_names:     ``{KEGG_Pathway_id:Name}`` dict
        kopath_keys:    ``{KO identifier:set[KEGG_Pathway_ids]}`` dict
        kopath_values:  ``{KEGG_Pathway_id:set[KO identifiers]}`` dict
        edger_files:    list of edgeR outputs paths
        deseq_diles:    list of DESeq outputs paths

    HARDCODED:
        Paths to files from SARTools:
            - edger: `'edger/*[pn].txt'`
            - deseq: `'deseq/*[pn].txt'`

    GLOBALS:
        - path to KO database:  PATH_KO_DB
    """
    cursor = connect_db(PATH_KO_DB)
    path_names = get_pathways(cursor)
    (kopath_keys, kopath_values,) = get_kopathways(cursor)
    edger_files = glob('edger/tables/[pn]*.txt')
    deseq_files = glob('deseq/tables/[pn]*.txt')
    return (path_names, kopath_keys, kopath_values, edger_files, deseq_files)


def run_ko_remap(deseq_files, edger_files, kopath_values, path_names):
    """Runs ``out_content(files, kopath_values, path_names (,'edgeR'))``
    for files from `edger_paths` and `deseq_paths`.

    Args:
        deseq_diles:    list of DESeq outputs paths
        edger_files:    list of edgeR outputs paths
        kopath_values:  ``{KEGG_Pathway_id: set[KO identifiers]}`` dict
        path_names:     ``{KEGG_Pathway_id: Name}`` dict
    """
    out_content(deseq_files, kopath_values, path_names)
    out_content(edger_files, kopath_values, path_names, 'edgeR')


def run_new_ko_remap(deseq_files, edger_files, kopath_values, all_conds, ref_cond):
    """Runs :func:`get_ko_fc`, :func:`low_change`, :func:`mapper` and
    :func:`mapper_write` in appropriate way for files from
    `deseq_files` and `edger_files`.

    Args:
        deseq_diles:    list of DESeq outputs paths
        edger_files:    list of edgeR outputs paths
        ref_cond:       Reference condition (group) - string
        kopath_values:  ``{KEGG_Pathway_id:set[KO identifiers]}`` dict
        all_conds:      list of conditions (list of strings)

    Returns:
        ko_dict_deseq:  ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict
        ko_dict_edger:  ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict

    HARDCODED:
        Output directories paths:
            - deseq: 'new_ko_remap/deseq/'
            - edger: 'new_ko_remap/edger/'
    """
    ko_dict_deseq = {}
    ko_dict_edger = {}
    for _file in deseq_files:
        if ref_cond in _file:
            ko_dict_deseq = get_ko_fc(ko_dict_deseq, ref_cond, _file, True)
    for _file in edger_files:
        if ref_cond in _file:
            ko_dict_edger = get_ko_fc(ko_dict_edger, ref_cond, _file)
    ko_dict_deseq = low_change(ko_dict_deseq, all_conds)
    ko_dict_edger = low_change(ko_dict_edger, all_conds)
    mapper_deseq = mapper(ko_dict_deseq, kopath_values)
    mapper_edger = mapper(ko_dict_edger, kopath_values)
    system('mkdir new_ko_remap')
    system('mkdir new_ko_remap/deseq new_ko_remap/edger')
    mapper_write(mapper_deseq, all_conds, 'new_ko_remap/deseq/')
    mapper_write(mapper_edger, all_conds, 'new_ko_remap/edger/')
    return (ko_dict_deseq, ko_dict_edger)


def run_ko_csv(ko_dict_deseq, ko_dict_edger, all_conds, kopath_keys, path_names, ref_cond):
    """For given ko_dicts writes CSV files with pathways and foldchanges

    Args:
        ko_dict:        ``{KO_id:{cond1:value1, cond2:value2...}...}`` dict
        all_conds:      list of conditions (list of strings)
        kopath_keys:    ``{KO identifier:set[KEGG_Pathway_ids]}`` dict
        path_names:     ``{KEGG_Pathway_id:Name}`` dict
        filepath:       output filepath

    Output files have following format (and header)::
        KO_id;Gene_name;paths ids;paths names;FC vs cond1;FC vs cond2;...;

    HARDCODED:
        Output files paths:
            - deseq: 'deseq.csv'
            - edger: 'edger.csv'
    """
    for touple in [(ko_dict_deseq, 'deseq.csv'),
                   (ko_dict_edger, 'edger.csv')]:
        filepath = touple[1]
        ko_dict = touple[0]
        with open(filepath, 'wb') as _file:
            l_header = ['KO_id', 'Gene_name', 'paths ids', 'paths names']
            for cond in all_conds:
                l_header.append('FC ' + ref_cond + ' vs ' + cond)
            header = ';'.join(l_header + ['\n'])
            _file.write(header)
            for Kid in ko_dict.keys():
                to_write = [Kid]                        # ['K01369']
                to_write.append(get_kegg_name(Kid))     # ['K01369','LGMN']
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


def progress(what, estimated_percentage=None, done=True):
    """Prints specially formatted information about progress.

    Args:
        what:
            a string with name of operation which was just performed,
            and should be reported to standard output as don or failed,

        estimated_percentage: (int)

            Percent should be calculated as part of whole execution;
            first and last 5 percent should be reserved for programs
            which runs 'metatranscriptomics', for pre- and postprocessing,

        done:
            informs whether the operation from 'what' argument failed or was
            successfully done.
    """
    state = 'DONE' if done else 'FAILED'
    print '{0}: {1}'.format(what, state)
    if estimated_percentage is not None:
        print 'progress={0}'.format(estimated_percentage)


def metatranscriptomics(opts):
    """Performs analyse of metagenomic data.

    See Also:
        For more information please refer to:
            - :func:`run_fastq_to_fasta`
            - :func:`run_rapsearch`
            - :func:`run_ko_map`
            - :func:`run_SARTools`
            - :func:`run_pre_ko_remap`
            - :func:`run_ko_remap`
            - :func:`run_new_ko_remap`
            - :func:`run_ko_csv`
    """
    assert opts.out_dir != 'in_situ'
    before_cwd = getcwd()
    tmp_dir = pjoin(opts.out_dir, '.meta_tmp_results')
    if not pexists(tmp_dir) or not opts.e:
        system('mkdir ' + tmp_dir)
    chdir(tmp_dir)
    for i in ('template_script_DESeq2.r', 'template_script_edgeR.r'):
        system('cp ' + pjoin(dirname(realpath(__file__)), i) + ' .')
    ref_cond, all_conds, fastqs = config_from_file(opts.metatr_config)
    progress('configuration file reading', 15)
    if (len(fastqs) > len(glob('*_R[12]_*.fasta'))) or not opts.e:
        run_fastq_to_fasta(fastqs)
    progress('fastq_to_fasta', 25)
    if (len(fastqs)/2 > len(glob('*.tmp.fasta'))) or not opts.e:
        run_cat_pairing()
    progress('cat', 35)
    if (len(fastqs)/2 > len(glob('*.txt.m8'))) or not opts.e:
        run_rapsearch(opts.threads)
    progress('rapsearch', 45)
    if (len(fastqs)/2 > len(glob('*.count'))) or not opts.e:
        run_ko_map()
    progress('KO mapping', 55)
    if ('edger' not in glob('*')) or not opts.e:
        system('mkdir deseq edger')
        run_SARTools()
    progress('SARTools', 65)
    path_names, kopath_keys, kopath_values, edger_files, deseq_files = run_pre_ko_remap()
    if opts.metatr_output_type != 'new':
        if (
            len(all_conds) * (len(all_conds) - 1) > len(glob('*.path_counts.csv')) or
            not opts.e
        ):
            run_ko_remap(deseq_files, edger_files, kopath_values, path_names)
    if opts.metatr_output_type != 'old':
        if ('edger.csv' not in glob('*')) or not opts.e:
            ko_dict_deseq, ko_dict_edger = run_new_ko_remap(
             deseq_files, edger_files, kopath_values, all_conds, ref_cond)
            progress('pathway mapping', 75)
            run_ko_csv(ko_dict_deseq, ko_dict_edger, all_conds, kopath_keys,
             path_names, ref_cond)
    progress('generating summative CSV', 85)
    if opts.metatr_output_type != 'new':
        old_path = opts.out_dir + '/old'
        system('mkdir ' + old_path)
        system('cp *path_counts.csv ' + old_path)
    if opts.metatr_output_type != 'old':
        system('cp deseq.csv edger.csv ' + opts.out_dir)
    system('cp edger/_report.html ' + opts.out_dir + '/edger_report.html')
    system('cp deseq/_report.html ' + opts.out_dir + '/deseq_report.html')
    chdir(before_cwd)
    progress('METATRANSCRIPTOMIC WORKFLOW', 95)
