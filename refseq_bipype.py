#!/usr/bin/python
from os import getcwd, walk, system, chdir, stat
from os.path import join as pjoin
from os.path import exists as pexists
from string import split, find, join
from itertools import combinations
from datetime import datetime
from copy import deepcopy
import cPickle as pickle
# These are probably not used anymore:
from re import match
from pprint import pprint
from collections import defaultdict


def cat_read(mode, fileext, paired_end=True):
    """Returns dictionary with paths to sequence files
    from current working directory.

    Dictionary has following format:
    {directory_path1:[file1_in_this_directory_path,
                      file2_in_this_directory_path,
                      file3_in_this_directory_path],
     directory_path2:[file4_in_this_directory_path,
                      file5_in_this_directory_path,
                      file6_in_this_directory_path]
    }
    Paths to files (dictionary values) are relative.
    
    Args:
        mode:       If mode is 'run', the gunzip command will be 
                    executed to extract compressed files.
        fileext:    Extension of sequence files, which will be putted into
                    the dictionary.  
        paired_end: If True, only paired-end reads are included 
                    in dictionary (default True).
    """
    double_ext = ['contigs.fa', 'txt.m8']
    if fileext in double_ext:
        cut_len = -2
    else:
        cut_len = -1
    seq_dict = {}
    if fileext == 'fa':
        fileext_l = ['fa', 'fasta', 'FASTA', 'FA', 'fas', 'fna', 'ffn', 'frn']
    else:
        fileext_l = [fileext]
    for root, folders, files in walk(getcwd()):
        for file_ in files:
            file_ext = '.'.join(split(file_, '.')[cut_len:])
            if file_ext in fileext_l:
                seq_dict = fastq_dict(seq_dict, root, file_)
            elif file_ext == 'gz':
                prop_ext = '.'.join(split(file_, '.')[cut_len-1:-1])
                if prop_ext in fileext_l:
                    to_gunzip = pjoin(root, file_)
                    print 'gunzip %s'%(to_gunzip)
                    if mode == 'run':
                        system('gunzip %s'%(to_gunzip))
                    seq_dict = fastq_dict(seq_dict, root, file_[:-3])
            else:
                pass
    if paired_end == True:
        seq_dict = paired_end_match(seq_dict)
    else:
        pass
    return seq_dict


def exist_check(program, names, todo):
    if program == 'refseq':
        if pexists(names['sam']):
            todo.remove('bowtie')
            print '%s found'%(names['sam'])
        if pexists(names['bam']):
            todo.remove('bam_make')
            print '%s found'%(names['bam'])
        if pexists(names['sorted.bam']):
            todo.remove('sort_index')
            print '%s found'%(names['sorted.bam'])
        if pexists(names['idxstats']):
            todo.remove('idxstat')
            print '%s found'%(names['idxstats'])
        if pexists(names['map_count']):
            todo.remove('perl')
            print '%s found'%(names['map_count'])
        if pexists(names['tax_count']):
            todo = []
            print '%s exists, to get intermediate files, omit -e option'%(names['tax_count'])
    if program == 'usearch':
        if pexists(names):
            todo = []
            print '%s found'%(names)
    if program == 'usearch_0':
        if stat(names).st_size == 0:
            todo = []
            print 'PANIC MODE: %s size = 0!'%(names)
    if program == 'MV':
        if pexists(names[0]) and pexists(names[1]):
            todo  = []
            print '%s and %s found'%(names[0], names[1])
        elif pexists(names[0]):
            todo = []
            print '%s exists, to get intermediate files, omit -e option'%(names[0])
        else:
            if pexists(names[2]) and pexists(names[3]) and pexists(names[4]):
                todo.remove('velveth')
                print '%s %s %s found'%(names[2], names[3], names[4])
            if pexists(names[5]) and pexists(names[6]) and pexists(names[7]) and pexists(names[8]):
                todo.remove('velvetg')
                print '%s %s %s %s found'%(names[5], names[6], names[7], names[8])
    if program == 'rapsearch':
        if pexists(names):
            todo = []
            print '%s found, rapsearch not executing'%(names)
    if program == 'cutadapt':
        if pexists(names[0]):
            todo.remove('cutadapt0')
        if pexists(names[1]):
            todo.remove('cutadapt1')
        if pexists(names[2]):
            todo.remove('flash')
        if pexists(names[3]):
            todo.remove('fq2fa')
    return todo


def tax_id_reader():
    """Returns {GI:TaxID} dictionary from file.
    
    Keys and values are integers.
    
    File has following format:
    13	9913
    15	9915
    16	9771
    17	9771
    where first column is a GI number, second is a TaxID.
    
    HARCODED:
    /home/pszczesny/workingdata/refseq/db/tax_id/gi_taxid_nucl.dmp   
    """
    print 'reading tax ids'
    gi_tax_path = '/home/pszczesny/workingdata/refseq/db/tax_id/gi_taxid_nucl.dmp'
    gi_tax_dict = {}
    curr_t = datetime.now()
    with open(gi_tax_path) as file_:
        counter = 0
        for line in file_:
            gi, tax = split(line)
            gi_tax_dict[int(gi)] = int(tax)
            counter += 1
            if counter%1000000 == 0:
                print "Added %i tax_ids"%(counter)
    file_.close()
    return gi_tax_dict


def tax_name_reader():
    """Returns {TaxID:scientific_name} dictionary from file.
    
    Lines with not scientific names are omitted.
    
    Keys are integers, values are strings.
    
    Example (appropriate file format included):
    
        File:
    2	|	prokaryotes	|	prokaryotes <Bacteria>	|	in-part	|
    6	|	Azorhizobium	|		|	scientific name	|
    6	|	Azorhizobium Dreyfus et al. 1988	|		|	synonym	|
    6	|	Azotirhizobium	|		|	equivalent name	|
    7	|	ATCC 43989	|		|	type material	|
    7	|	Azorhizobium caulinodans	|		|	scientific name	|
     
        Output:
            {6:'Azotirhizobium', 7:'Azorhizobium caulinodans'}
         
    HARCODED:
    /home/pszczesny/workingdata/refseq/db/tax_id/names.dmp  
    """
    print 'reading tax names'
    tax_name_path = '/home/pszczesny/workingdata/refseq/db/tax_id/names.dmp'
    tax_name_dict = {}
    curr_t1 = datetime.now()
    with open(tax_name_path) as file_:
        counter = 0
        for line in file_:
            line_l = split(line, '\t|\t')
            if line_l[3] == 'scientific name\t|\n':
                tax_id = line_l[0]
                if len(line_l[2]) == 0:
                    tax_name = line_l[1]
                else:
                    tax_name = line_l[2]
                tax_name_dict[int(tax_id)] = tax_name
            counter += 1
            if counter%100000 == 0:
                print "Processed %i tax_name lines"%(counter)
    file_.close()
    return tax_name_dict


def idx_reader(file_path):
    """Reads file in samtools idxstats output format.
    
    Arg:    
        file_path: path to samtools idxstats output file
                   'File is TAB-delimited with each line consisting of
                    reference sequence name, sequence length, 
                    # mapped reads and # unmapped reads.' 
        
    Returns:
        {GI:#_mapped_reads} dictionary (keys and values are integers).  
    """
    idx_file_dict = {}
    for line in file(file_path).readlines():
        line_l = split(line)
        reads = int(line_l[2])
        if reads != 0:
            gi = int(split(line_l[0], '|')[1])
            idx_file_dict[gi] = reads
    return idx_file_dict


def idx_map(mode, file_, tax_name_dict, tax_id_dict, outfile):
    '''Parses and writes data from samtools idxstats output file.
        
    Firstly, using idx_reader(file_), create {GI:#_mapped_reads} 
    dictionary.
    Secondly, replace every GI number (key) with TaxID if appropriate
    one is available in tax_id_dict.
    Than, replace every GI number/TaxID (key) with scientific name if
    appropriate one is available in tax_name_dict.
    Finally, write data to outfile in key;value format, where
        key   is GI number/TaxID/scientific name
        value is number of mapped reads 
        
    Args:
        mode:          if (mode == 'run') function do mentioned things.
                       elif (mode == 'test') function prints
                            file_ and outfile arguments.
                       else: pass
                       
        file_:         Path to samtools idxstats output file. 
                       For more information refer to idx_reader()
                       
        tax_name_dict: {TaxID:scientific_name} dictionary.
                       For more information refer to tax_name_reader()
                                             
        tax_id_dict:   {GI:TaxID},        
                       For more information refer to tax_id_reader()
                       
        outfile:       Path to output file.
    '''
    from operator import itemgetter
    tax_count_dict = {}
    if mode == 'run':
        idx_file_dict = idx_reader(file_)
        for gi in idx_file_dict:
            try:
                tax_id = tax_id_dict[gi]
                if tax_id in tax_count_dict:
                    tax_count_dict[tax_id] += idx_file_dict[gi]
                else:
                    tax_count_dict[tax_id] = idx_file_dict[gi]
            except KeyError:
                    tax_count_dict[gi] = idx_file_dict[gi]
        tmp_dict = deepcopy(tax_count_dict)
        for tax_gi_id in tax_count_dict:
            try:
                tmp_dict[tax_name_dict[tax_gi_id]] = tmp_dict.pop(tax_gi_id)
            except KeyError:
                pass
        sorted_dict = sorted(tmp_dict.iteritems(), key=itemgetter(1), reverse=True)
        tax_write = open(outfile, 'w')
        for tax_pair in sorted_dict:
            tax_write_string = '%s;%i\n'%(tax_pair[0], tax_pair[1])
            tax_write.write(tax_write_string)
        tax_write.close()
    elif mode == 'test':
        print 'IDX_MAP %s %s'%(file_, outfile)
    else:
        pass

def pair_uni_name(file_pair):
    """ Returns one name for tuple contains pair of files.

Args:

	file_pair: tuple of paired_end read
		

Example:

	input: ('Amp15_BFk_B_p_CGTACG_L001_R1_001.fastq','Amp15_BFk_B_p_CGTACG_L001_R2_001.fastq')
	output: 'Amp15_BFk_B_p_CGTACG_L001_001'
"""
    file_body = '.'.join(split(file_pair[0], '.')[:-1])
    if 'R1' in file_body:
        R_loc = find(file_body, 'R1')
        out_name = file_body[:R_loc] + file_body[R_loc+3:]
    else:
        R_loc = find(file_body, 'R2')
        out_name = file_body[:R_loc] + file_body[R_loc+3:]
    return out_name


def refseq_ref_namespace(katalog, seq, postfix, out_dir='in_situ', map_dir='in_situ'):
    """ Return dictionary which keys are names of file extensions and values are paths to file with corresponding extension

Arg:
	(katalog, seq, postfix, out_dir='in_situ', map_dir='in_situ')
	seq: sequence file or pair of such files
	
Example:

	for seq: abcd.fastq we have output_dictionary['fastq']= pjoin(katalog, abcd.fastq)

Other keys: 'sam','sam2', 'bam', 'sorted','sorted.bam', 'idxstats', 'tax_count', 'map_count'
"""
    if out_dir == 'in_situ':
        out_dir = katalog
    ref_namespace = {}
    if type(seq) == tuple:
        sample_name = pair_uni_name(seq)
        ref_namespace['fastq'] = (pjoin(katalog, seq[0]), pjoin(katalog, seq[1]))
    else:
        ref_namespace['fastq'] = pjoin(katalog, seq)
        sample_name = split(seq, '.')[0]
    if postfix != '':
        sample_name = sample_name + '_' + postfix
    else:
        pass
    ref_namespace['sam'] = sample_name + '.sam'
    ref_namespace['sam2'] = sample_name + '_2' + '.sam'
    ref_namespace['bam'] = sample_name + '.bam'
    ref_namespace['sorted'] = sample_name + '.sorted'
    ref_namespace['sorted.bam'] = sample_name + '.sorted.bam'
    ref_namespace['idxstats'] = sample_name + '.idxstats'
    ref_namespace['tax_count'] = sample_name + '.tax_count'
    for element in ref_namespace:
        if element != 'fastq':
            ref_namespace[element] = pjoin(out_dir, ref_namespace[element])
    if map_dir == 'in_situ':
        ref_namespace['map_count'] = pjoin(out_dir, sample_name + '.map_count')
    else:
        ref_namespace['map_count'] = pjoin(map_dir, sample_name + '.map_count')
    return ref_namespace



def fastq_dict(seq_dict, root, file_): # Why not defaultdict from collections? 
    """Adds file_ to list available under seq_dict[root]

    If key 'root' is not present, function will create it.
    """
    if root in seq_dict:
        seq_dict[root].append(file_)
    else:
        seq_dict[root] = [file_]
    return seq_dict


def paired_end_match(seq_dict):
    """Returns dictionary with paths to parired-end reads only.
    
    Dictionary has following format:
    {directory_path1:[(paired-end_read1R1_path, paired-end_read1R2_path),
                      (paired-end_read2R1_path, paired-end_read2R2_path),
                      (paired-end_read3R1_path, paired-end_read3R2_path)],
     directory_path2:[(paired-end_read4R1_path, paired-end_read4R2_path)]
    }
    """ 
    pe_dict = {}
    for directory in seq_dict:
        for seq1, seq2 in combinations(seq_dict[directory], 2):
            seq1_s = set(split(seq1, '_'))
            seq2_s = set(split(seq2, '_'))
            if seq1_s^seq2_s == set(['R1', 'R2']):
                if directory in pe_dict:
                    pe_dict[directory].append((seq1, seq2))
                else:
                    pe_dict[directory] = [(seq1, seq2)]
    return pe_dict


def bowtie2_run(mode, proc, ref, out, inp1, inp2=False):
    if inp2:
        print 'bowtie2 -p %i -x %s -1 %s -2 %s -S %s'%(proc, ref, inp1, inp2, out)
        if mode == 'run':
            system('bowtie2 -p %i -x %s -1 %s -2 %s -S %s'%(proc, ref, inp1, inp2, out))
    else:
        print 'bowtie2 -p %i -x %s -f -U %s -S %s'%(proc, ref, inp1, out)
        if mode == 'run':
            system('bowtie2 -p %i -x %s -f -U %s -S %s'%(proc, ref, inp1, out))


def sam_merge(sam1, sam2):
    """ Merges two files in SAM format into one.

Args:
	sam1: SAM file
	sam2: SAM file

Returns:
	sam1 file in SAM format which contains now lines from both 		files.Function also removes sam2 file.

"""
    f1 = open(sam1, 'r')
    f2 = open(sam2, 'r')

    f1c = f1.readlines()
    f2c = f2.readlines()

    f1.close()
    f2.close()

    f1c_idx = 0
    for line in f1c:
        line = line.rstrip()
        if line[0] == '@':
            f1c_idx +=1

    f2c_idx = 0
    for line in f2c:
        line = line.rstrip()
        if line[0] == '@':
            f2c_idx += 1

    common = f1c[0:f1c_idx-1] + f2c[1:f2c_idx-1] + [f1c[f1c_idx-1]]
    cont = f1c[f1c_idx:]+ f2c[f2c_idx:]

    f3 = open(sam1, 'w')
    f3.write(join(common,''),)
    f3.write(join(cont, ''),)
    f3.close()
    remove_sam2 = 'rm %s'%(sam2)
    system(remove_sam2)


def bam_make(mode, sam, bam): #multi
    """ If mode=="run" and @SQ lines are present in the header, function converts  SAM to BAM

Args:
	mode: string type
	sam: SAM file
	bam: BAM file
"""

    print 'samtools view -bS %s > %s'%(sam, bam)
    if mode == 'run':
        system('samtools view -bS %s > %s'%(sam, bam))


def bam_sorting(mode, bam, sorted_name): #multi
    """ If mode=="run" function reads bam (file in BAM format),sort it by aligned read position and write it out to BAM file whose name is: sorted_name.

Args: 
	mode: string type
	bam: BAM file
	sorted_name: name of output file
"""
    print 'samtools sort %s %s'%(bam, sorted_name)
    if mode == 'run':
        system('samtools sort %s %s'%(bam, sorted_name))


def bam_indexing(mode, sorted_bam): #multi
    """ If mode=="run" creates an index file sorted_bam.bam.bai for the sorted_bam.bam file.

Args:

	mode: string type
	sorted_bam: BAM file
"""
    print 'samtools index %s'%(sorted_bam)
    if mode == 'run':
        system('samtools index %s'%(sorted_bam))


def bam_idxstating(mode, sorted_bam, idxstats): #multi
    """Launch samtools idxstats.
    
    Args:
        mode:       if mode=='run', launch program.
        sorted_bam: path to file in BAM format, which is the input
                    to samtools.
        idxstats:   path to file, where the output from samtools will be
                    writed.
    """
    print 'samtools idxstats %s > %s'%(sorted_bam, idxstats)
    if mode == 'run':
        system ('samtools idxstats %s > %s'%(sorted_bam, idxstats))
 
                    
def idxstat_perling(mode, idxstats, map_count): #multi
    """Counts and writes to map_count both
       (sum of all (numbers of mapped reads) in idxstats) and 
       (sum of all (numbers of unmapped reads) in idxstats).
    
    Function launches short "one-liner" in perl.
    
    Output has following format:
         #mapped - #unmapped
    For example: 123 - 456  
    
    Args:
        mode:      If mode=='run' "one-liner" is launched.
        idxstats:  Path to samtools idxstats output file.
                   Please refer to idx_reader() for more information.
        map_count: Path to file, where difference will be writed.
    """
    perl_command = """perl -e 'while(<>){chomp;@a=split "\t", $_; $b+=$a[2]; $c+=$a[3];} print "$b - $c\n";' %s > %s"""%(
        idxstats, map_count
        )
    print perl_command
    if mode == 'run':
        system(perl_command)


def refseq_mapping(mode, e, katalog, pair, postfix, refseq, tax_name_dict, tax_id_dict, threads, map_dir, refseq_2=False):
    todo = ['bowtie', 'bam_make', 'sort_index', 'idxstat', 'perl', 'idx_map']
    ref_namespace = refseq_ref_namespace(katalog, pair, postfix, 'in_situ', map_dir)
    if e:
        todo = exist_check('refseq', ref_namespace, todo)
    if 'bowtie' in todo:
        if type(pair) == tuple:
            print "\nAnalysis of %s %s vs %s\n"%(pair[0], pair[1], refseq)
            bowtie2_run(mode, threads, refseq, ref_namespace['sam'], ref_namespace['fastq'][0],
                    ref_namespace['fastq'][1])
            if refseq_2:
                print "\nAnalysis of %s %s vs %s\n"%(pair[0], pair[1], refseq_2)
                bowtie2_run(
                    mode, threads, refseq_2, ref_namespace['sam2'], ref_namespace['fastq'][0],
                    ref_namespace['fastq'][1]
                    )
        else:
            print "\nAnalysis of %s vs %s\n"%(pair, refseq)
            bowtie2_run(mode, threads, refseq, ref_namespace['sam'], ref_namespace['fastq'])
            if refseq_2:
                print "\nAnalysis of %s %s vs %s\n"%(pair[0], pair[1], refseq_2)
                bowtie2_run(mode, threads, refseq_2, ref_namespace['sam2'], ref_namespace['fastq'])
    if 'bam_make' in todo:
        if refseq_2:
            print 'Merging %s and %s'%(ref_namespace['sam'], ref_namespace['sam2'])
            if mode == 'run':
                sam_merge(ref_namespace['sam'], ref_namespace['sam2'])
        bam_make(mode, ref_namespace['sam'], ref_namespace['bam'])
    if 'sort_index' in todo:
        bam_sorting(mode, ref_namespace['bam'], ref_namespace['sorted'])
    if 'idxstat' in todo:
        bam_indexing(mode, ref_namespace['sorted.bam'])
        bam_idxstating(mode, ref_namespace['sorted.bam'], ref_namespace['idxstats'])
    if 'perl' in todo:
        idxstat_perling(mode, ref_namespace['idxstats'], ref_namespace['map_count'])
    if 'idx_map' in todo:
        idx_map(mode, ref_namespace['idxstats'], tax_name_dict, tax_id_dict, ref_namespace['tax_count'])


def ins_len_read(pair, cat):
    sample_name = pair[0]
    zeros= find(sample_name, '00') #WARNING HARDCODE
    try:
        ins_len = int(sample_name[zeros-1:zeros+2])
    except:
        aqq = len(open(pjoin(cat,sample_name), 'r').readlines()[1].rstrip())
        if aqq > 200:
            ins_len = 500
        else:
            ins_len = 200
    return ins_len


def gzip_MV(MV_dir):
    gzip_list = [
        'Sequences' ,'Roadmaps', 'PreGraph', 'Graph2', 'contigs.fa' , 'stats.txt',  'LastGraph',
        'meta-velvetg.Graph2-stats.txt', 'meta-velvetg.LastGraph', 'meta-velvetg.LastGraph-stats.txt',
        'meta-velvetg.split-stats.txt'
        ]
    to_gzip = ''
    for element in gzip_list:
        to_gzip += ' %s'%(pjoin(MV_dir, element))
        to_gzip_comm = 'tar --remove-files -czf %s %s'%(pjoin(MV_dir, 'intermediates.tgz'), to_gzip)
        system(to_gzip_comm)


def rapsearch(mode, e, contig_loc, rap_out, KEGG=None):
    rap_log = rap_out + '.log'
    rap_err = rap_out + '.err'
    rap_loc = '/home/pszczesny/soft/RAPSearch2.12_64bits/bin/rapsearch'
    if KEGG == 'masl28282828282828282828282828282828282828282828282828282828':
        ref_prot_loc = '/home/pszczesny/storage/workingdata/rapsearch/masl28'
        rap_com = '%s -q %s -d %s -o %s -z 12 -v 20 -b 1 -t n -a t 1> %s 2> %s'%(
            rap_loc, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    elif KEGG == 'KO':
        ref_prot_loc = '/home/pszczesny/soft/KEGG/ko.pep.rapsearch.db'
        rap_com = '%s -q %s -d %s -o %s -z 12 -v 20 -b 1 -t n -a t 1> %s 2> %s'%(
            rap_loc, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    else:
        ref_prot_loc = '/home/pszczesny/workingdata/refseq/protein/refseq_protein'
        rap_com = '%s -q %s -d %s  -o %s -z 10 -e 0.001 -b 100 -v 100 -g T -a T 1> %s 2> %s'%(
            rap_loc, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    print rap_com
    todo = [rap_com]
    if e:
        todo = exist_check('rapsearch', rap_out, todo)
    if (mode == 'run') and (len(todo) != 0):
        system(rap_com)


def MV(mode, e, k_mers, cat, pair, ins_len, rap=False):
    fileext = '.'.join(split(pair[0], '.')[1:])
    if ins_len == 9999:
        ins_len = ins_len_read(pair, cat)
    out_dir_local = pair_uni_name(pair) + '_velvh_out'
    out_dir = pjoin(cat, out_dir_local)
    if len(k_mers) == 1:
        k_min = int(k_mers[0])
        k_max = int(k_min)+1
        k_step = 2
    elif len(k_mers) == 2:
        k_min = int(k_mers[0])
        k_max = int(k_mers[1])
        k_step = 2
    elif len(k_mers) == 3:
        k_min = int(k_mers[0])
        k_max = int(k_mers[1])
        k_step = int(k_mers[2])
    else:
        k_min = 31
        k_max = 32
        k_step = 1
    if k_min > k_max:
        k_min, k_max = k_max, k_min
    elif k_min == k_max:
        k_max += 1
    else:
        pass
    for idx in xrange(k_min, k_max, k_step):
        todo = ['velveth', 'velvetg', 'meta']
        tmp_out_dir = out_dir + '_k-mer_' + str(idx)
        log_loc = tmp_out_dir + '/logfile'
        velveth_run = 'velveth %s %i -%s -shortPaired %s %s'%(
            tmp_out_dir, idx, fileext, pjoin(cat, pair[0]), pjoin(cat, pair[1])
            )
        velvetg_run = 'velvetg %s -exp_cov auto -ins_length %i'%(tmp_out_dir, ins_len)
        meta_run = 'meta-velvetg %s -ins_length %i | tee %s'%(tmp_out_dir, ins_len, log_loc)
        print velveth_run, '\n', velvetg_run, '\n', meta_run
        if e:
            filelist = [pjoin(tmp_out_dir, 'meta-velvetg.contigs.fa')]
            filelist.append(pjoin(tmp_out_dir, 'intermediates.tgz'))
            filelist.append(pjoin(tmp_out_dir, 'Sequences'))
            filelist.append(pjoin(tmp_out_dir, 'Roadmaps'))
            filelist.append(pjoin(tmp_out_dir, 'Log'))
            filelist.append(pjoin(tmp_out_dir, 'contigs.fa'))
            filelist.append(pjoin(tmp_out_dir, 'Graph2'))
            filelist.append(pjoin(tmp_out_dir, 'LastGraph'))
            filelist.append(pjoin(tmp_out_dir, 'stats.txt'))
            todo = exist_check('MV', filelist, todo)
        if mode == 'run':
            if 'velveth' in todo:
                system(velveth_run)
            if 'velvetg' in todo:
                system(velvetg_run)
            if 'meta' in todo:
                system(meta_run)
            gzip_MV(tmp_out_dir)
        if rap:
            rap_in = pjoin(tmp_out_dir, 'meta-velvetg.contigs.fa')
            rap_f =  pair_uni_name(pair) + '.rapsearch'
            rap_out = pjoin(tmp_out_dir, rap_f)
            rapsearch(mode, e, rap_in,  rap_out)


def usearch(mode, e, search_type, infile, database, outdir, threads):
    usearch_loc = '/home/pszczesny/soft/usearch'
    outfile = infile + '.usearch_' + search_type
    usearch_command = '%s -usearch_local %s -db %s -evalue 0.01 -id 0.9 -blast6out %s -strand both -threads %i'%(
        usearch_loc, infile, database, outfile, threads
        )
    todo = ['usearch']
    print usearch_command
    if e:
        todo = exist_check('usearch', outfile, todo)
    todo = exist_check('usearch_0', infile, todo)
    if (mode == 'run') and (len(todo) != 0):
        system(usearch_command)


def adapter_read_bck(adapter_file, filename):
    """ Function read all the lines of adapter_file and separate words in each line. If first word in line is equal filename, second and third are returned as tuple.

Arg: adapter_file, filename

Return: fin_adap_1, fin_adap_2

"""
    with open(adapter_file, 'r') as plik:
        content = plik.readlines()
        for linia in content:
            if linia:
                fastqname, adap1, adap2 = split(linia)
                if fastqname == filename:
                    fin_adap_1, fin_adap_2 = adap1, adap2
                else:
                    print filename, fastqname
                    pass
    return fin_adap_1, fin_adap_2


def adapter_read(filename):
    """ Replace 'NNNNNN' part of adp_2 by this part of the file name which contains only letters that are symbols of nucleotides.

Example:

	input: Amp18_BFp_B_p_GACGAC_L001_R1_001.fastq
	output: ('AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG')
"""
    adp_1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    adp_2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
    nameparts = split(filename, '_')
    nucleotides = set('ACTG')
    identifier = [item for item in nameparts if not set(item).difference(nucleotides)][0]
    adp_2 = adp_2.replace('NNNNNN', identifier)
    return adp_1, adp_2


def cutadapt(mode, e, cat, R1_file, R2_file, adapter_file, usearch_16S=False, usearch_ITS=False, threads=False):
    R1_fastq = pjoin(cat, R1_file)
    R2_fastq = pjoin(cat, R2_file)
    if adapter_file == 'use_filenames': #WARNING HARDCODE
        adapter_1, adapter_2 = adapter_read(R1_file)
    else:
        adapter_1, adapter_2 = adapter_read_bck(adapter_file, R1_file)
    outname_1 = '.'.join(split(R1_fastq, '.')[:-1]) + '.cutadapt.fastq'
    outname_2 = '.'.join(split(R2_fastq, '.')[:-1]) + '.cutadapt.fastq'
    outname_uni_fastq_preflash = pair_uni_name([R1_file])+ '.amplicons.cutadapt.flash.merged.fastq'
    outname_uni_fastq_postflash = outname_uni_fastq_preflash + '.extendedFrags.fastq'
    outname_uni_fasta = outname_uni_fastq_postflash[:-1] + 'a'
    names = [(outname_1, R1_fastq), (outname_2, R2_fastq)]
    outnames = [outname_1, outname_2, outname_uni_fastq_postflash, outname_uni_fasta]
    todo = ['cutadapt0', 'cutadapt1', 'flash', 'fq2fa']
    if e:
        todo = exist_check('cutadapt', outnames, todo)
    for idx in xrange(len(names)):
        cut_com = 'cutadapt -b %s -b %s -o %s %s'%(adapter_1, adapter_2, names[idx][0], names[idx][1])
        print cut_com
        if (mode == 'run') and ('cutadapt' + str(idx) in todo):
            system(cut_com)
    flash_com = 'flash %s %s -o %s'%(outname_1, outname_2, outname_uni_fastq_preflash)
    print flash_com
    if 'flash' in todo and mode == 'run':
        system(flash_com)
    fq2fa_com = '/usr/local/bioinformatics/fastx_toolkit/bin/fastq_to_fasta < %s > %s -Q33'%(outname_uni_fastq_postflash, outname_uni_fasta)
    print fq2fa_com
    if 'fq2fa' in todo and mode == 'run':
        system(fq2fa_com)
    if usearch_16S:
        usearch(mode, e, '16S', outname_uni_fasta, usearch_16S, cat, threads)
    if usearch_ITS:
        usearch(mode, e, 'ITS', outname_uni_fasta, usearch_ITS, cat, threads)


def auto_tax_read(db_loc):
    '''Reads {GI:TaxID} & {TaxID:scientific_name} dictionaries from pickle
    
    Args:
        db_loc: path to pickle file
        
    Returns:
        Two dictionaries:
          {GI:TaxID}              Keys are integers, values are strings
          {TaxID:scientific_name} Keys and values are integers
    '''
    with open(db_loc, 'rb') as fp:
        tax_names = pickle.load(fp)
        tax_id = pickle.load(fp)
    fp.close()
    return tax_id, tax_names


def taxa_read(read_mode, db_loc=None):
    '''Returns {GI:TaxID} & {TaxID:scientific_name} dictionaries.
    
    Args:
        read_mode: if (read_mode == 'manual'):
                     tax_id_reader() and tax_name_reader() functions
                     will be used and data is read from text files.
                     Refer to mentioned functions for more information.
                   else:
                     auto_tax_read(db_loc) function will be used 
                     and data is read from pickle file.
                     Refer to mentioned function for more information.
                       
        db_loc: path to pickle file, unused in 'manual' mode
                (Default = None).
        
    Returns:
        Two dictionaries:
          {GI:TaxID}              Keys are integers, values are strings
          {TaxID:scientific_name} Keys and values are integers
    '''
    if read_mode == 'manual':
        tax_id_start = datetime.now()
        tax_id_dict = tax_id_reader()
        tax_name_start = datetime.now()
        print 'tax_ids read in ', str(tax_name_start - tax_id_start)
        tax_name_dict = tax_name_reader()
        tax_name_end = datetime.now()
        print 'tax_names read in ', str(tax_name_end - tax_name_start)
    else:
        tax_id_dict, tax_name_dict = auto_tax_read(db_loc)
    return tax_id_dict, tax_name_dict


def reconstruct(mode, thr, e, pair, cat, prefix, rec_db_loc):
    tmp_loc = '/tmp/' + '_'.join(str(datetime.now()).split())
    mkdir = "mkdir %s"%(tmp_loc)
    fq1 = pjoin(cat, pair[0])
    fq2 = pjoin(cat, pair[1])
    reads_fastq = pjoin(tmp_loc, 'reads.fastq')
    merge = "cat %s %s > %s"%(fq1, fq2, reads_fastq)
    reads_sai = pjoin(tmp_loc, 'reads.sai')
    bwa_aln_com = "bwa aln -t %i %s %s > %s"%(thr, rec_db_loc, reads_fastq, reads_sai)
    alnsam = pjoin(tmp_loc, 'aln.sam')
    bwa_samse_com = "bwa samse %s %s %s > %s"%(rec_db_loc, reads_sai, reads_fastq, alnsam)
    alnbam = pjoin(tmp_loc, 'aln.bam')
    alnsrt = pjoin(tmp_loc, 'aln.sorted')
    alnsrtbam = pjoin(tmp_loc, 'aln.sorted.bam')
    if prefix[-1] != '_':
        prefix = prefix+'_'
    alnfastq = pjoin(cat, prefix + 'aln.fastq')
    pile_com = "samtools mpileup -uf %s %s |bcftools view -cg - |vcfutils.pl vcf2fq > %s"%(rec_db_loc, alnsrtbam, alnfastq)
    rm_tmp = "rm -rf %s"%(tmp_loc)

    print "%s"%(mkdir)
    system("%s"%(mkdir))
    print "%s"%(merge)
    print "%s"%(bwa_aln_com)
    print "%s"%(bwa_samse_com)
    if mode == 'run':
        system("%s"%(merge))
        system("%s"%(bwa_aln_com))
        system("%s"%(bwa_samse_com))
    bam_make(mode, alnsam, alnbam)
    bam_sorting(mode, alnbam, alnsrt)
    bam_indexing(mode, alnsrtbam)
    print "%s"%(pile_com)
    if mode == 'run':
        system("%s"%(pile_com))
    print "%s"%(rm_tmp)
    system("%s"%(rm_tmp))


def humann(mode, e, m8_dict, typ='m8'):
    for path in m8_dict:
        flag = True
        if typ == 'm8':
            hum_loc = pjoin(path, 'humann-0.99')
        else:
            hum_loc = path
            flag = None
        if e:
            if pexists(hum_loc):
                print 'Found %s'%(hum_loc)
                flag = None
        if flag:
            get_humann = 'cp -r /home/pszczesny/soft/humann-0.99 %s'%(path)
            print get_humann
            if mode == 'run':
                system(get_humann)
        input_loc = pjoin(hum_loc, 'input')
        if typ == 'm8':
            chdir(hum_loc)
        else:
            curr_loc = hum_loc[:-6]
            chdir(curr_loc)
        print 'I am in here: %s'%(getcwd())
        if flag:
            mv_m8 = 'mv %s %s'%(pjoin(path, m8_dict[path][0]), pjoin(input_loc, m8_dict[path][0][:-3].replace('_', '-')))
            print mv_m8
            clean_comm = 'rm -rf %s/*'%(input_loc)
            print 'Executing: %s'%(clean_comm)
        if typ == 'm8':
            hmp_get_comm = 'cp /home/pszczesny/soft/humann-0.99/input/hmp_metadata.dat %s/'%(input_loc)
        else:
            hmp_get_comm = 'cp /home/pszczesny/soft/humann-0.99/input/hmp_metadata.dat %s/'%(hum_loc)
        print 'Getting hmp: %s'%(hmp_get_comm)
        humann_comm = 'scons 1>log 2>errlog'
        print 'Executing: %s'%(humann_comm)
        if mode == 'run':
            if flag:
                system(clean_comm)
                system(mv_m8)
            system(hmp_get_comm)
            system(humann_comm)


def sample(opts):
    if opts.to_calculate == None:
        opts.to_calculate = []
    refseq_test = len(set(['f', 'b', 'p'])&set(opts.to_calculate))
    if (refseq_test > 0) and opts.mode == 'run':
        if pexists(opts.db_NCBI_taxonomy):
            tax_id_dict, tax_name_dict = taxa_read('auto', opts.db_NCBI_taxonomy)
        else:
            tax_id_dict, tax_name_dict = taxa_read('manual')
    else:
        tax_id_dict = tax_name_dict = {}
    fastq_dict = cat_read(opts.mode, 'fastq')
    for cat in fastq_dict:
        for pair in fastq_dict[cat]:
            if opts.reconstruct:
                reconstruct(opts.mode, opts.threads, opts.e, pair, cat, opts.reconstruct, opts.db_reconstruct)
            if opts.MV != None:
                if 'rap_prot' in opts.to_calculate:
                    MV(opts.mode, opts.e, opts.MV, cat, pair, opts.ins_len, 1)
                else:
                    MV(opts.mode, opts.e, opts.MV, cat, pair, opts.ins_len)
            if ('f' in opts.to_calculate) or ('b' in opts.to_calculate):
                postfix = opts.postfix + 'fungi'
                refseq_mapping(
                    opts.mode, opts.e, cat, pair, postfix, opts.db_refseq[0], tax_name_dict, tax_id_dict,
                    opts.threads, opts.out_dir
                    )
            if ('p' in opts.to_calculate) or ('b' in opts.to_calculate):
                postfix = opts.postfix + 'plant'
                refseq_mapping(
                    opts.mode, opts.e, cat, pair, postfix, opts.db_refseq[1], tax_name_dict,
                    tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq[2]
                    )
            if opts.cutadapt != '':
                if ('both' in opts.cutadapt) or (('16S' in opts.cutadapt) and ('ITS' in opts.cutadapt)):
                    cutadapt(opts.mode, opts.e, cat, pair[0], pair[1], opts.cutadapt[0],
                        opts.db_16S, opts.db_ITS, opts.threads
                        )
                elif '16S' in opts.cutadapt:
                    cutadapt(opts.mode, opts.e, cat, pair[0], pair[1], opts.cutadapt[0],
                        opts.db_16S, False, opts.threads
                        )
                elif 'ITS' in opts.cutadapt:
                    cutadapt(opts.mode, opts.e, cat, pair[0], pair[1], opts.cutadapt[0],
                        False, opts.db_ITS, opts.threads
                        )
                else:
                    cutadapt(opts.mode, opts.e, cat, pair[0], pair[1], opts.cutadapt[0])

    if opts.MV != None:
        contig_dict = cat_read(opts.mode, 'fa', False)
        for cat in contig_dict:
            for contig in contig_dict[cat]:
                if ('f' in opts.to_calculate) or ('b' in opts.to_calculate):
                    postfix = opts.postfix + 'fungi'
                    refseq_mapping(
                        opts.mode, opts.e, cat, contig, postfix, opts.db_refseq[0], tax_name_dict,
                        tax_id_dict, opts.threads, opts.out_dir
                        )
                if ('p' in opts.to_calculate) or ('b' in opts.to_calculate):
                    postfix = opts.postfix + 'plant'
                    refseq_mapping(
                        opts.mode, opts.e, cat, contig, postfix, opts.db_refseq[1], tax_name_dict,
                        tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq[2]
                        )
    else:
        if 'rap_prot' in opts.to_calculate:
            contig_dict = cat_read(opts.mode, 'contigs.fa', False)
            for cat in contig_dict:
                for contig in contig_dict[cat]:
                    rap_in = pjoin(cat, contig)
                    rap_out = pjoin(cat, split(cat, '/')[-1] + '.rapsearch')
                    rapsearch(opts.mode, opts.e, rap_in,  rap_out)

    if ('rap_KEGG' in opts.to_calculate) or ('rap_KO' in opts.to_calculate):
        fastq_dict = cat_read(opts.mode, 'fastq')
        fasta_dict = {}
        for path in fastq_dict:
            fasta_dict[path] = []
            for fastq_pair in fastq_dict[path]:
                for fastq in fastq_pair:
                    fasta = fastq[:-1] + 'a'
                    flag = True
                    if opts.e:
                        if pexists(pjoin(path, fasta)):
                            print 'Found %s'%(pjoin(fasta))
                            flag = None
                    if flag:
                        fq2fa_com = '/usr/local/bioinformatics/fastx_toolkit/bin/fastq_to_fasta < %s > %s -Q33'%(pjoin(path, fastq), pjoin(path, fasta))
                        print fq2fa_com
                        if opts.mode == 'run':
                            system(fq2fa_com)
                    fasta_dict[path].append(fasta)
        rap_paths = []
        for path in fasta_dict:
            # WARNING!!! HARDCODE
            if 'seq' in path:
                new_path = '_'.join(split(path,'_')[:-1])
                if pexists(new_path):
                    pass
                else:
                    system('mkdir %s'%(new_path))
            else:
                new_path = path
            rap_in =  pjoin(new_path, 'tmp.fasta')
            rap_out = pjoin(new_path, split(new_path, '/')[-1]+'.txt')
            rap_paths.append((rap_in, rap_out))
            for fasta in fasta_dict[path]:
                flag = True
                if opts.e:
                    if pexists(rap_in):
                        print 'Found %s'%(rap_in)
                        flag = False
                if flag:
                    cat_com = 'cat %s >> %s'%(pjoin(path, fasta), rap_in)
                    print cat_com
                    if opts.mode == 'run':
                        system(cat_com)
        if 'rap_KEGG' in opts.to_calculate:
            kegg_db = 'masl'
        else:
            kegg_db = 'KO'
        for pair in rap_paths:
            rapsearch(opts.mode, opts.e, pair[0], pair[1], kegg_db)

    if 'humann' in opts.to_calculate:
        m8_dict = cat_read(opts.mode, 'txt.m8', False)
        if len(m8_dict) > 0:
            humann(opts.mode, opts.e, m8_dict)
        else:
            txt_dict = cat_read(opts.mode, 'txt', False)
            tmp_dict = txt_dict.copy()
            for key in txt_dict:
                if 'input' not in key:
                    del tmp_dict[key]
            humann(opts.mode, opts.e, tmp_dict, 'nom8')

    if '16S' in opts.to_calculate or 'ITS' in opts.to_calculate:
        fasta_dict = cat_read(opts.mode, 'fa', False)
        for cat in fasta_dict:
            for fasta in fasta_dict[cat]:
                if '16S' in opts.to_calculate:
                    usearch(opts.mode, opts.e, '16S', pjoin(cat, fasta), opts.db_16S, cat, opts.threads)
                if 'ITS' in opts.to_calculate:
                    usearch(opts.mode, opts.e, 'ITS', pjoin(cat, fasta), opts.db_ITS, cat, opts.threads)


def SSU_read(loc, typ=None):
    """Extracts from specially formatted FASTA file taxonomical data
    and returns them as hierarchically organised dict.

    Args:
        loc: location - path to the fasta file database

        typ: None or '16S' TODO

    Returns:
        A dict with taxonomical data in format:

    """
    # TODO: possible memory leak: file loc might be not closed.
    # TODO: It might be better with "with"
    tax_dict = {}
    for linia in open(loc, 'r').readlines():
        if linia[0] == '>':
            linia = linia.rstrip()
            if typ:
                linia = linia.split()
                tax_idx = linia[0][1:]
                tax_dict[tax_idx] = linia[1].split(';')
            else:
                line = linia.split('|')
                tax_idx = line[0][1:]
                if linia.count('Fungi') == 2:
                    tax_line = line[-1].replace(" ", "")
                elif linia.count('Fungi') == 1:
                    tax_line = line[2].replace(" ", "")
                    if len(tax_line.split(';')) == 1:
                        tax_line = line[-1].replace(" ", "")
                if tax_line == '-':
                    tax_line = 'Fungi'+';'+linia.split('|')[1].replace(" ", "_")
                tax_dict[tax_idx] = tax_line.split(';')
    return tax_dict


def tuple_to_dict(tuple_dict):
    fin_dict = {}
    for key in tuple_dict.keys():
        curr_dict = fin_dict
        val = tuple_dict[key]
        for i in xrange(len(key)):
            if key[i] in curr_dict:
                if i == len(key)-1:
                    if 'subsum' in curr_dict[key[i]]:
                        curr_dict[key[i]]['subsum'] += val
                    else:
                        curr_dict[key[i]]['subsum'] = val
                else:
                    curr_dict = curr_dict[key[i]]
            else:
                if i == len(key)-1:
                    curr_dict[key[i]] = {'subsum': val}
                else:
                    curr_dict[key[i]] = {}
                    curr_dict = curr_dict[key[i]]
    return fin_dict


def tuple_to_xml_dict(tuple_dict):
    fin_dict = {}
    for taxonomy in tuple_dict:
        for level in taxonomy:
            if level in fin_dict.keys():
                fin_dict[level] += tuple_dict[taxonomy]
            else:
                fin_dict[level] = tuple_dict[taxonomy]
    return fin_dict


def dict_purify(bac_dict):
    treshold = 10 #WARNING HARDCODE
    all_keys = bac_dict.keys()
    to_remove = []
    for key in all_keys:
        if bac_dict[key] < treshold:
            to_remove.append(key)
    for klucz in to_remove:
        del bac_dict[klucz]
    return bac_dict


def file_analysis(typ, name, SSU=None):
    import operator
    if not pexists(name):
        return 'NA', 'NA'
    else:
        with open(name, 'r') as content:
            linie = content.readlines()
            if typ in ['fmc', 'pmc']:
                mapped, unmapped =  linie[0].strip().split('-')
                try:
                    mapped = int(mapped)
                except:
                    mapped = 'NA'
                try:
                    unmapped = int(unmapped)
                except:
                    unampped = 'NA'
                if mapped == 'NA' or unmapped == 'NA':
                    total = 'NA'
                else:
                    total = mapped + unmapped
                return mapped, total
            if typ == 'log':
                data = linie[-1].split(' ')
                try:
                    n50 = int(data[8][:-1])
                    total = int(data[12][:-1])
                    using =split(data[14], '/')[0]
                except:
                    n50 = total = using = 'NA'
                return n50, total, using
            if typ == '16S':
                bac_arch = [0, 0]
                bac_dict = {}
                for linia in linie:
                    tax_id = linia.split('\t')[1]
                    cult_control = 0
                    spec_idx = -1
                    while cult_control == 0:
                        if SSU[tax_id][0] == 'Eukaryota':
                            cult_control = 1
                            spec = False
                        elif 'uncultured' in SSU[tax_id][spec_idx] or 'unidentified' in SSU[tax_id][spec_idx]:
                            spec_idx -= 1
                        elif 'Candidate' in SSU[tax_id][spec_idx]:
                            spec = tuple(SSU[tax_id][:spec_idx+1])
                            cult_control = 1
                        else:
                            spec=tuple(SSU[tax_id][:spec_idx+1])
                            cult_control = 1
                    if SSU[tax_id][0] == 'Bacteria':
                        bac_arch[0] += 1
                    elif SSU[tax_id][0] == 'Archaea':
                        bac_arch[1] += 1
                    else:
                        pass
                    if spec:
                        if spec == 'Phaseolus_acutifolius_(tepary_bean)':#HARDCODE!!!
                            pass
                        else:
                            if spec not in bac_dict:
                                bac_dict[spec] = 1
                            else:
                                bac_dict[spec] += 1
                    else:
                        pass
                cut_dict = dict_purify(bac_dict)
                full_bac_dict = tuple_to_dict(cut_dict)
                tax_bac_dict = tuple_to_xml_dict(cut_dict)
                return full_bac_dict, tax_bac_dict
            if typ == 'ITS':
                fun_dict = {}
                for linia in linie:
                    tax_id = linia.split('\t')[1].split('|')[0]
                    spec_name = tuple(SSU[tax_id])
                    if spec_name not in fun_dict:
                        fun_dict[spec_name] = 1
                    else:
                        fun_dict[spec_name] += 1
                cut_dict = dict_purify(fun_dict)
                full_fun_dict = tuple_to_dict(cut_dict)
                tax_fun_dict = tuple_to_xml_dict(cut_dict)
                return full_fun_dict, tax_fun_dict


def input_locations(mode, out_types):
    """

    Args:

    Returns:


    """
    in_dict = {}
    for out_type in out_types:
        if out_type in ['ITS', '16S']:
            in_dict[out_type] = cat_read(mode, 'usearch_' + out_type, False)
        else:
            in_dict[out_type] = cat_read(mode, out_type, False)
    return in_dict


def dict_prepare(typ, indict, SSU):
    """

    Args:
        out_type: an "output type" - typically: 'ITS', '16S' or 'txt'

        indict:

        SSU: a dict with taxonomic information in format the same as generated by
            SSU_read() function, it est:
                TODO


    Returns:


    """
    all_dicts = {}
    all_tax_dict = {}
    for path in indict.keys():
        for plik in indict[path]:
            plik_id = '_'.join([typ, split(plik,'.')[0]])
            analysed_dict, tax_dict = file_analysis(typ, plik, SSU)
            all_dicts[plik_id] = analysed_dict
            all_tax_dict[plik_id] = tax_dict
    return all_dicts, all_tax_dict


def dict_sum(dicto, val):
    """Sums the values in the given dict,
    using recurrence to handle nested dicts.

    Args:
        dicto: a dict, where values are of any type with defined + operation
        val: an int - initial value

    Returns:
        val: an int - sum of all values in given dict plus
            value given as second argument
    """
    for key in dicto.keys():
        if isinstance(dicto[key], dict):
            val = dict_sum(dicto[key], val)
        else:
            val += dicto[key]
    return val


def update_dict(tax_tree, curr_tax):
    """Updates the dict with keys form another one.
    If during walking the dict, the 'subsum' string is found (as key or value),
    recurrence stops and nodes deeper inside that branch are ignored.
    The dicts have to contain only other dicts or 'subsum'!

    Args:
        tax_tree: a dict to update
        curr_tax: a dictionary to be added into tax_tree

    Returns:
        tax_tree: updated dict

    Example:
        input: {'a':{'a':{}}},{'a':{'b':{}},'b':{},'c':'subsum','subsum':{}}
        output: {'a': {'a': {}, 'b': {}}, 'b': {}, 'c': {}}
    """
    for key in curr_tax:
        if key != 'subsum':
            if key not in tax_tree:
                tax_tree[key] = {}
                if curr_tax[key] != 'subsum':
                    tax_tree[key] = update_dict(tax_tree[key], curr_tax[key])
            else:
                tax_tree[key] = update_dict(tax_tree[key], curr_tax[key])
    return tax_tree


def tree_of_life(full_dict):
    """Rewrites the data from a dict containing information arranged by
    type and then by file into two dicts: full_tree and file_total_count.

    Args:
        full_dict: a dict of dicts, with structure:
            {type_1: {file_1: {}, file_2: {}}, type_2: {file_3: {}}}

    Returns:
        full_tree: A dict without information about type and file.
        file_total_count: A dict where:
            - information about files are kept,
            - values of every child dict are summed into a single value
    """
    file_total_count = {}
    full_tree = {}
    for typ in full_dict:
        for plik in full_dict[typ]:
            file_total_count[plik] = dict_sum(full_dict[typ][plik], 0)
            curr_tax = full_dict[typ][plik]
            full_tree = update_dict(full_tree, curr_tax)
    return full_tree, file_total_count


def xml_name_parse(full_dic):
    """

    Args:

    Returns:


    """
    name_set = set()
    for typ in full_dic:
        for plik in full_dic[typ]:
            display_name = split(plik, '.')[0]
            try:
                nameparts = split(plik, '_')
                nucleotides = set('ACTG')
                identifier = [item for item in nameparts if not set(item).difference(nucleotides)][0]
                display_name = '_'.join(nameparts[0:nameparts.index(identifier+1)])
            except:
                pass
            name_set.add(display_name)
    name_list = list(name_set)
    return name_list


def xml_vals(xml_names, tax_dict):
    """

    Args:

    Returns:


    """
    all_tax_set = set()
    for plik in tax_dict:
        for level in tax_dict[plik]:
            all_tax_set.add(level)
    xml_dict = {}
    for tax in all_tax_set:
        xml_dict[tax] = []
        for name in xml_names:
            for plik in tax_dict:
                if name in plik:
                    if tax in tax_dict[plik]:
                        xml_dict[tax].append(tax_dict[plik][tax])
                    else:
                        xml_dict[tax].append(0)
    return xml_dict


def prettify(elem):
    """Creates a string representation of an XML element by calling
    xml.etree.ElementTree.tostring() and then, parses this string again
    to obtain pretty-printed XML string, where indents are four spaces long.

    Args:
        elem: instance of class Element from xml.etree.ElementTree

    Returns:
        A string containing pretty-printed XML of <elem>
    """
    import xml.etree.ElementTree as ET
    from xml.dom import minidom
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxmltoprettyxml(indent="    ")


def dict_depth(d, depth=0):
    # TODO: needed only if line "recurse_depth = dict_depth(tax_tree)"
    # TODO: from xml_prepare will be kept.
    if not isinstance(d, dict) or not d:
        return depth
    return max(dict_depth(v, depth+1) for k, v in d.iteritems())


def deunique(node):
    """Extracts name of node from string representing taxonomic position
    of this node. Compare with linia_unique().

    Example:
        input: 'lvl1_____lvl2_____lvl_3_____lvln'
        output: 'lvln'

    Args:
        node: a string with a list of ancestors of the node
            and the node itself, separated by '_____'.

    Returns:
        A string with pure name of the node.
    """
    clean_node = node.split('_____')[-1]
    return clean_node


def xml_prepare(xml_names, xml_dict, tax_tree, name_total_count, unit='reads'):
    #TODO
    """
    Args:
        xml_names: identificators, obtained by modifying filenames,
            given by input_d dict.
            Example:
            [identificator_1, identificator_2, ..., identificator_x]

        xml_dict: A dict. Contains lists of numbers of occurrences of nodes by node. Keys are names of nodes.
            Values are lists of constant length, where every item on position *X* means:
            number of occurrences of nodes with prefix equal to *name of node* in file *X*.
            *file X* is that file, which is on position *X* on list xml_names.
             Note, that for all i: len(xml_dict[i]) is equal to len(xml_names).
            Example:
            {node: [count_1, count_2, ... count_x], node_2: [count_1, count_2, ... count_x]}

        tax_tree: A dict of dicts, etc. In form of nested dicts, it represents phylogenetic tree.
            Example:
            {a:{aa:{}, ab:{aba:{}, abb:{abba:{}}}}}
            Explicit example in description of tax_tree_graphlan, look for: total_tax_tree.

        name_total_count: A dict. Contains summed numbers of occurrences of nodes by file.
            Keys are identificators (derived from filenames - look for xml_names), values ale sums.
            Example:
            {identificator_1: count_1, identificator_2: count_2}

        unit: A string. Defines units to be coded in Krona XML. Default = 'reads'


    Returns:
    """
    import xml.etree.ElementTree as ET
    krona = ET.Element('krona')
    attributes = ET.SubElement(krona, 'attributes', magnitude=unit)
    attribute = ET.SubElement(attributes, 'attribute', display=unit)
    attribute.text = unit
    datasets = ET.SubElement(krona, 'datasets')
    for name in xml_names:
        dataset = ET.SubElement(datasets, 'dataset')
        dataset.text = name
    node = ET.SubElement(krona, 'node', name='%s count'%(unit))
    reads = ET.SubElement(node, unit)
    for name in xml_names:
        val = ET.SubElement(reads, 'val')
        val.text = str(name_total_count[name])
    recurse_depth = dict_depth(tax_tree)    # TODO: this line is not needed
    # WARNING HARDCODE - however it ignores all absent levels!
    # which is actually kind of cool
    for lvl_1 in tax_tree:
        node1 = ET.SubElement(node, 'node', name=deunique(lvl_1))
        reads1 = ET.SubElement(node1, unit)
        for element in xml_dict[lvl_1]:
            val = ET.SubElement(reads1, 'val')
            val.text = str(element)
        tmp_2_d = tax_tree[lvl_1]
        for lvl_2 in tmp_2_d:
            node2 = ET.SubElement(node1, 'node', name=deunique(lvl_2))
            reads2 = ET.SubElement(node2, unit)
            for element in xml_dict[lvl_2]:
                val = ET.SubElement(reads2, 'val')
                val.text = str(element)
            tmp_3_d = tmp_2_d[lvl_2]
            for lvl_3 in tmp_3_d:
                node3 = ET.SubElement(node2, 'node', name=deunique(lvl_3))
                reads3 = ET.SubElement(node3, unit)
                for element in xml_dict[lvl_3]:
                    val = ET.SubElement(reads3, 'val')
                    val.text = str(element)
                tmp_4_d = tmp_3_d[lvl_3]
                for lvl_4 in tmp_4_d:
                    node4 = ET.SubElement(node3, 'node', name=deunique(lvl_4))
                    reads4 = ET.SubElement(node4, unit)
                    for element in xml_dict[lvl_4]:
                        val = ET.SubElement(reads4, 'val')
                        val.text = str(element)
                    tmp_5_d = tmp_4_d[lvl_4]
                    for lvl_5 in tmp_5_d:
                        node5 = ET.SubElement(node4, 'node', name=deunique(lvl_5))
                        reads5 = ET.SubElement(node5, unit)
                        for element in xml_dict[lvl_5]:
                            val = ET.SubElement(reads5, 'val')
                            val.text = str(element)
                        tmp_6_d = tmp_5_d[lvl_5]
                        for lvl_6 in tmp_6_d:
                            node6 = ET.SubElement(node5, 'node', name=deunique(lvl_6))
                            reads6 = ET.SubElement(node6, unit)
                            for element in xml_dict[lvl_6]:
                                print 'lvl6', lvl_6, xml_dict[lvl_6]
                                val = ET.SubElement(reads6, 'val')
                                val.text = str(element)
                            tmp_7_d = tmp_6_d[lvl_6]
                            for lvl_7 in tmp_7_d:
                                node7 = ET.SubElement(node6, 'node', name=deunique(lvl_7))
                                reads7 = ET.SubElement(node7, unit)
                                for element in xml_dict[lvl_7]:
                                    print 'lvl7', lvl_7, xml_dict[lvl_7]
                                    val = ET.SubElement(reads7, 'val')
                                    val.text = str(element)
                                tmp_8_d = tmp_7_d[lvl_7]
                                for lvl_8 in tmp_8_d:
                                    node8 = ET.SubElement(node7, 'node', name=deunique(lvl_8))
                                    reads8 = ET.SubElement(node8, unit)
                                    for element in xml_dict[lvl_8]:
                                        print 'lvl8', lvl_8, xml_dict[lvl_8]
                                        val = ET.SubElement(reads8, 'val')
                                        val.text = str(element)
                                    tmp_9_d = tmp_8_d[lvl_8]
                                    for lvl_9 in tmp_9_d:
                                        node9 = ET.SubElement(node8, 'node', name=deunique(lvl_9))
                                        reads9 = ET.SubElement(node9, unit)
                                        for element in xml_dict[lvl_9]:
                                            val = ET.SubElement(reads9, 'val')
                                            val.text = str(element)

    return prettify(krona)


def name_total_reduction(xml_names, file_total_count):
    """

    Args:

    Returns:


    """
    name_total_count = {}
    for name in xml_names:
        count = 0
        for plik in file_total_count:
            if name in plik:
                count += file_total_count[plik]
            name_total_count[name] = count
    return name_total_count


def xml_format(full_dict, tax_dict):
    """

    Args:
        full_dict: a dict of dicts with structure like:
            {
                type_1:
                    {file_1: {}, file_2: {}},
                type_2:
                    {file_3: {}}}
            }
            Type is an "output type" and typically is 'ITS', '16S' or 'txt'.

    Returns:


    """
    tax_tree, file_total_count = tree_of_life(full_dict)
    xml_names = xml_name_parse(full_dict)
    xml_dict = xml_vals(xml_names, tax_dict)
    name_total_count = name_total_reduction(xml_names, file_total_count)
    return xml_names, xml_dict, tax_tree, name_total_count


def out_namespace(curr, out_type):
    """"Generates paths to xml and html files, where krona xml and krona html
    files are or will be placed. Path composition varies, according to given
    arguments and always contains some string representation of <out_type>.
    Path may starts in current working directory or in any give directory.
    A basename is generated accordingly to the following table:

        out_type                  basename (filename without extension)
        'txt'                    'humann-graphlan'
        a str other than 'txt'   out_type
        a list of strings        '_'.join(out_type)


    Args:
        curr: A string - path to directory where the files are/will be placed.
              If == 'in_situ', the path will start in current working directory.

        out_type: A variable that determines the basename of file.

            In practice, we expect here: 'ITS', '16S', 'txt'
            or some combination of first two in form of a list.

    Returns:
        A tuple (out_xml, out_html):

        out_xml: a string with path to file. The filename is derived from
        (out_type) and file has '.krona' extension.

        out_html: a string with path to file. The filename is derived from
        (out_type) and file has '.html' extension.
    """
    if curr == 'in_situ':
        curr = getcwd()
    if type(out_type) == str:
        if out_type != 'txt':
            out_xml = pjoin(curr, out_type+'.krona')
        else:
            out_xml = pjoin(curr, 'humann-graphlan.krona')
    elif type(out_type) == list:
        out_xml = pjoin(curr, '_'.join(out_type)+'.krona')
    else:
        print 'Unknown type of output type:', out_type, ' is ', type(out_type)
    out_html = out_xml.replace('krona', 'html')
    return out_xml, out_html


def outprint(xml_string, out_xml):
    """Writes xml_string into the file of name out_xml.
    If the file does not exist, creates a new file.

    Args:
        xml_string: a string with content to write

        out_xml: a string with path to file
    """
    krona_write = open(out_xml, 'w')
    krona_write.write(xml_string)
    krona_write.close()


def txt_dict_clean(dicto):
    """Cleans given file structure by removing files
    which do not contain words "graplan" or "tree".
    Emptied directories are also removed.

    Args:
        dicto (dict of str: str): A dict where:
            keys are names of directories,
            values are lists containing filenames that are in the directory.
            Example:
            {'dict_name_1': ['file_1','file_2'], 'dict_name_2': ['file_1']}

    Returns:
        A dict in the same format as given one (cleaned from unwanted files)
    """
    to_remove = set()
    for directory in dicto:
        for plik in dicto[directory]:
            if ('graphlan' not in plik) or ('tree' not in plik):
                to_remove.add((directory, plik))
    for pair in to_remove:
        dicto[pair[0]].remove(pair[1])
    dir_to_del = []
    for directory in dicto:
        if len(dicto[directory]) == 0:
            dir_to_del.append(directory)
    for directory in dir_to_del:
        del dicto[directory]
    return dicto

def xml_names_graphlan(input_d):
    """Extracts values from given dict, puts the values on a list
     and then, removes from every value those parts,
     which are prefixes and suffixes, common for all items on the list.

    Args:
        input_d (dict of str: str): A dict, where:
            keys are names of directories,
            values are lists containing filenames that are inside the directory
    Returns:
        A list with all filenames with were given on input,
        trimmed by common prefixes and suffixes.

    """
    from os.path import commonprefix as CP

    # create a set (and then a list) with all names of files given on input
    allnames = set()
    for path in input_d:
        files = set(input_d[path])
        allnames = allnames | files
    allnames = list(allnames)

    comm_pref = CP(allnames)

    allrev = []
    for plik in allnames:
        allrev.append(plik[::-1])

    comm_suff = CP(allrev)[::-1]

    # TODO: a co jeśli comm_pref lub comm_suff występują też wewnątrz stringu? W sumie to rare-case, ale
    # TODO: np pliki: ./axxbxx i ./abxxx
    # remove common prefix and suffix from every filename
    for idx in xrange(len(allnames)):
        file_name = allnames[idx]
        file_name = file_name.replace(comm_pref, '')
        file_name = file_name.replace(comm_suff, '')
        allnames[idx] = file_name

    return allnames


def tax_tree_extend(tax_tree, linia):
    """Recursively adds elements form list <linia> into the tree tax_tree.

    Args:
        tax_tree: tree represented as dict of dicts.

        linia: a list in order: from the oldest ancestor to the youngest descendant.
            Elements should be formatted to contain information
            about ancestors (like it does linia_unique function).

    Returns:
        Extended tree in form of dict of dicts.
    """
    if linia[0] not in tax_tree:
        tax_tree[linia[0]] = {}
    if len(linia) > 1:
        tax_tree[linia[0]] = tax_tree_extend(tax_tree[linia[0]], linia[1:])
    return tax_tree

def linia_unique(linia):
    """Creates list where every element includes
    information about previous elements, in order, separated by "____".

    Compare with: deunique

    Example:
        input: ['1','2','3'],
        output: ['1', '1_____2', '1_____2_____3']

    Args:
        linia: a list.

    Returns:
        A list.
   """
    line = linia[:]
    for idx in xrange(len(linia)):
        line[idx] = '_____'.join(linia[:idx+1])
    return line

def tax_tree_graphlan(input_d):
    # TODO: Currently Graphlan files have to contain only taxonomic data (formatting data after \t are not allowed)
    # TODO: Possible workaround: "line = line.split('\t')[0]" after line "for line in kos:". But - is it really desired?
    """Reads and interprets data of taxonomic relations from files in Graphlan-like format.

    Args:
        input_d (dict of str: str): A dict, where:
            keys are names of directories,
            values are lists containing filenames that are inside the directory

        Example of featured file:
            Listeriaceae.Listeria.Lgrayi
            Listeriaceae.Listeria.Linnocua

    Returns:
        A tuple (total_tax_tree, per_file_tax_tree, multi_flat_tax_tree):

        total_tax_tree: A dict of dicts, etc. In form of nested dicts, represents
            phylogenetic trees of entities described by files featured in the input.
            This is a total tree - nodes from all files are merged here.
            Example:
            {
                'Bacillaceae':
                {
                    'Bacillaceae_____Anoxybacillus': {}
                },
                'Listeriaceae':
                {
                    'Listeriaceae_____Listeria':
                    {
                        'Listeriaceae_____Listeria_____Lgrayi': {},
                        'Listeriaceae_____Listeria_____Linnocua': {},
                    }
                }
            }

        per_file_tax_tree: A dict of dicts, etc. Keys are filenames.
            In form of nested dicts, represents phylogenetic trees of entities,
            described by files featured in the input.
           Example:
            {
                'annot_1.txt':
                {
                    'Bacillaceae':
                    {
                        'Bacillaceae_____Anoxybacillus': {}
                    }
                }
                'annot_2.txt':
                {
                    'Listeriaceae':
                    {
                        'Listeriaceae_____Listeria':
                        {
                            'Listeriaceae_____Listeria_____Lgrayi': {},
                            'Listeriaceae_____Listeria_____Linnocua': {},
                        }
                    }
                }
            }

        multi_flat_tax_tree: A dict of dicts.
            First level: keys are filenames, values are dicts.
            Second level: keys are names of nodes,
                          values are numbers of nodes with names starting from "node_name"
            Example:
            {
                'annot_1.txt':
                {
                    'Bacillaceae': 8,
                    'Bacillaceae_____Anoxybacillus': 8
                },
                'annot_2.txt':
                {
                    'Listeriaceae': 16,
                    'Listeriaceae_____Listeria': 16,
                    'Listeriaceae_____Listeria_____Lgrayi': 8,
                    'Listeriaceae_____Listeria_____Linnocua': 8
                 }
            }
    """
    total_tax_tree = {}
    per_file_tax_tree = {}
    multi_flat_tax_tree = {}
    for path in input_d:
        for plik in input_d[path]:
            multi_flat_tax_tree[plik] = {}
            kos = open(pjoin(path, plik), 'r').readlines()
            file_tax_tree = {}
            for line in kos:
                line = line.rstrip()
                linia = line.split('.')
                linia = linia_unique(linia)
                file_tax_tree = tax_tree_extend(file_tax_tree, linia)
                for node in linia:
                    try:
                        multi_flat_tax_tree[plik][node] += 1
                    except KeyError:
                        multi_flat_tax_tree[plik][node] = 1
            total_tax_tree.update(file_tax_tree)
            per_file_tax_tree[plik] = file_tax_tree
    return total_tax_tree, per_file_tax_tree, multi_flat_tax_tree

def xml_counts_graphlan(tax_tree, per_file_tax_tree, xml_names, multi_flat_tax_tree):
    # TODO: tax_tree is not used here
    """Groups and/or sums numbers of nodes (obtained from multi_flat_tax_tree):
        1. by node_name - grouping
        2. by identificator (obtained from filename) - summing

    Args:
        tax_tree: A dict of dicts, etc. Taxonomic tree represented by nested dicts.
            Example: {a:{aa:{}, ab:{aba:{}, abb:{abba:{}}}}}
            Explicit example in description of tax_tree_graphlan(), look for: total_tax_tree.

        per_file_tax_tree: A dict of dicts, etc. Top-level keys are filenames.
            Taxonomic tree represented by nested dicts. Similar to tax_tree.
            Example included in description of tax_tree_graphlan().

        xml_names: identificators, obtained by modifying filenames from input_d dict.
            Example: [identificator_1, identificator_2, ...]

        multi_flat_tax_tree: A dict of dicts.
            First level: keys are filenames, values are dicts.
            Second level: keys are names of nodes,
                          values are numbers of nodes with names starting from "node_name"
            Explicit example in description of tax_tree_graphlan().

    Returns:
        A tuple (xml_dict, name_total_count):

        xml_dict: A dict. Contains lists of numbers of occurrences of nodes by node. Keys are names of nodes.
            Values are lists of constant length, where every item on position *X* means:
            number of occurrences of nodes with prefix equal to *name of node* in file *X*.
            *file X* is that file, which is on position *X* on list xml_names.
            Example: {node: [count_1, count_2], node_2: [count_1, count_2]}

        name_total_count: A dict. Contains summed numbers of occurrences of nodes by file.
            Keys are identificators (derived from filenames - look for xml_names), values ale sums.
            Example: {identificator_1: count_1, identificator_2: count_1}
    """
    all_nodes = set()
    xml_dict = {}
    name_total_count = {}
    for plik in multi_flat_tax_tree:
        all_nodes = all_nodes | set(multi_flat_tax_tree[plik].keys())
    for node in all_nodes:
        xml_dict[node] = []
        for name in xml_names:
            for plik in multi_flat_tax_tree:
                if name in plik:
                    if node in multi_flat_tax_tree[plik]:
                        xml_dict[node].append(multi_flat_tax_tree[plik][node])
                    else:
                        xml_dict[node].append(0)
    for name in xml_names:
        for plik in multi_flat_tax_tree:
            if name in plik:
                name_total_count[name] = 0
                for node in per_file_tax_tree[plik]:
                    name_total_count[name] += multi_flat_tax_tree[plik][node]
    return xml_dict, name_total_count


def graphlan_to_krona(input_d):
    """
    Modify files (specified in inupt_d) created by Graphlan to create
    a set of xml strings, which may be assembled into a input file for Krona.

    Args:
        input_d (dict of str: str): A dict, where:
            keys are names of directories,
            values are lists containing filenames,
            which are inside the directory

    Returns:
        A tuple (xml_names, xml_dict, tax_tree, name_total_count):

        xml_names: identificators, obtained by modifying filenames,
            given by input_d dict.
            Example:
            [identificator_1, identificator_2, ..., identificator_x]

        xml_dict: A dict. Contains lists of numbers of occurrences of nodes by node. Keys are names of nodes.
            Values are lists of constant length, where every item on position *X* means:
            number of occurrences of nodes with prefix equal to *name of node* in file *X*.
            *file X* is that file, which is on position *X* on list xml_names.
            Note, that for all i: len(xml_dict[i]) is equal to len(xml_names).
            Example:
            {node: [count_1, count_2, ... count_x], node_2: [count_1, count_2, ... count_x]}

        tax_tree: A dict of dicts, etc. In form of nested dicts, it represents phylogenetic tree.
            Example:
            {a:{aa:{}, ab:{aba:{}, abb:{abba:{}}}}}
            Explicit example in description of tax_tree_graphlan, look for: total_tax_tree.

        name_total_count: A dict. Contains summed numbers of occurrences of nodes by file.
            Keys are identificators (derived from filenames - look for xml_names), values ale sums.
            Example:
            {identificator_1: count_1, identificator_2: count_2}

    """
    input_d = txt_dict_clean(input_d)
    xml_names = xml_names_graphlan(input_d)
    tax_tree, per_file_tax_tree, multi_flat_tax_tree = tax_tree_graphlan(input_d)
    xml_dict, name_total_count = xml_counts_graphlan(tax_tree, per_file_tax_tree, xml_names, multi_flat_tax_tree)
    return xml_names, xml_dict, tax_tree, name_total_count


def aftershave(opts):
    """TODO

    Args:
        opts: namespace with options explained in file "bipype"
    """
    SSU = {}
    metag_flag = 0
    if '16S' in opts.output_type:
        SSU['16S'] = SSU_read('/home/pszczesny/soft/bipype/SSU_candidate_db.fasta', '16S')
        metag_flag = 1
    if 'ITS' in opts.output_type:
        SSU['ITS'] = SSU_read('/home/pszczesny/soft/bipype/UNITE_public_from_27.01.13.fasta')
        metag_flag = 1
    input_dic = input_locations(opts.mode, opts.output_type)
    analysed_dict = {}
    tax_dict = {}
    pure_tax = {}
    for out_type in opts.output_type:
        try:
            analysed_dicto, tax_dicto = dict_prepare(out_type, input_dic[out_type], SSU[out_type])
            analysed_dict[out_type] = analysed_dicto
            tax_dict[out_type] = tax_dicto
            pure_tax.update(tax_dicto)
        except KeyError:
            if len(input_dic) == 1:
                input_dicto = input_dic['txt']
    if metag_flag == 1:
        krona_xml_name, krona_html_name = out_namespace(opts.out_dir, opts.output_type)
        xml_names, xml_dict, tax_tree, name_total_count = xml_format(analysed_dict, pure_tax)
        krona_unit = 'reads'
    else:
        if opts.output_type[0] == 'txt':
            krona_xml_name, krona_html_name = out_namespace(opts.out_dir, opts.output_type[0])
            xml_names, xml_dict, tax_tree, name_total_count = graphlan_to_krona(input_dicto)
            krona_unit = 'processes'
    xml_string = xml_prepare(xml_names, xml_dict, tax_tree, name_total_count, krona_unit)
    outprint(xml_string, krona_xml_name)
    krona_to_html_comm = 'ktImportXML -o %s %s'%(krona_html_name, krona_xml_name)
    if opts.mode == 'run':
        system(krona_to_html_comm)
