#!/usr/bin/python

from os import getcwd, walk, system, chdir, stat
from os.path import join as pjoin
from os.path import exists as pexists
from string import split, find, join
from itertools import combinations
from datetime import datetime
from copy import deepcopy
import cPickle as pickle
from settings_bipype import *


def cat_read(mode, fileext, paired_end=True):
    """Returns a dict with paths to sequence files
    from current working directory.

    The dict has following format:
    {directory_path1:[file1_in_this_directory_path,
                      file2_in_this_directory_path,
                      file3_in_this_directory_path],
     directory_path2:[file4_in_this_directory_path,
                      file5_in_this_directory_path,
                      file6_in_this_directory_path]
    }
    Paths to files (values of dict) are relative.

    Args:
        mode:       If mode is 'run', the gunzip command will be
                    executed to extract compressed files.
        fileext:    Extension of sequence files, which will be putted into
                    the dict.
        paired_end: If True, only paired-end reads are included
                    in the dict (default True).
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
    if paired_end:
        seq_dict = paired_end_match(seq_dict)
    else:
        pass
    return seq_dict


def exist_check(program, names, todo):
    """Checks if a part of program is actually done.

    Function tries to find output files from different parts of program
    and basing on the results, removes corresponding tokens from
    todo list.

    In one special case (program=='usearch_0') function checks size of
    the file specified by <names> argument in bits. If size==0: todo=[]
    and PANIC MODE information is printed.

    Args:
        program:  One of the following strings representing programs
                  (de facto functions):
                        'refseq'            (&)
                        'usearch'           (@)
                        'usearch_0'         (@)
                        'MV'                (^)
                        'rapsearch'         (@)
                        'cutadapt'          (^)


        names:    There are three possibilities:
                    (&) A dict in following format:
        {file extension : path to file with corresponding extension}
                    (@) A string (one path)
                    (^) A list of paths
                  Paths are hypothetical - if one really exists,
                  appropriate token is removed from todo list.

        todo:     List of tokens representing different parts of
                  function specified by 'program' argument.

    Please, pay attention to mutual compatibility of arguments.
    For more information please refer to adequate one of the following
    functions code:
        - refseq_ref_namespace() & refseq_mapping()
        - usearch()
        - MV()
        - rapsearch()
        - cutadapt()

    Result:
        todo:     Checked list of tokens.
    """
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
            todo = []
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
    """Returns a dict with {GI:TaxID} from file.

    Keys and values are integers.

    File has following format:
    13  9913
    15  9915
    16  9771
    17  9771
    where first column is a GI, second is a TaxID.

    GLOBLAS:
        PATH_GI_TAX
    """
    print 'reading tax ids'
    gi_tax_dict = {}
    #curr_t = datetime.now()
    with open(PATH_GI_TAX) as file_:
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
    """Returns a dict with {TaxID:scientific_name} from file.

    Lines without scientific names are omitted.

    Keys are integers, values are strings.

    Example (appropriate file format included):

        File:
    2   |       prokaryotes     |       prokaryotes <Bacteria>  |       in-part |
    6   |       Azorhizobium    |               |       scientific name |
    6   |       Azorhizobium Dreyfus et al. 1988        |               |       synonym |
    6   |       Azotirhizobium  |               |       equivalent name |
    7   |       ATCC 43989      |               |       type material   |
    7   |       Azorhizobium caulinodans        |               |       scientific name |

        Output:
            {6:'Azotirhizobium', 7:'Azorhizobium caulinodans'}

    GLOBALS:
        PATH_TAX_NAME
    """
    print 'reading tax names'
    tax_name_dict = {}
    #curr_t1 = datetime.now()
    with open(PATH_TAX_NAME) as file_:
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
        a dict with {GI:#_mapped_reads} (keys and values are integers).
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
    """Parses and writes data from samtools idxstats output file.

    Firstly, using idx_reader(file_), creates {GI:#_mapped_reads}
    (a dict).
    Secondly, replaces every GI (key) with TaxID if appropriate
    one is available in tax_id_dict.
    Than, replaces every GI /TaxID (key) with scientific name if
    appropriate one is available in tax_name_dict.
    Finally, writes data (sorted by numbers of mapped reads in
    decreasing order) to outfile in key;value format, where:
        key   is GI/TaxID/scientific name
        value is number of mapped reads

    Args:
        mode:          If (mode == 'run'), function do mentioned things
                       Elif (mode == 'test') function prints
                            file_ and outfile arguments.
                       Else function do nothing.

        file_:         Path to samtools idxstats output file.
                       For more information refer to idx_reader()

        tax_name_dict: {TaxID:scientific_name} dict.
                       For more information refer to tax_name_reader()

        tax_id_dict:   {GI:TaxID} dict,
                       For more information refer to tax_id_reader()

        outfile:       Path to output file.
    """
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
    """Returns one name for given tuple containing pair of filenames.

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


def refseq_ref_namespace(directory, seq, postfix, out_dir='in_situ', map_dir='in_situ'):
    """Returns a dict within:
        - keys are following file extensions: fastq, sam, sam2, bam,
          sorted, sorted.bam, idx_stats, tax_count, map_count
        - values are paths to file with corresponding extension.
        Filenames have following format:
        sample_name + '_' + postfix + extension where sample_name is
        basename of seq (if seq is a file) or
        output of pair_uni_name(seq) (if seq is a tuple)

    Args:

        directory: Path to directory, where file with .fastq extension
                   will be written

        seq:       Tuple of paired_end read or sequence file.

        postfix:   String added to the end of file basenames.

        out_dir:   Path to directory, where files with .sam, .sam2,
                   .bam,.sorted, .sorted.bam, .idx_stats and .tax_count
                   extensions will be written.
                   If out_dir is 'in_situ' (default), out_dir=directory

        map_dir:   Path to directory, where file with
                   .map_count extension will be written. If map_dir is
                   'in_situ'(default),map_dir=out_dir.
"""
    if out_dir == 'in_situ':
        out_dir = directory
    ref_namespace = {}
    if type(seq) == tuple:
        sample_name = pair_uni_name(seq)
        ref_namespace['fastq'] = (pjoin(directory, seq[0]), pjoin(directory, seq[1]))
    else:
        ref_namespace['fastq'] = pjoin(directory, seq)
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


def fastq_dict(seq_dict, root, file_):  # Why not defaultdict from collections?
    """Adds file_ to list available under seq_dict[root]

    If key 'root' is not present, function will create it.
    """
    if root in seq_dict:
        seq_dict[root].append(file_)
    else:
        seq_dict[root] = [file_]
    return seq_dict


def paired_end_match(seq_dict):
    """Returns a dict with paths to paired-end reads only.
    Files will be matched as pair if their names differs only by parts:
        - R1 and R2, if the parts of filename are separated by underscores ('_')
        - 1 and 2, if the parts of filename are separated by dots ('.')

    Arg:
        seq_dict: {directory_path:file_path}

    Returns:
        Dict in following format:
    {directory_path1:[(paired-end_read1_R1_path, paired-end_read1_R2_path),
                      (paired-end_read2_R1_path, paired-end_read2_R2_path),
                      (paired-end_read3_R1_path, paired-end_read3_R2_path)],
     directory_path2:[(paired-end_read4_R1_path, paired-end_read4_R2_path)]
    }

    Paths to files (both in input and output dict) are relative.
    """
    pe_dict = {}
    for directory in seq_dict:
        for seq1, seq2 in combinations(seq_dict[directory], 2):
            # for format with underscores
            seq1_u = set(split(seq1, '_'))
            seq2_u = set(split(seq2, '_'))
            # for format with dots, like output of metAMOS preprocessing
            seq1_d = set(split(seq1, '_'))
            seq2_d = set(split(seq2, '_'))
            if seq1_u ^ seq2_u == {'R1', 'R2'} or seq1_d ^ seq2_d == {'1', '2'}:
                if directory in pe_dict:
                    pe_dict[directory].append((seq1, seq2))
                else:
                    pe_dict[directory] = [(seq1, seq2)]
    return pe_dict


def bowtie2_run(mode, proc, ref, out, inp1, inp2=False):
    """If mode=='run', launches Bowtie 2.

    Args:
        mode

        proc: [-p in Bowtie 2] 'Number of alignment threads to launch.'

        ref:  [-x in Bowtie 2] 'The basename of the index
                                for the reference genome.'

        out:  [-S in Bowtie 2] 'File to write SAM alignments to.'

        inp1: [-1 in Bowtie 2] 'Comma-separated list of files
                                containing mate 1s (filename usually
                                includes _1), e.g.
                                -1 flyA_1.fq,flyB_1.fq  '

        inp2: [-2 in Bowtie 2]  'Comma-separated list of files
                                containing mate 2s (filename usually
                                includes _2), e.g.
                                -2 flyA_2.fq,flyB_2.fq  '
              If inp2=False: -2 argument is not given in Bowtie 2.
                     (False by default)
              Else: inp2 string is passed to Bowtie 2 as -2 argument.
    """
    if inp2:
        print 'bowtie2 -p %i -x %s -1 %s -2 %s -S %s'%(proc, ref, inp1, inp2, out)
        if mode == 'run':
            system(PATH_BOWTIE2+' -p %i -x %s -1 %s -2 %s -S %s'%(proc, ref, inp1, inp2, out))
    else:
        print 'bowtie2 -p %i -x %s -f -U %s -S %s'%(proc, ref, inp1, out)
        if mode == 'run':
            system(PATH_BOWTIE2+' -p %i -x %s -f -U %s -S %s'%(proc, ref, inp1, out))


def sam_merge(sam1, sam2):
    """Merges two files in SAM format into one.

    Args:
        sam1: SAM file
        sam2: SAM file

    Returns:
       sam1 file in SAM format which contains now lines from both files.
       Function also removes sam2 file.
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
            f1c_idx += 1

    f2c_idx = 0
    for line in f2c:
        line = line.rstrip()
        if line[0] == '@':
            f2c_idx += 1

    common = f1c[0:f1c_idx-1] + f2c[1:f2c_idx-1] + [f1c[f1c_idx-1]]
    cont = f1c[f1c_idx:]+ f2c[f2c_idx:]

    f3 = open(sam1, 'w')
    f3.write(join(common, ''),)
    f3.write(join(cont, ''),)
    f3.close()
    remove_sam2 = 'rm %s' % (sam2)
    system(remove_sam2)


def bam_make(mode, sam, bam): # multi
    """If mode=="run" and @SQ lines are present in the header,
    function converts SAM to BAM

    Args:
        mode: string type
        sam: SAM file
        bam: BAM file
    """

    print 'samtools view -bS %s > %s'%(sam, bam)
    if mode == 'run':
        system('samtools view -bS %s > %s'%(sam, bam))


def bam_sorting(mode, bam, sorted_name): # multi
    """If mode=="run" function reads bam (file in BAM format), sort it by aligned
    read position and write it out to BAM file whose name is: <sorted_name>.

    Args:
        mode: string type
        bam: BAM file
        sorted_name: name of output file
    """
    print 'samtools sort %s %s'%(bam, sorted_name)
    if mode == 'run':
        system('samtools sort %s %s'%(bam, sorted_name))


def bam_indexing(mode, sorted_bam): # multi
    """If mode=="run" creates an index file
    <sorted_bam>.bam.bai for the <sorted_bam>.bam file.

    Args:
        mode: string type
        sorted_bam: BAM file
"""
    print 'samtools index %s'%(sorted_bam)
    if mode == 'run':
        system('samtools index %s'%(sorted_bam))


def bam_idxstating(mode, sorted_bam, idxstats): # multi
    """Launch samtools idxstats.

    Args:
        mode:       if mode=='run', launch program.
        sorted_bam: path to file in BAM format, which is the input
                    to samtools.
        idxstats:   path to file, where the output from samtools will be
                    written.
    """
    print 'samtools idxstats %s > %s'%(sorted_bam, idxstats)
    if mode == 'run':
        system ('samtools idxstats %s > %s'%(sorted_bam, idxstats))


def idxstat_perling(mode, idxstats, map_count): # multi
    """Counts and writes to map_count both:
       - (sum of all (numbers of mapped reads) in <idxstats>),
       - (sum of all (numbers of unmapped reads) in <idxstats>).

    Function launches short "one-liner" in perl.

    Output has following format:
                    #mapped - #unmapped
    For example:        123 - 456

    Args:
        mode:      If mode=='run' "one-liner" is launched.
        idxstats:  Path to samtools idxstats output file.
                   Please refer to idx_reader() for more information.
        map_count: Path to file, where differences will be written.
    """
    perl_command = """perl -e 'while(<>){chomp;@a=split "\t", $_; $b+=$a[2]; $c+=$a[3];} print "$b - $c\n";' %s > %s"""%(
        idxstats, map_count
        )
    print perl_command
    if mode == 'run':
        system(perl_command)


def refseq_mapping(mode, e, directory, pair, postfix, refseq, tax_name_dict, tax_id_dict, threads, map_dir, refseq_2=False):
    """Aligns reads to reference sequence(s) using Bowtie 2,
       than parses and writes data to files.

    Firstly, Bowtie 2 align reads to reference sequence (or sequences,
        if refseq_2 is not False and merge output SAM files).
    Secondly, BAM files are made, sorted and indexed.
    Thirdly, Function launches 'samtools idxstats'.
    In the next step, perl 'one-liner' counts and writes to file
        sums of mapped reads and unmapped reads.
    Finally, idx_map() function is called:
        Parse data from samtools idxstats output file and
        writes data to outfile in key; value format, where:
        key     is    GI/TaxID/scientific name
        value   is    number of mapped reads.

    Args:
        mode:            If mode=='run', function operate on data
                         and print information about it.
                         Else, prints information, without operating
                         on data (it is kind of test).

        e:               If True, checks if a part of workflow is
                         actually done and doesn't duplicate this jobs.
                         For more information, please refer to
                         exist_check() function.

        directory:       Path to directory, where files will be written.

        pair:            Tuple of paths to file (paired-end reads) OR
                         String: path to sequence file.

        postfix:         String added to the end of files basenames.
                         Argument for refseq_ref_namespace().

        refseq:          'The basename of the index for the reference
                         genome.' Argument for bowtie2_run().

        tax_name_dict:   {TaxID:scientific_name} (a dict).
                         For more information refer to tax_name_reader()
                         Argument for idx_map().

        tax_id_dict:     {GI:TaxID} (a dict),
                         For more information refer to tax_id_reader()
                         Argument for idx_map().

        threads          Number of threads for Bowtie 2 calculations.

        map_dir          Directory where sums of mapped and unmapped
                         reads will be written. For more information
                         refer to idxstat_perling().

        refseq_2:        'The basename of the index for the reference
                         genome.' Argument for bowtie2_run().
                         If selected, launches Bowtie 2 on it, then
                         merge output with output from Bowtie 2
                         launched on 'refseq'.
                         If refseq_2==False, Bowtie 2 is working only
                         on 'refseq' argument.
                         (False by default).

    For more information take a look at:  refseq_ref_namespace()
                                          exist_check()
                                          bowtie2_run()
                                          sam_merge()
                                          bam_make()
                                          bam_sorting()
                                          bam_indexing()
                                          bam_idxstating()
                                          idxstat_perling()
                                          idx_map()
    """
    todo = ['bowtie', 'bam_make', 'sort_index', 'idxstat', 'perl', 'idx_map']
    ref_namespace = refseq_ref_namespace(directory, pair, postfix, 'in_situ', map_dir)
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
    """Returns 200 in case, when read length
    is less then 200 or 500 in other case.

    Args:
        pair: tuple of paired_end read
        cat: name of folder with the sample files

    Example:
        input: ('Amp15_BFk_B_p_CGTACG_L001_R1_001.fastq','Amp15_BFk_B_p_CGTACG_L001_R2_001.fastq'), 'catalog_name'
        output: 500
    """

    sample_name = pair[0]
    zeros = find(sample_name, '00') #WARNING HARDCODE
    try:
        ins_len = int(sample_name[zeros-1:zeros+2])
    # CHANGED
    except ValueError:
        aqq = len(open(pjoin(cat, sample_name), 'r').readlines()[1].rstrip())
        if aqq > 200:
            ins_len = 500
        else:
            ins_len = 200
    return ins_len


def gzip_MV(MV_dir):
    """
    Returns archive with catalogs Sequences, Roadmaps, PreGraph, Graph2,
    LastGraph and files with statistics.

    Args:
        MV_dir: folder with Velvet results

"""
    gzip_list = [
        'Sequences', 'Roadmaps', 'PreGraph', 'Graph2', 'contigs.fa', 'stats.txt', 'LastGraph',
        'meta-velvetg.Graph2-stats.txt', 'meta-velvetg.LastGraph', 'meta-velvetg.LastGraph-stats.txt',
        'meta-velvetg.split-stats.txt'
        ]
    to_gzip = ''
    for element in gzip_list:
        to_gzip += ' %s'%(pjoin(MV_dir, element))
        to_gzip_comm = 'tar --remove-files -czf %s %s'%(pjoin(MV_dir, 'intermediates.tgz'), to_gzip)
        system(to_gzip_comm)


def rapsearch(mode, e, contig_loc, rap_out, KEGG=None):
    """Runs RAPSearch with using KEGG databases for similarity search.

    Args:
        mode: if mode="run", then program runs rapsearch
        e: if e=True, then function runs exist_check function
        contig_loc: a query file
        rap_out: output file name
        KEGG: default is None, if KEGG= KO,
            then ko.pep.rapsearch.db is chosen as protein database

    HARDCODED:
        if KEGG='masl'

    GLOBALS:
        - path to RAPSearch program: PATH_RAPSEARCH
        - paths to data for RAPSearch:
            PATH_REF_PROT_MASL28,
            PATH_REF_PROT_KO,
            PATH_REF_PROT_ELSE
    """
    rap_log = rap_out + '.log'
    rap_err = rap_out + '.err'
    if KEGG == 'masl':
        ref_prot_loc = PATH_REF_PROT_MASL28
        rap_com = '%s -q %s -d %s -o %s -z 12 -v 20 -b 1 -t n -a t 1> %s 2> %s'%(
            PATH_RAPSEARCH, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    elif KEGG == 'KO':
        ref_prot_loc = PATH_REF_PROT_KO
        rap_com = '%s -q %s -d %s -o %s -z 12 -v 20 -b 1 -t n -a t 1> %s 2> %s'%(
            PATH_RAPSEARCH, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    else:
        ref_prot_loc = PATH_REF_PROT_ELSE
        rap_com = '%s -q %s -d %s  -o %s -z 10 -e 0.001 -b 100 -v 100 -g T -a T 1> %s 2> %s'%(
            PATH_RAPSEARCH, contig_loc, ref_prot_loc, rap_out, rap_log, rap_err
            )
    print rap_com
    todo = [rap_com]
    if e:
        todo = exist_check('rapsearch', rap_out, todo)
    if (mode == 'run') and (len(todo) != 0):
        system(rap_com)


def MV(mode, e, k_mers, cat, pair, ins_len, rap=False):
    """Runs Velvet. Runs function gzip_MV on folder with Velvet results,
        which name was created with command pair_uni_name(pair) + '_velvh_out'
    If parameter rap is true, then the function also run function rapsearch.

    Args:
        mode: if mode=="run" program runs velveth or velvetg or metavelvetg
        e: boolean parameter. If e==True, then program changes todo list with exist_check function
        k_mers: k-length nucleotides reads list
        cat: name of current folder
        pair: tuple of paired_end reads
        ins_len: If ins_len==9999, then ins_len is the output from ins_len_read function
        rap: If rap==True, then program runs rapsearch function. Default rap=False

    GLOBALS:
        PATH_VELVETH
        PATH_VELVETG
        PATH_METAVELVETG
    """
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
        velveth_run = PATH_VELVETH+' %s %i -%s -shortPaired %s %s'%(
            tmp_out_dir, idx, fileext, pjoin(cat, pair[0]), pjoin(cat, pair[1])
            )
        velvetg_run = PATH_VELVETG+' %s -exp_cov auto -ins_length %i'%(tmp_out_dir, ins_len)
        meta_run = PATH_METAVELVETG+' %s -ins_length %i | tee %s'%(tmp_out_dir, ins_len, log_loc)
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


def MH(mode, e, t, cat, pair, rap=False, presets='meta-large'):
    """Runs Megahit. If rap=True, function also run rapsearch().

    Args:
        mode:    if mode!="run" program prints commands without running them
        e:       boolean parameter. If e==True, then program changes todo list
                   with exist_check function. NOT YET !!!!!!!!!!!
        t:       number of threads
        cat:     name of current folder
        pair:    tuple of paired_end reads
        rap:     if rap==True (default False), program runs rapsearch function.
        presets: argument for megahit (default 'meta-large')

    GLOBAL:
        PATH_MEGAHIT
    """
    out_dir = pjoin(getcwd(), pair_uni_name(pair) + '_megahit_out')
    mega_hit_run = (PATH_MEGAHIT + ' -t ' + str(t) + ' -1 ' + pair[0] + ' -2 ' +
        pair[1] + ' -o ' + out_dir + ' --presets ' + presets)
    print mega_hit_run
    if mode=="run":
        system(mega_hit_run)
        system('rm -r ' + pjoin(out_dir,'intermediate_contigs'))
    if rap:
        rap_in = pjoin(out_dir, 'final.contigs.fa')
        rap_f =  pair_uni_name(pair) + '.rapsearch'
        rap_out = pjoin(out_dir, rap_f)
        rapsearch(mode, e, rap_in,  rap_out)


def usearch(mode, e, search_type, infile, database, threads):
    """Launches Usearch with -usearch_local command.

    HARDCODED:
        - Usearch options:
            - evalue = 0.01
            - id = 0.9
            - strand = both
    Args:
        mode:        If mode=='run' Usearch is launched.
        e:           If True, checks if a part of workflow is
                        actually done and doesn't duplicate this
                        jobs. For more information, please
                        refer to exist_check() function.
        search_type &
        infile:      outfile = infile+'.usearch_'+search_type
                        outfile is passed to [-blast6out in Usearch]
        database:    [-db in Usearch]
        threads:     [-threads in Usearch] Number of threads used in
                     calculations.

    GLOBALS:
        PATH_USEARCH
    """
    outfile = infile + '.usearch_' + search_type
    usearch_command = '%s -usearch_local %s -db %s -evalue 0.01 -id 0.9 -blast6out %s -strand both -threads %i'%(
        PATH_USEARCH, infile, database, outfile, threads
        )
    todo = ['usearch']
    print usearch_command
    if e:
        todo = exist_check('usearch', outfile, todo)
    todo = exist_check('usearch_0', infile, todo)
    if (mode == 'run') and (len(todo) != 0):
        system(usearch_command)


def adapter_read_bck(adapter_file, filename):
    """Function reads all lines of adapter_file and separates words in each line.
    If first word in a line is equal to filename,
    then second and third words are returned as a tuple.

    Args:
        adapter_file, filename

    Returns:
        A tuple (fin_adap_1, fin_adap_2)

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
    """Replace 'NNNNNN' part of adp_2 by this part of the filename,
     which contains only letters that are symbols of nucleotides.

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
    """Searches for the adapters in reads from input files, removes them,
    when it finds them and writes to output files which have .cutadapt.fastq
    extension.
    Then function runs FLASH (software tool to merge paired-end reads)
    with .cutadapt.fastq
    files as input and .amplicons.cutadapt.flash.merged.fastq files as output.
    These results are input for fastq_to_fasta which converts file format from
    fastq to fasta.

    Args:
        mode:             if mode=='run', function operate on data
        e:                if True, checks if a part of the workflow is actually
                          done
                          and omits these parts, to avoid duplicating this job.
        cat:              name of folder with R_1 and R_2 files
        R1_file, R2_file: input files
        adapter_file:     If(adapter_file='use_filenames'), function will use
                          adapters returned by adapter_read(R1_file)
                          Else,function will use adapters returned by
                          adapter_read_bck(adapter_file, R1_file)
        usearch_16S:      if True runs usearch(mode, e, '16S',
                          outname_uni_fasta,   usearch_16S, threads)
                          where outname_uni_fasta is cutadapt output.
        usearch_ITS:      if True runs usearch(mode, e, 'ITS',
                          outname_uni_fasta, usearch_ITS, threads)
                          where outname_uni_fasta is cutadapt output.
     For more information please refer to:
        - adapter_read_bck()
        - adapter_read()

    GLOBALS:
        PATH_FQ2FA
    """
    R1_fastq = pjoin(cat, R1_file)
    R2_fastq = pjoin(cat, R2_file)
    if adapter_file == 'use_filenames': #WARNING HARDCODE
        adapter_1, adapter_2 = adapter_read(R1_file)
    elif adapter_file == 'use_paths':
        adapter_1, adapter_2 = adapter_read(R1_fastq)
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
    fq2fa_com = '%s < %s > %s -Q33'%(PATH_FQ2FA, outname_uni_fastq_postflash, outname_uni_fasta)
    print fq2fa_com
    if 'fq2fa' in todo and mode == 'run':
        system(fq2fa_com)
    if usearch_16S:
        usearch(mode, e, '16S', outname_uni_fasta, usearch_16S, threads)
    if usearch_ITS:
        usearch(mode, e, 'ITS', outname_uni_fasta, usearch_ITS, threads)


def auto_tax_read(db_loc):
    """Reads {GI:TaxID} & {TaxID:scientific_name} dicts from pickle

    Args:
        db_loc: path to pickle file

    Returns:
        Two dicts:
          {GI:TaxID}              Keys are integers, values are strings
          {TaxID:scientific_name} Keys and values are integers
    """
    with open(db_loc, 'rb') as fp:
        tax_names = pickle.load(fp)
        tax_id = pickle.load(fp)
    return tax_id, tax_names


def taxa_read(read_mode, db_loc=None):
    """Returns {GI:TaxID} & {TaxID:scientific_name} dicts.

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
        Two dicts:
          {GI:TaxID}              Keys are integers, values are strings
          {TaxID:scientific_name} Keys and values are integers
    """
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
    """Runs bwa. Find the SA coordinates of the input reads.
    Generate alignments in the SAM format given single-end reads.
    Repetitive hits will be randomly chosen.

    Args:
        mode: if mode="run", then commands "bwa aln" and "bwa samse" were run
        thr: Number of threads (multi-threading mode)
        e: (!!!! Wasn`t used !!!!)
        pair: tuple of paired_end read
        cat: name of folder with the sample files
        prefix: prefix for output files names
        rec_db_loc: input database with fasta files
    """
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
    """
Copies humann to the current directory, moves input (*.m8) files to the input directory,
copies hmp_metadata.dat file to the input directory
 and runs humann

GLOBALS:
        - path to humann program: PATH_HUMANN
        - path to data for humann: PATH_HUMANN_DATA

Args:
    mode: if mode="run", then humann will be copied to the current directory
    e: if e=True, then function checks existion of humann-0.99 results folder
    m8_dict: the similarty search results folders
    typ: default typ="m8", in that case new catalog humann-0.99 will be created in rapsearch result folder;
         in other case humann analysis results will be added in rapsearch result folder
"""
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
            get_humann = 'cp -r %s %s'%(PATH_HUMANN, path)
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
            hmp_get_comm = 'cp %s %s/'%(PATH_HUMANN_DATA, input_loc)
        else:
            hmp_get_comm = 'cp %s %s/'%(PATH_HUMANN_DATA, hum_loc)
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
    """Runs programs such, as Velvet, bwa, Bowtie 2, humann,
    cuts adapters or run RAPSearch if these options were chosen.

    Args:
        opts: A namespace, where:

            opts.to_calculate:
                Input sequences type (Plant or fungi)

            opts.mode:
                A mode in which the program will be run.

                One of: ['test', 'run']

            opts.db_NCBI_taxonomy:
                Path to pickle file

            opts.reconstruct:
                Boolean attribute, which allows to run reconstruct function if True.

            opts.assembler:
                k-length nucleotids reads list or None

            opts.db_reconstruct:
                input database with fasta files

            opts.threads:
                Number of threads (multi-threading mode)

            opts.e:
                Boolean attribute, which checks file existing

            opts.ins_len:
                200 in case, when read length is less,
                then 200 and 500 in other case.

            opts.db_refseq_fungi:
                A list of paths to refseq databases to use in fungi analysis.
                Up to two paths are allowed.
                Warning:
                If there are multiple databases with filenames like
                <your_path><second_part_of_name>, all will be loaded.

            opts.db_refseq_plant:
                A list of paths to refseq databases to use in plants analysis.
                Up to two paths are allowed.
                Also check warning in opts.db_refseq_fungi description.

            opts.our_dir:
                Directory where sums of mapped and unmapped
                reads will be written

            opts.postfix:


            opts.cutadapt:
                List of adapter types: 16S, ITS or 'both'.

            opts.db_16S:
                boolean parameter; if True, then Usearch is runned for 16S adapters
                

            opts.db_ITS:
                boolean parameter; if True, then Usearch is runned for ITS adapters

    GLOBALS:
        PATH_FQ2FA
    """
    if len(opts.db_refseq_plant) == 1:
        opts.db_refseq_plant += [False]
    if len(opts.db_refseq_fungi) == 1:
        opts.db_refseq_fungi += [False]

    if opts.to_calculate == None:
        opts.to_calculate = []
    refseq_test = len(set(['f', 'b', 'p']) & set(opts.to_calculate))
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
            if opts.assembler == ["MH"]:
                MH(opts.mode, opts.e, opts.threads, cat, pair, ('rap_prot' in opts.to_calculate))
            elif opts.assembler != None:
                MV(opts.mode, opts.e, opts.assembler, cat, pair, opts.ins_len, ('rap_prot' in opts.to_calculate))
            if ('f' in opts.to_calculate) or ('b' in opts.to_calculate):
                postfix = opts.postfix + 'fungi'
                refseq_mapping(
                    opts.mode, opts.e, cat, pair, postfix, opts.db_refseq_plant[0], tax_name_dict,
                    tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq_plant[1]
                    )
            if ('p' in opts.to_calculate) or ('b' in opts.to_calculate):
                postfix = opts.postfix + 'plant'
                refseq_mapping(
                    opts.mode, opts.e, cat, pair, postfix, opts.db_refseq_fungi[0], tax_name_dict,
                    tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq_fungi[1]
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

    if opts.assembler != None:
        contig_dict = cat_read(opts.mode, 'fa', False)
        for cat in contig_dict:
            for contig in contig_dict[cat]:
                if ('f' in opts.to_calculate) or ('b' in opts.to_calculate):
                    postfix = opts.postfix + 'fungi'
                    refseq_mapping(
                        opts.mode, opts.e, cat, contig, postfix, opts.db_refseq_plant[0], tax_name_dict,
                        tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq_plant[1]
                        )
                if ('p' in opts.to_calculate) or ('b' in opts.to_calculate):
                    postfix = opts.postfix + 'plant'
                    refseq_mapping(
                        opts.mode, opts.e, cat, contig, postfix, opts.db_refseq_fungi[0], tax_name_dict,
                        tax_id_dict, opts.threads, opts.out_dir, opts.db_refseq_fungi[1]
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
                        fq2fa_com = '%s < %s > %s -Q33'%(PATH_FQ2FA, pjoin(path, fastq), pjoin(path, fasta))
                        print fq2fa_com
                        if opts.mode == 'run':
                            system(fq2fa_com)
                    fasta_dict[path].append(fasta)
        rap_paths = []
        for path in fasta_dict:
            # WARNING!!! HARDCODE
            if 'seq' in path:
                new_path = '_'.join(split(path, '_')[:-1])
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
                    usearch(opts.mode, opts.e, '16S', pjoin(cat, fasta), opts.db_16S, opts.threads)
                if 'ITS' in opts.to_calculate:
                    usearch(opts.mode, opts.e, 'ITS', pjoin(cat, fasta), opts.db_ITS, opts.threads)


def SSU_read(loc, typ=None):
    """Extracts from specially formatted FASTA file taxonomic data
    and returns them as hierarchically organised dict.

    Headers of FASTA files should follow schemas presented by following examples:

        - Default format (used in UNITE database):
        >DQ233785|uncultured ectomycorrhizal fungus|Fungi|Thelephora terrestris|Fungi; Basidiomycota; Agaricomycotina; Agaricomycetes; Incertae sedis; Thelephorales; Thelephoraceae; Thelephora; Thelephora terrestris

        - Alternative format (a.k.a '16S'):
        >AF093247.1.2007 Eukaryota;Amoebozoa;Mycetozoa;Myxogastria;;Hyperamoeba_sp._ATCC50750

    Args:
        loc: location - path to the FASTA file

        typ: tells function which format of header is present in given file.
            If specified, indicates use of alternative format. Default: None.

    Returns:
        A dict with taxonomic data, where:
            - keys are sequence identifiers;
            - values are lists of taxonomic terms, in order:
                from the most generic to the most specific one.
        Example:
        {
            'DQ482017':
            [
                'Fungi',
                'Basidiomycota',
                'Agaricomycotina',
                'Agaricomycetes',
                'Incertaesedis',
                'Thelephorales',
                'Thelephoraceae',
                'Tomentella',
                'Tomentellasublilacina'
            ],
            'EF031133':
            [
                'Fungi',
                'Basidiomycota',
                'Agaricomycotina',
                'Agaricomycetes',
                'Incertaesedis',
                'Thelephorales',
                'Thelephoraceae',
                'Thelephora',
                'Thelephoraterrestris'
            ]
        }

    """
    tax_dict = {}
    with open(loc, 'r') as fasta_file:
        for linia in fasta_file.readlines():
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
                    else:
                        tax_line = ""
                    if tax_line == '-':
                        tax_line = 'Fungi'+';'+linia.split('|')[1].replace(" ", "_")
                    tax_dict[tax_idx] = tax_line.split(';')
    return tax_dict


def tuple_to_dict(tuple_dict):
    """ Input dict has tuples as keys and int type values.
    In the dict returned by function:
       - keys are elements from input dict's tuples,
       - values are dicts.


    Args:
        tuple_dict: a dict.
        Example: {('first','second','third'): 10, ('first','second'): 30, ...}

    Example will show why this function can be useful:

    input:
        {
            ('Animalia','Mammalia','Canis lupus familiaris'): 10,
            ('Bacteria','Enterobacteriaceae','Escherichia coli'): 20,
            ('Animalia','Mammalia','Felis catus'):30
        }
    output:
        {
            'Animalia':
            {
                'Mammalia':
                {
                    'Canis lupus familiaris': {'subsum': 10},
                    'Felis catus': {'subsum': 30}
                }
            }
            'Bacteria':
            {
                'Enterobacteriaceae':
                {
                    'Escherichia coli': {'subsum': 20}
                }
            }
        }
    """
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
    """Input dict has tuples as keys and int type values.
    In dict returned by function all elements from input dict's
    tuples are single keys but their current values
    are the sum of the values of all tuples, in which they were.

    Arg: tuple_dict
    example: {
                ('first','second','third'):10,
                ('first','third','fifth':30)
                ...
            }


    Example:
        input:
            {
                ('Animalia','Mammalia','Canis lupus familiaris'):10,
                ('Bacteria','Enterobacteriaceae','Escherichia coli'):20,
                ('Animalia','Mammalia','Felis catus'):30
            }
        output:
            {
                'Animalia': 40,
                'Felis catus': 30,
                'Enterobacteriaceae': 20,
                'Mammalia': 40,
                'Escherichia coli': 20,
                'Bacteria': 20,
                'Canis lupus familiaris': 10
            }
    """
    fin_dict = {}
    for taxonomy in tuple_dict:
        for level in taxonomy:
            if level in fin_dict.keys():
                fin_dict[level] += tuple_dict[taxonomy]
            else:
                fin_dict[level] = tuple_dict[taxonomy]
    return fin_dict


def dict_purify(bac_dict):
    """Removes from given dict values below the hardcoded threshold.

    Args:
        bac_dict: a dict.

    Returns:
        A new dict without entries, within values were lower than threshold.

    HARDCODED:
        threshold = 10
    """
    threshold = 10
    all_keys = bac_dict.keys()
    to_remove = []
    for key in all_keys:
        if bac_dict[key] < threshold:
            to_remove.append(key)
    for key in to_remove:
        del bac_dict[key]
    return bac_dict


def file_analysis(typ, name, SSU=None):
    """Using given SSU as database, performs 'statistical' analysis
    of taxonomy from file <name>. It counts occurrences of different
    species in given file and returns result in form of two dicts.

    The threshold hardcoded in dict_purify is used,
    to remove neglectable entities.

    Args:
        typ: an "output type" - typically: 'ITS', '16S' or 'txt'.

        name: name of file to analise
            TODO: There was a suggestion, that there might be something
            wrong with this variable, when passed in only file_analysis
            call present in this file. The issue is linked to "full path
            or basename?" question. It should be investigated, in fully
            functional testing environment.

        SSU: A dict with taxonomic data, where:
                keys are sequence identifiers; values are lists of taxonomic
                terms, in order: from the most generic to the most specific one.
            Explicit examples included in description of SSU_read() function.

    Returns:
        A 2-tuple:

            If file of name <name> was not found:
                ('NA', 'NA')
            If data of type '16S' was successfully analysed:
                (full_bac_dict, tax_bac_dict)
            If data of type 'ITS' was successfully analysed:
                (full_fun_dict, tax_fun_dict)

        full_bac_dict / full_fun_dict: a dict of dicts.
            It represents taxonomic tree. Generally in this dict:
                - keys are names of taxa,
                - values are nested dicts of the same type.
            Values in the deepest levels ("leafs") are counts of occurrences.

            Example:
            {
                'Bacillaceae':
                {
                    'Anoxybacillus': {'subsum': 15}
                },
                'Listeriaceae':
                {
                    'Listeria':
                    {
                        'Lgrayi': {'subsum': 22},
                        'Linnocua': {'subsum': 11}
                    }
                }
            }

        tax_bac_dict / tax_fun_dict: a dict, with taxonomic statistics:
            - keys are names of taxa
            - values are counts of occurrences

            Example:
            {
                'Bacillaceae': 15,
                'Anoxybacillus': 15,
                'Listeriaceae': 33,
                'Listeria': 33,
                'Lgrayi': 22,
                'Linnocua': 11
            }


    HARDCODED:
        In '16S' analysis, if spec is 'Phaseolus_acutifolius_(tepary_bean)',
        then it is not counted neither as bacteria nor archaea.
    """
    if not pexists(name):
        return 'NA', 'NA'
    else:
        with open(name, 'r') as content:
            linie = content.readlines()
            # I have commented following block as a result of
            # discussion on slack with siwek (Michal)
            # if typ in ['fmc', 'pmc']:
            #     mapped, unmapped = linie[0].strip().split('-')
            #     try:
            #         mapped = int(mapped)
            #     except:
            #         mapped = 'NA'
            #     try:
            #         unmapped = int(unmapped)
            #     except:
            #         unampped = 'NA'
            #     if mapped == 'NA' or unmapped == 'NA':
            #         total = 'NA'
            #     else:
            #         total = mapped + unmapped
            #     return mapped, total
            # if typ == 'log':
            #     data = linie[-1].split(' ')
            #     try:
            #         n50 = int(data[8][:-1])
            #         total = int(data[12][:-1])
            #         using = split(data[14], '/')[0]
            #     except:
            #         n50 = total = using = 'NA'
            #     return n50, total, using
            if typ == '16S':
                bac_arch = [0, 0]   # counts: 0 - of bacteria 1 - of archaea
                bac_dict = {}   # counts species/strains of bacteria
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
                            spec = tuple(SSU[tax_id][:spec_idx+1])
                            cult_control = 1
                    if SSU[tax_id][0] == 'Bacteria':
                        bac_arch[0] += 1
                    elif SSU[tax_id][0] == 'Archaea':
                        bac_arch[1] += 1
                    else:
                        pass
                    if spec:
                        if spec == 'Phaseolus_acutifolius_(tepary_bean)':
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
    """Generates lists of locations where input files are located.
    Lists are grouped by directories and later, by "output types".

    As input files considered are only files,
    which meet all the following conditions:
        - they are located in current working directory
          or in subdirectories of current working directory,
        - have suffix of filename equal to <out_type>,
        - if <out_type> is 'ITS' or '16S', filenames also
          have to contain 'usearch_' before the <out_type> in name.

    Args:
        mode: a parameter passed to cat_read function.
            If mode is 'run', the gunzip command will be executed,
            to extract compressed files.

        out_types: a list with "output types",
            typically consisting  of: 'ITS', '16S' or 'txt'.

    Returns:
        A dict where keys are "output types", values are dicts describing
        file locations returned by cat_read() function.

        Example:
        {
            'ITS':
            {
                directory_path1:
                    [file1_in_this_directory_path,
                    file2_in_this_directory_path,
                    file3_in_this_directory_path],
                directory_path2:
                    [file4_in_this_directory_path,
                    file5_in_this_directory_path,
                    file6_in_this_directory_path]
            }
        }
    """
    in_dict = {}
    for out_type in out_types:
        if out_type in ['ITS', '16S']:
            in_dict[out_type] = cat_read(mode, 'usearch_' + out_type, False)
        else:
            in_dict[out_type] = cat_read(mode, out_type, False)
    return in_dict


def dict_prepare(typ, indict, SSU):
    """Runs files analysis and then, creates two dicts with summarized
    taxonomic data. The resulting dicts keep information about files
    from which given taxa comes.

    Args:
        typ: an "output type" - typically: 'ITS', '16S' or 'txt'.

        indict: A dict with file structure, where:
                keys are directories,
                values are lists of files.
            Example in description od output of cat_read() function.

        SSU: A dict with taxonomic data, where:
                keys are sequence identifiers;
                values are lists of taxonomic terms, in order:
                from the most generic to the most specific one.
            Example in description of output of SSU_read() function.

    Returns:
        A tuple: (all_dicts, all_tax_dict)

        all_dicts: A dict, where:
            keys are file identifiers,
            values are taxonomic trees with numbers of occurrences
            of particular species placed inside leafs

            Example:
            {
                'filename_1.ext':
                {
                    'Bacillaceae': {'Anoxybacillus': {'subsum': 15}}
                },
                'filename_2.ext':
                {
                    'Listeriaceae':
                    {
                        'Listeria':
                        {
                            'Lgrayi': {'subsum': 22},
                            'Linnocua': {'subsum': 11}
                        }
                    }
                }
            }

        all_tax_dict: A dict, where:
            keys are file identifiers,
            values are dicts, where:
                keys are names of taxa,
                values are numbers of occurrences of particular taxon.

            Example:
            {
                'filename_1.ext':
                {
                    'Bacillaceae': 5, 'Anoxybacillus': 5
                },
                'filename_2.ext':
                {
                    'Listeriaceae': 19, 'Listeria': 19, 'Lgrayi': 16, 'Linnocua': 3
                }
            }
    """
    all_dicts = {}
    all_tax_dict = {}
    for path in indict.keys():
        for plik in indict[path]:
            plik_id = '_'.join([typ, split(plik, '.')[0]])
            analysed_dict, tax_dict = file_analysis(typ, pjoin(path,plik), SSU)
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
    If during walking the dict, the 'subsum' string is found,
    recurrence stops and nodes beneath are ignored.
    The dicts have to contain only other dicts or 'subsum'!

    Args:
        tax_tree: a dict to update
        curr_tax: a dict to be merged into the <tax_tree>

    Returns:
        tax_tree: updated dict

    Example:
        input (tax_tree, curr_tax):
            (
                {
                    'Listeriaceae':
                    {
                        'Listeria': { 'Lgrayi': {} }
                    }
                },
                {
                    'Listeriaceae':
                    {
                        'Listeria': { 'Linnocua': {'subsum': 11} }
                    }
                }
            )
        output:
            {
                'Listeriaceae':
                {
                    'Listeria':
                    {
                        'Lgrayi': {},
                        'Linnocua': {}
                    }
                }
            }
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
    type and then by file into two dicts:
        - first represents full taxonomic tree
        - second presents how many species are inside particular files

    Args:
        full_dict: a dict of dicts, with structure like:
            {
                'output_type_1':
                {
                    'filename_1.ext':
                    {
                        'Bacillaceae': {'Anoxybacillus': {'subsum': 15}}
                    }
                }
            }
            Longer example available in description of xml_format().

    Returns:
        A tuple (full_tree, file_total_count)

        full_tree: A dict without information about type and file.
            Example:
            {
                'Bacillaceae':
                {
                    'Anoxybacillus': {}
                },
                'Listeriaceae':
                {
                    'Listeria':
                    {
                        'Lgrayi': {}, 'Linnocua': {}
                    }
                }
            }

        file_total_count: A dict within information about files are kept
            and values of every leaf are summed into values representing
            total count of species inside particular file.
            Example:
            {'filename_1.ext': 15, 'filename_2.ext': 33}
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
    """Creates a list of simplified filenames (full filenames comes
    from <full_dic>). The list is guaranteed to not have any duplicates.

    If filename contains sequence of nucleotides from set {A, C, T, G},
    separated from other parts of the name by underscore (_), then
    simplified filename will be substring of filename from beginning
    to this sequence (inclusive). This set of nucleotides in most cases
    will represent "index" from naming schema for fastq files, so
    since this step, many files will belong to one alias.

    Otherwise, the first chunk of filename will be used. In this case,
    filename will be splitted with dot '.' as delimiter. It also
    indicates grouping files with similar filenames since this step.

    Keep in mind, that this step groups files referring to the same
    dataset under single name, and technically is vulnerable for errors:
    If filenames given to this functions do not follow appropriate
    naming schema or when these filenames are not named in prefix-code
    convention, then some false-positives might be generated in
    some of functions which use results of xml_name_parse().

    Example:
        Following filenames:
            '16S_ArchV3V4_M_BF_02_TAGCTT_L001_001.amplicons.cutadapt.flash.merged.fastq.extendedFrags.fasta.usearch_16S'
            '16S_ArchV3V4_M_BF_02_TAGCTT_L001_001.cutadapt.amplicons.cutadapt.flash.merged.fastq.extendedFrags.fasta.usearch_16S'
            'Amp45_BFp_B_CAAAAG_L001_R12.fasta.usearch_ITS'
            'meta-velvetg.contigs.fa.usearch_ITS'
        will become:
            '16S_ArchV3V4_M_BF_02_TAGCTT' # note: this groups first two files
            'Amp45_BFp_B_CAAAAG'
            'meta-velvetg'

    Args:
        full_dic: a dict of dicts, with structure like:
            {
                'output_type_1':
                {
                    'filename_1.ext':
                    {
                        'Bacillaceae': {'Anoxybacillus': {'subsum': 15}}
                    }
                }
            }
            Longer example avaliable in description of xml_format().

    Returns:
        A list with simplified filenames.
        The list is guaranteed to not have any duplicates.

        Example:
        ['filename_1', 'filename_2']

    """
    # Make set for names.
    # It allows to get rid of duplicates, but order of elements is lost.
    name_set = set()
    for typ in full_dic:
        for plik in full_dic[typ]:
            # Remove extension from filename
            display_name = split(plik, '.')[0]
            try:
                nameparts = split(plik, '_')
                nucleotides = set('ACTG')
                # let identifier be the first part of filename, that is
                # comprised only of nucleotides ACTG.
                identifier = [item for item in nameparts if not set(item).difference(nucleotides)][0]
                # CHANGED
                # Due to occuring error:
                #   TypeError: cannot concatenate 'str' and 'int' objects in original line:
                # display_name = '_'.join(nameparts[0:nameparts.index(identifier+1)])
                # and with acceptance from siwiak, I made a small correction:
                display_name = '_'.join(nameparts[0:nameparts.index(identifier)+1])
            except:
                pass
            name_set.add(display_name)
    name_list = list(name_set)
    return name_list


def xml_vals(xml_names, tax_dict):
    """Reformats information extracted from tax_dict into another dict,
    to allow use of this data in xml-generation process.

    Keep in mind, that this step - along with generation of xml_names
    in xml_name_parse - groups files referring to the same dataset under
    a single name, and technically is vulnerable for some errors:
    If filenames given to these functions do not follow appropriate
    naming schema or when these filenames are not named in prefix-code
    convention, then some false-positives might be generated when
    grouping results and this will influence results.

    Args:
        xml_names: readable identifiers of groups of files referring to
            common dataset, derived from filenames.
            Example ['filename_1', 'filename_2']

        tax_dict: A dict, where:
            keys are file identifiers,
            values are dicts, where:
                keys are names of taxa,
                values are numbers of occurrences of particular taxon.

            Example:
            {
                'filename_1.ext':
                {
                    'Bacillaceae': 5, 'Anoxybacillus': 5
                },
                'filename_2.ext':
                {
                    'Bacillaceae': 1, 'Anoxybacillus': 1
                    'Listeriaceae': 19, 'Listeria': 19, 'Lgrayi': 16, 'Linnocua': 3
                }
            }

    Returns:
        A dict, where:
            - keys are names of taxa,
            - values are lists with counters of occurrences of particular
              taxon in subsequent groups of files. Order on this list
              is defined by order of names in xml_names list.

        Example:
        {
            'Anoxybacillus': [5, 1],
            'Bacillaceae': [5, 1],
            'Lgrayi': [0, 16],
            'Linnocua': [0, 3],
            'Listeria': [0, 19],
            'Listeriaceae': [0, 19]
        }
    """
    # create a set of all taxa from tax_dict
    all_tax_set = set()
    for plik in tax_dict:
        for level in tax_dict[plik]:
            all_tax_set.add(level)

    xml_dict = {}
    for tax in all_tax_set:
        xml_dict[tax] = []
        for name in xml_names:
            for plik in tax_dict:
                # Following if might generate false positives,
                # but it is partially desired. If there is a lot of
                # lanes, reads, etc for given dataset then everything
                # will be grouped under a single "name", thanks
                # to how the names are generated in xml_name_parse()
                # function, and thanks to this line:
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
    return reparsed.toprettyxml(indent="    ")


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
    """Generates xml string, dedicated to use with Krona.
    Xml is generated with use of default python xml module.

    Args:

        xml_names: a list of readable identifiers of datasets.
            Example: ['filename_1', 'filename_2']

        xml_dict: a dict with lists of occurrences of nodes by node.
            Keys are names of nodes. Values are lists of constant length,
            where every item on position *X* means: number of occurrences
            of node in dataset on position *X* on list xml_names.

            Note, that for all i: len(xml_dict[i]) is equal to len(xml_names).

            Example:
            {node_1: [count_1, count_2, ... count_x],
             node_2: [count_1, count_2, ... count_x]}

        tax_tree: a dict of dicts, etc. In form of nested dicts,
            represents phylogenetic tree.
            Example:
            {a:{aa:{}, ab:{aba:{}, abb:{abba:{}}}}}

        name_total_count: a dict with numbers of occurrences of nodes
            by dataset. Keys are names of dataset from xml_names, values are counts.
            Example:
            {identifier_1: count_1, identifier_2: count_2}

        unit: A string. Defines units to be used in the Krona.
            Default = 'reads'


    Returns:
        A single string with generated xml.

    HARDCODED:
        Number of levels to iterate in tax_tree is hardcoded to 9.
        As original comment states:
        'however it ignores all absent levels! which is actually kind of cool'
    """
    # Create a xml
    import xml.etree.ElementTree as ET
    krona = ET.Element('krona')

    # Set units
    attributes = ET.SubElement(krona, 'attributes', magnitude=unit)
    attribute = ET.SubElement(attributes, 'attribute', display=unit)
    attribute.text = unit

    # Put dataset names into the xml
    datasets = ET.SubElement(krona, 'datasets')
    for name in xml_names:
        dataset = ET.SubElement(datasets, 'dataset')
        dataset.text = name

    # Put counts adequate to corresponding names into the xml
    node = ET.SubElement(krona, 'node', name='%s count'%(unit))
    reads = ET.SubElement(node, unit)
    for name in xml_names:
        val = ET.SubElement(reads, 'val')
        val.text = str(name_total_count[name])

    # CHANGED
    # I have commented this line, after consultation on slack (Michal)
    # recurse_depth = dict_depth(tax_tree)

    # Translate tex_tree into xml
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
                                val = ET.SubElement(reads6, 'val')
                                val.text = str(element)
                            tmp_7_d = tmp_6_d[lvl_6]
                            for lvl_7 in tmp_7_d:
                                node7 = ET.SubElement(node6, 'node', name=deunique(lvl_7))
                                reads7 = ET.SubElement(node7, unit)
                                for element in xml_dict[lvl_7]:
                                    val = ET.SubElement(reads7, 'val')
                                    val.text = str(element)
                                tmp_8_d = tmp_7_d[lvl_7]
                                for lvl_8 in tmp_8_d:
                                    node8 = ET.SubElement(node7, 'node', name=deunique(lvl_8))
                                    reads8 = ET.SubElement(node8, unit)
                                    for element in xml_dict[lvl_8]:
                                        val = ET.SubElement(reads8, 'val')
                                        val.text = str(element)
                                    tmp_9_d = tmp_8_d[lvl_8]
                                    for lvl_9 in tmp_9_d:
                                        node9 = ET.SubElement(node8, 'node', name=deunique(lvl_9))
                                        reads9 = ET.SubElement(node9, unit)
                                        for element in xml_dict[lvl_9]:
                                            val = ET.SubElement(reads9, 'val')
                                            val.text = str(element)

    # Push created xml structure to beautician and return resulting string
    return prettify(krona)


def name_total_reduction(xml_names, file_total_count):
    """Rewrites given <file_total_count> dict, to group total counts
    by names (datasets), instant of grouping by filenames.

    Keep in mind, that this step - along with generation of xml_names
    in xml_name_parse - groups files referring to the same dataset under
    single name, and technically is vulnerable for some errors.
    For more information, check description of xml_vals() function.

    Args:
        xml_names: readable identifiers of groups of files referring to
            common dataset; derived from filenames.
            Example: ['filename_1', 'filename_2']

        file_total_count: A dict within information about files are kept
            and values of every leaf are summed into values (usually
            representing total count of species inside particular file)
            Example:
                {
                    'filename_1.ext': 15,
                    'filename_2.ext': 33,
                    'filename_2_part_2.ext': 33
                }

    Returns:
        A dict, where keys are names as present in xml_names and values
        are summed total counts (usually of occurrences of different species).
        Example: {'filename_1': 15, 'filename_2': 66}
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
    """Creates a set of dicts, which are useful to generate xml files from.

    Args:
        full_dict: a dict of dicts with structure like:
            {
                'output_type_1':
                {
                    'filename_1.ext':
                    {
                        'Bacillaceae': {'Anoxybacillus': {'subsum': 15}}
                    },
                    'filename_2.ext':
                    {
                        'Listeriaceae':
                        {
                            'Listeria':
                            {
                                'Lgrayi': {'subsum': 22},
                                'Linnocua': {'subsum': 11}
                            }
                        }
                    }
                },
                'output_type_2':
                {
                    'filename_2.ext':
                    {
                        'Bacillaceae': {'Anoxybacillus': {'subsum': 15}}
                    }
                }
            }
            "output type" typically will be 'ITS', '16S'.

        tax_dict: A dict, where:
            keys are file identifiers,
            values are dicts, where:
                keys are names of taxa,
                values are numbers of occurrences of particular taxon.

            Example:
            {
                'filename_1.ext':
                {
                    'Bacillaceae': 5, 'Anoxybacillus': 5
                },
                'filename_2.ext':
                {
                    'Bacillaceae': 1, 'Anoxybacillus': 1
                    'Listeriaceae': 19, 'Listeria': 19, 'Lgrayi': 16, 'Linnocua': 3
                }
            }


    Returns:
        A tuple: (xml_names, xml_dict, tax_tree, name_total_count)

        xml_names: readable identifiers of groups of files referring to
            common dataset; derived from filenames.
            Example: ['filename_1', 'filename_2']

        xml_dict: a dict, where:
            - keys are names of taxa,
            - values are lists with counters of occurences of particular
              taxon in subsequent groups of files. Order on this list
              is defined by order of names in xml_names list.

            Example:
            {
                'Anoxybacillus': [5, 1],
                'Bacillaceae': [5, 1],
                'Lgrayi': [0, 16],
                'Linnocua': [0, 3],
                'Listeria': [0, 19],
                'Listeriaceae': [0, 19]
            }

        tax_tree: a dict without information about type and file.
            Example:
            {
                'Bacillaceae':
                {
                    'Anoxybacillus': {}
                },
                'Listeriaceae':
                {
                    'Listeria':
                    {
                        'Lgrayi': {}, 'Linnocua': {}
                    }
                }
            }


        name_total_count: a dict, where keys are names as present in
            xml_names and values are summed total counts (usually of
            occurrences of different species).
            Example: {'filename_1': 5, 'filename_2': 20}
    """
    # Rewrites the data from a dict containing information arranged by
    # type and then by file into two dicts:
    #   - first represents full taxonomic tree
    #   - second presents how many species are inside particular files
    tax_tree, file_total_count = tree_of_life(full_dict)

    # Creates a list of simplified filenames (full filenames comes
    # from <full_dic>). The list is guaranteed to not have duplicates.
    xml_names = xml_name_parse(full_dict)

    # Reformats information extracted from tax_dict into another dict,
    # to allow use of this data in xml-generation process.
    xml_dict = xml_vals(xml_names, tax_dict)

    # Rewrites given <file_total_count> dict, to group total counts
    # by names (datasets), instant of grouping by filenames.
    name_total_count = name_total_reduction(xml_names, file_total_count)

    return xml_names, xml_dict, tax_tree, name_total_count


def out_namespace(curr, out_type):
    """"Generates paths to places, where krona xml and krona html files
    are or will be placed. Path composition varies, according to given
    arguments and always contains some string representation of <out_type>.
    Paths may start in current working directory or in directory specified by <curr>.

    A basename is generated accordingly to the following table:

        out_type                |   basename (filename without extension)

        'txt'                   |  'humann-graphlan'
        a str other than 'txt'  |  out_type
        a list of strings       |  '_'.join(out_type)

    Examples of generated basenames, basing on <out_type>:

        for 'ITS':                ('ITS.krona', 'ITS.html')
        for ['ITS', '16S']:       ('ITS_16S.krona', 'ITS_16S.html')
        for 'txt':                'humann-graphlan'


    Args:
        curr: A string - path to directory where the files are/will be placed.
              If curr is 'in_situ', the path will start in current working directory.

        out_type: A variable that determines the basename of file.

            In practice, we expect here: 'ITS', '16S', 'txt'
            or some combination of first two in form of a list.

    Returns:
        A tuple (out_xml, out_html):

        out_xml: a string with path to file. The filename is derived
            from (out_type) and file has '.krona' extension.

        out_html: a string with path to file. The filename is derived
            from (out_type) and file has '.html' extension.
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


def prepare_taxonomy_stats(opts):
    """Performs statistical analysis of taxonomy from appropriate files
    from current working directory: counts occurrences of different taxa
    and prepares the results to be presented in HTML format.
    Results will be converted to HTML (with the Krona program),
    but only when <opts.mode> is set to 'run'.

    Args:
        opts: A namespace, where:

            opts.output_type:
                Allows to choice on which files the analysis will be
                performed and also determines basenames of output files.

                One of: ['ITS', '16S']

            opts.mode:
                A mode in which the program will be run.

                One of: ['test', 'run']

            opts.out_dir:
                Indicates, where the output files should be located.
                To specify current working directory, use 'in_situ'.

                One of: ['in_situ', a_string_with_path_to_dir]

            opts.db_taxonomy_16S:
                Path to 16S database used in taxonomy classification
                (fasta with specially formatted headers)

            opts.db_taxonomy_ITS:
                Path to ITS database used in taxonomy classification
                (fasta with specially formatted headers)

    Input:

        Analysis will be performed on files, meeting all the following conditions:
            - files are located in current working directory
              or in subdirectories of current working directory,
            - suffix of filename is equal to <opts.output_type>,
            - if <opts.output_type> is 'ITS' or '16S', filenames
             contains 'usearch_' before the <opts.output_type> in name.

    Output:

        Output files will be placed in <opts.out_dir> directory, with
        basenames depending on <opts.output_type>.

        Examples of filenames:
            for 'ITS':                ('ITS.krona', 'ITS.html')
            for ['ITS', '16S']:       ('ITS_16S.krona', 'ITS_16S.html')

    """
    SSU = {}
    # Extracts from specially formatted FASTA file taxonomical data
    # and returns them as hierarchically organised dict
    if '16S' in opts.output_type:
        SSU['16S'] = SSU_read(opts.db_taxonomy_16S, '16S')
    if 'ITS' in opts.output_type:
        SSU['ITS'] = SSU_read(opts.db_taxonomy_ITS)
    # Generates list of locations were input files are located.
    input_dic = input_locations(opts.mode, opts.output_type)
    analysed_dict = {}
    tax_dict = {}
    pure_tax = {}
    for out_type in opts.output_type:
        # Runs files analysis on data which were interpreted by SSU_read
        # and then, creates dicts with summarized taxonomic data
        analysed_dicto, tax_dicto = dict_prepare(out_type, input_dic[out_type], SSU[out_type])
        analysed_dict[out_type] = analysed_dicto
        tax_dict[out_type] = tax_dicto
        pure_tax.update(tax_dicto)
    krona_xml_name, krona_html_name = out_namespace(opts.out_dir, opts.output_type)
    xml_names, xml_dict, tax_tree, name_total_count = xml_format(analysed_dict, pure_tax)
    krona_unit = 'reads'
    chdir('../')
    xml_string = xml_prepare(xml_names, xml_dict, tax_tree, name_total_count, krona_unit)
    # Writes xml_string into the file given by out_namespace
    outprint(xml_string, krona_xml_name)
    if opts.mode == 'run':
        krona_to_html_comm = 'ktImportXML -o %s %s'%(krona_html_name, krona_xml_name)
        system(krona_to_html_comm)
