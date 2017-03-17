#!/usr/bin/python
import os
import argparse
import shutil
from metatranscriptomics_bipype import metatranscriptomics
from refseq_bipype import sample, prepare_taxonomy_stats
from settings_bipype import *


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        fromfile_prefix_chars='@',
        description='bipype stands for BioInformatics-PYthon-PipE',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='All commands may be presented in a configuration file, '
               'fed to bipype with @ prefix as "bipype @my_conf_file", '
               'my_conf_file should contain all desired commands and their'
               'options (if applicable) one per line')

    general = parser.add_argument_group(
        'general',
        'performance and I/O related options')
    general.add_argument(
        "-t", "--threads",
        help='number of threads to be used',
        type=int,
        default=8)
    general.add_argument(
        "-m", "--mode",
        help='available modes: test, run',
        choices=['test', 'run'],
        type=str,
        default='test')
    general.add_argument(
        "--out_dir", "-o",
        type=str,
        metavar='OUTPUT_DIRECTORY',
        default='in_situ',
        help='Directory for output files, '
             'default value is not usable for metatranscriptomics')
    general.add_argument(
        "--input", "-i",
        nargs='*',
        default=None)
    general.add_argument(
        "--ins_len",
        default='9999',
        help='insert length - be advised - you better use it for single run',
        type=int)
    general.add_argument(
        "-postfix",
        type=str,
        help='alphanumerical postfix of processed file',
        default='')
    general.add_argument(
        "-e",
        help="use existing files",
        action='store_true')

    taxonomy_stats = parser.add_argument_group(
        'taxonomy_stats options',
        'options related to process of preparing taxonomy results')
    taxonomy_stats.add_argument(
        "-ot", "--output_type",
        nargs='*',
        default=['ITS', '16S', 'txt'],
        help='Choice of files searched for an analysis, coded: '
             '16S and ITS on usearches - file .usearch_ITS, .usearch_16S')

    dataclean = parser.add_argument_group(
        'input cleaning',
        'methods of cleaning input from noise')
    dataclean.add_argument(
        "-ic", "--initial_cleaning",
        choices=['usearch', 'fastx'],
        type=str,
        default='',
        help='Choice of initial cleaning method')
    dataclean.add_argument(
        "--cutadapt",
        type=str,
        default='',
        nargs=2,
        help='-cutadapt ADAPTER_FILE search_options, '
             'location of file with adapters to be used by cutadapt '
             '(possible "use_filenames" to determine adapters from hardcode), '
             'and list of usearches to be run on created files - '
             'possible options are 16S, ITS, both. '
             'Please note, that mapping options -16S, -ITS are completely '
             'irrelevant if you use cutadapt. '
             'Other note - this is !!!IMPORTANT!!! to present location of file '
             'with adapters as first option of this argument')

    mappings = parser.add_argument_group(
        'mapping',
        'mappings to be done during the run')
    mappings.add_argument(
        "-ITS",
        help='usearch ITS database',
        dest='to_calculate',
        action='append_const',
        const='ITS')
    mappings.add_argument(
        "-16S",
        help='usearch 16S database',
        dest='to_calculate',
        action='append_const',
        const='16S')
    mappings.add_argument(
        "-refseq",
        help='map samples on refseq: p - plant, f - fungi, b - both',
        nargs='?',
        choices=['f', 'p', 'b'],
        dest='to_calculate',
        action='append',
        const='f',
        default='f')
    mappings.add_argument(
        "-MEGAN",
        help='generate .rma files using MEGAN',
        dest='to_calculate',
        action='append_const',
        const='MEGAN')

    contigs = parser.add_argument_group(
        'contig reconstruction',
        'methods of contig reconstruction to be used during the run')
    contigs.add_argument(
        "-assembler",
        default=None,
        type=str,
        help='choose assembler: MH for Megahit, MV for MetaVelvet',
        choices=['MH', 'MV', None])
    contigs.add_argument(
        "-MV",
        default=None,
        help="MetaVelvet's options:"
             "[initial k-mer size - default = 31,"
             "[final k-mer size, [step - default=2 (odd numbers)]]]",
        nargs='*')
    contigs.add_argument(
        "-reconstruct",
        default=None,
        const="rec_",
        help='reconstruct relating to database, options: database_loc, prefix',
        nargs='?')
    contigs.add_argument(
        "-humann",
        help='mapping rapsearch using humann',
        dest='to_calculate',
        action='append_const',
        const='humann')
    contigs.add_argument(
        "-rapsearch",
        help='RAPsearch on protein database',
        nargs='?',
        choices=['rap_KEGG', 'rap_prot', 'rap_KO'],
        dest='to_calculate',
        action='append',
        const='rap_prot',
        default='rap_prot')

    db_loc = parser.add_argument_group(
        'databases',
        'locations of databases to be used for searches')
    db_loc.add_argument(
        "--db_16S",
        type=str,
        default=PATH_X16S_DB,
        help='location of 16S database used in usearch (bowtie indexed')
    db_loc.add_argument(
        "--db_ITS",
        type=str,
        default=PATH_ITS_DB,
        help='location of ITS database used in usearch (bowtie indexed)')
    db_loc.add_argument(
        "--db_refseq_fungi",
        type=str,
        default=[PREF_PATH_REF_FUNGI],
        nargs='*',
        help='location of refseq database to use for fungi analysis. '
             'Up to two paths allowed. If there are multiple databases with'
             'filenames like <path_which_you_specified><second_part_of_name>, '
             'all will be loaded.')
    db_loc.add_argument(
        "--db_refseq_plant",
        type=str,
        default=[PREF_PATH_REF_PLANT_1, PREF_PATH_REF_PLANT_2],
        nargs='*',
        help='location of refseq database to use for plants analysis. '
             'Up to two paths allowed. If there are multiple databases with'
             'filenames like <path_which_you_specified><second_part_of_name>, '
             'all will be loaded.')
    db_loc.add_argument(
        "--db_NCBI_taxonomy",
        type=str,
        default=PATH_NCBI_TAXA_DB,
        help='location of cPickled NCBI_tax_id, '
             'NCBI tax names and NCBI tax ids')
    db_loc.add_argument(
        "--db_reconstruct",
        type=str,
        default=PATH_RECONSTRUCT_DB,
        help='location of database for reconstruction')
    db_loc.add_argument(
        "--db_taxonomy_16S",
        type=str,
        default=PATH_16S_DATABASE,
        help='location of 16S database used in taxonomy classification '
             '(FASTA with specially formatted headers)')
    db_loc.add_argument(
        "--db_taxonomy_ITS",
        type=str,
        default=PATH_ITS_DATABASE,
        help='location of ITS database used in taxonomy classification '
             '(FASTA with specially formatted headers)')

    metatr = parser.add_argument_group(
        'metatranscriptomics',
        'parameters for metatranscriptomic pipe')
    metatr.add_argument(
        "--metatr_config",
        type=str,
        help='configuration file which is necessary to launch '
             'metatranscriptomic part of bipype'),
    metatr.add_argument(
        "-mot", "--metatr_output_type",
        type=str,
        choices=['old', 'new', 'both'],
        default='both',
        help='determines type of the one of the outputs from '
             'metatranscriptomic part of bipype')
    return parser.parse_args(args)


if __name__ == '__main__':
    opts = parse_arguments()
    if opts.metatr_config:
        metatranscriptomics(opts)
    else:
        if opts.input:
            temp_path = os.getcwd()
            if not os.path.exists(temp_path):
                os.makedirs(temp_path)
            for input_file in opts.input:
                shutil.copy(input_file, os.path.join(temp_path))
        sample(opts)
        prepare_taxonomy_stats(opts)
