#!/usr/bin/python
import fileinput
from os import path, getcwd

common_header = """\
set useParseTextTaxonomy=false;
load taxGIFile='/home/pszczesny/workingdata/mapping/gi_taxid_prot.dmp';
set taxUseGIMap=true;
"""

template = """\
import blastFile='{0}'
fastaFile='{1}'
meganFile='{2}'
"""

common_options = "\
maxMatches=100\
 minScore=50.0\
 maxExpected=0.01\
 topPercent=10.0\
 minSupport=50\
 minComplexity=0.44\
 useMinimalCoverageHeuristic=false\
 useSeed=true\
 useCOG=true\
 useKegg=true\
 paired=false\
 useIdentityFilter=false\
 blastFormat=RapSearch\
 mapping='Taxonomy:GI_MAP=true'"

options_heavy = common_options + 'textStoragePolicy=inRMA;\n'
options_light = common_options + 'textStoragePolicy=inOriginal;\n'

del common_options

cwd = getcwd()

path_heavy = path.join(cwd, 'commands.rma')
path_light = path.join(cwd, 'commands.little.rma')

with open(path_heavy, 'w') as file_heavy, open(path_light, 'w') as file_light:

    both = (file_heavy, file_light)
    write_to_both = lambda contents: map(lambda f: f.write(contents), both)

    write_to_both(common_header)

    for line in fileinput.input():

        full_path = path.join(cwd, line.rstrip())
        directory = path.dirname(full_path)
        basename = path.basename(full_path)

        fasta_path = path.join(directory, 'meta-velvetg.contigs.fa')
        megan_path = path.join(cwd, basename + '.rma')

        write_to_both(template.format(full_path, fasta_path, megan_path))

        file_heavy.write(options_heavy)
        file_light.write(options_light)
