# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 08:59:11 2015

@author: jeff

Dependencies, these should all be in your path and callable as listed:
    mothur, FastTreeMP, taxtastic, seqmagick, pplacer

call as python genome_finder_place_it.py [query] [ref] for analysis or
python genome_finder_place_it.py [ref] to generate a reference package.  
[ref] or [query] includes the entire file name except .fasta.

"""
##### set user variables #####

cpus = str(24)                                                                   # number of cpus for mothur to use
ref_dir = '/volumes/hd1/ref_genome_database_v2/'                                   # location of/for reference package
align_ref = '/volumes/deming/databases/silva.combined/silva.combined.fasta'     # pathway and name of reference alignment
executable = '/bin/bash'                                                        # shell for executing commands, change for windows

##### end set user variables #####

import re
import subprocess
import sys
from Bio import SeqIO

def clean_name(file_name):
    bad_character = re.compile('[-@!\.%,;\(\):\'\"\s]')
    with open(file_name + '.clean.fasta', 'w') as fasta_out:
        for record in SeqIO.parse(file_name + '.fasta', 'fasta'):       
            record.name = re.sub(bad_character, '_', str(record.name))
            print >> fasta_out, '>' + record.name
            print >> fasta_out, record.seq
    
if len(sys.argv) == 2:    
    
    ref = sys.argv[1]
    clean_name(ref_dir + ref)

    mothur_commands = 'mothur "#align.seqs(candidate=' + ref_dir + ref + '.clean.fasta, flip=t, processors=' + cpus + ', template=' + align_ref + ');' \
    'screen.seqs(minlength=1200);' \
    'filter.seqs(trump=.,vertical=T)"'
      
    mothur = subprocess.Popen(mothur_commands, shell = True, executable = executable)
    mothur.communicate()
    
    deunique = subprocess.Popen('seqmagick mogrify --deduplicate-sequences ' + ref_dir + ref + '.clean.good.filter.fasta', shell = True, executable = executable)
    deunique.communicate()    
    
    fasttree = subprocess.Popen('FastTreeMP -nt -gtr -log ' + ref_dir + ref + '.log ' + ref_dir + ref + '.clean.good.filter.fasta > ' + ref_dir + ref + '.clean.good.filter.tre', shell = True, executable = executable)
    fasttree.communicate()
    
    taxit = subprocess.Popen('taxit create -l 16S_rRNA -P ' + ref_dir + ref + '.refpkg --aln-fasta ' + ref_dir + ref + '.clean.good.filter.fasta --tree-stats ' + ref_dir + ref + '.log --tree-file ' + ref_dir + ref + '.clean.good.filter.tre', shell = True, executable = executable)
    taxit.communicate()
    
elif len(sys.argv) == 3:
    
    ref = sys.argv[2]
    query = sys.argv[1]
            
    combine = subprocess.Popen('cat ' + ref_dir + ref + '.clean.fasta ' + query + '.fasta > ' + query + '_' + ref + '.fasta', shell = True, executable = executable)
    combine.communicate()
    
    clean_name(query + '_' + ref)
    
    mothur_commands = 'mothur "#align.seqs(candidate=' + query + '_' + ref + '.clean.fasta, flip=t, processors=' + cpus + ', template=' + align_ref + ');' \
    'filter.seqs(hard=' + ref_dir + ref + '.filter);' \
    'screen.seqs(minlength=50)"'
    mothur = subprocess.Popen(mothur_commands, shell = True, executable = '/bin/bash')
    mothur.communicate()
    
    final_clean = subprocess.Popen('tr "." "-" < ' + query + '_' + ref + '.clean.filter.good.fasta > ' + query + '_' + ref + '.pplacer.filter.fasta', shell = True, executable = executable)
    final_clean.communicate()
    
    pplacer = subprocess.Popen('pplacer -p --keep-at-most 1 -c ' + ref_dir + ref + '.refpkg ' + query + '_' + ref + '.pplacer.filter.fasta', shell = True, executable = executable)
    pplacer.communicate()
    
else:
    print 'wrong number of positional arguments!'
    
    
