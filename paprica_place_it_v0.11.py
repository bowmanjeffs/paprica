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

cpus = str(1)                                                                   # number of cpus for mothur to use
ref_dir = '/home/user/genome_finder/ref_genome_database_a/'                                   # location of/for reference package
align_ref = '/home/user/genome_finder/silva.seed_v119.align'     # pathway and name of reference alignment
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
            
def check_overlap(query, ref, ref_dir):
    with open('bad_queries.txt', 'w') as exclude:
        ref_names = set()
        position_dict = {}
        
        for record in SeqIO.parse(ref_dir + ref + '.clean.fasta', 'fasta'):
            ref_names.add(record.name)
        
        ## determine the number of reference bases in each position of the alignment
        
        for record in SeqIO.parse(query + '.' + ref + '.pplacer.filter.fasta', 'fasta'):
            if record.name in ref_names:
                for i,j in enumerate(record.seq):
                    if j == '-':
                        continue
                    else:
                        try:
                            position_dict[i] = position_dict[i] + 1
                        except KeyError:
                            position_dict[i] = 1
                            
        ## determine which query seqs don't share any positions with reference alignment
                            
        for record in SeqIO.parse(query + '.' + ref + '.pplacer.filter.fasta', 'fasta'):
            if record.name not in ref_names:            
                shared = []
                
                for i,j in enumerate(record.seq):
                    if j == '-':
                        continue
                    else:
                        try:
                            temp = position_dict[i]
                        except KeyError:
                            temp = 0
                        shared.append(temp)
                        
                if sum(shared) == 0:
                    print >> exclude, record.name
                    
        ## eliminate any queries that share 0 positions with alignment
                    
    make_exclusion = subprocess.Popen('seqmagick mogrify --exclude-from-file bad_queries.txt ' + query + '.' + ref + '.pplacer.filter.fasta', shell = True, executable = executable)
    make_exclusion.communicate()
    
if len(sys.argv) == 2:    
    
    ref = sys.argv[1]
    clean_name(ref_dir + ref)
    
    ## first build a tree with appropriate filtering to generate a good tree
    ## this is the tree that taxit will use to build the reference pacakage

    mothur_commands = 'mothur "#align.seqs(candidate=' + ref_dir + ref + '.clean.fasta, flip=t, processors=' + cpus + ', template=' + align_ref + ');' \
    'screen.seqs(minlength=1200);' \
    'filter.seqs(trump=.,vertical=T,soft=50)"'
      
    mothur = subprocess.Popen(mothur_commands, shell = True, executable = executable)
    mothur.communicate()
    
    deunique = subprocess.Popen('seqmagick mogrify --deduplicate-sequences ' + ref_dir + ref + '.clean.good.filter.fasta', shell = True, executable = executable)
    deunique.communicate()    
    
    fasttree = subprocess.Popen('FastTreeMP -nt -gtr ' + ref_dir + ref + '.clean.good.filter.fasta > ' + ref_dir + ref + '.clean.good.filter.tre', shell = True, executable = executable)
    fasttree.communicate()
    
    ## then build a tree without filtering to obtain rate categories
    ## the log from this tree will be used by taxit to build the reference package
    ## this is also the alignment that will be used to build the reference package
    ## the actual tree is labeled "bad" and should not be used
    
    ## first we need to get a list of seq ids from the alignment used to build the tree
    ## so that we can make sure other seqs are not included, as the different filtering
    ## will have made them differentially unique
    
    get_ids = subprocess.Popen('seqmagick extract-ids ' + ref_dir + ref + '.clean.good.filter.fasta > deunique.ids', shell = True, executable = executable)
    get_ids.communicate()
    
    mothur_commands = 'mothur "#align.seqs(candidate=' + ref_dir + ref + '.clean.fasta, flip=T, processors=' + cpus + ', template=' + align_ref + ');' \
    'screen.seqs(minlength=1200);' \
    'filter.seqs(trump=.,vertical=T)"'
      
    mothur = subprocess.Popen(mothur_commands, shell = True, executable = executable)
    mothur.communicate()
    
    deunique = subprocess.Popen('seqmagick mogrify --include-from-file deunique.ids ' + ref_dir + ref + '.clean.good.filter.fasta', shell = True, executable = executable)
    deunique.communicate()    
    
    fasttree = subprocess.Popen('FastTreeMP -nt -gtr -log ' + ref_dir + ref + '.log ' + ref_dir + ref + '.clean.good.filter.fasta > ' + ref_dir + ref + '.bad.tre', shell = True, executable = executable)
    fasttree.communicate()
    
    ## generate the reference package using the good tree and good log file
    
    taxit = subprocess.Popen('taxit create -l 16S_rRNA -P ' + ref_dir + ref + '.refpkg --aln-fasta ' + ref_dir + ref + '.clean.good.filter.fasta --tree-stats ' + ref_dir + ref + '.log --tree-file ' + ref_dir + ref + '.clean.good.filter.tre', shell = True, executable = executable)
    taxit.communicate()
    
elif len(sys.argv) == 3:
    
    ref = sys.argv[2]
    query = sys.argv[1]
            
    combine = subprocess.Popen('cat ' + ref_dir + ref + '.clean.fasta ' + query + '.fasta > ' + query + '.' + ref + '.fasta', shell = True, executable = executable)
    combine.communicate()
    
    clean_name(query + '.' + ref)
    
    mothur_commands = 'mothur "#align.seqs(candidate=' + query + '.' + ref + '.clean.fasta, flip=T, processors=' + cpus + ', template=' + align_ref + ');' \
    'filter.seqs(vertical=T);' \
    'screen.seqs(minlength=50)"'
    mothur = subprocess.Popen(mothur_commands, shell = True, executable = '/bin/bash')
    mothur.communicate()
    
    final_clean = subprocess.Popen('tr "." "-" < ' + query + '.' + ref + '.clean.filter.good.fasta > ' + query + '.' + ref + '.pplacer.filter.fasta', shell = True, executable = executable)
    final_clean.communicate()
    
    check_overlap(query, ref, ref_dir)
    
    pplacer = subprocess.Popen('pplacer -p --keep-at-most 20 -c ' + ref_dir + ref + '.refpkg ' + query + '.' + ref + '.pplacer.filter.fasta', shell = True, executable = executable)
    pplacer.communicate()
    
else:
    print 'wrong number of positional arguments!'
    
    
