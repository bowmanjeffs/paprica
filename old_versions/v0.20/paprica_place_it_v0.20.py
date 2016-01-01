# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 08:59:11 2015

@author: jeff

REQUIRES:
    Files:
        bacterial_ssu.cm
        
    Programs:
        infernal
        taxtastic
        seqmagick,
        pplacer,
        raxmlHPC-PTHREADS-AVX2
        
    Python modules:
        Bio

CALL AS:
    python genome_finder_place_it.py [query] [ref] for analysis or
    python genome_finder_place_it.py [ref] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (.fasta).

"""
executable = '/bin/bash' # shell for executing commands

import re
import subprocess
import sys
import os

## Read in profile.  Required variables are ref_dir and cutoff. ###

variables = {}

def get_variable(line, variables):
    
    line = line.rstrip()
    line = line.split('=')
    variable = line[0]
    value = line[1]
    
    variables[variable] = value
    return variables

with open('paprica_profile.txt', 'r') as profile:
    for line in profile:
        if line.startswith('#') == False:
            if line != '\n':
                get_variable(line, variables)

## Define function to clean record names.

from Bio import SeqIO

def clean_name(file_name):
    bad_character = re.compile('[=-@!%,;\(\):\'\"\s]') # probably need to add \.
    with open(file_name + '.clean.fasta', 'w') as fasta_out:
        for record in SeqIO.parse(file_name + '.fasta', 'fasta'): 
            record.name = re.sub(bad_character, '_', str(record.description))
            print >> fasta_out, '>' + record.name
            print >> fasta_out, record.seq
            
## Execute main program.

if len(sys.argv) < 3:
    
    ## Add a dummy name for diagnostic testing if necessary.  This assumes
    ## that no command line arguments are provided, as if run with execfile().
    
    if len(sys.argv) == 1:   
        ref = 'combined_16S.tax'
    
    elif len(sys.argv) == 2:
        ref = sys.argv[1]
        
    clean_name(variables['ref_dir'] + ref)
    
    ## Conduct alignment with Infernal (cmalign) against the bacteria profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    
    ## Degap first, just in case.
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + variables['ref_dir'] + ref + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()

    infernal_commands = 'cmalign -o ' + variables['ref_dir'] + ref + '.clean.align.sto --outformat Pfam ' + variables['cm'] + ' ' + variables['ref_dir'] + ref + '.clean.fasta'      
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + variables['ref_dir'] + ref + '.clean.align.sto ' + variables['ref_dir'] + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()

    transcribe = subprocess.Popen('seqmagick mogrify --transcribe rna2dna ' + variables['ref_dir'] + ref + '.clean.align.fasta', shell = True, executable = executable)
    transcribe.communicate()       
    
    rm = subprocess.call('rm ' + variables['ref_dir'] + '*ref.tre', shell = True, executable = executable)
    raxml1 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + variables['cpus'] + ' -m GTRGAMMA -s ' + variables['ref_dir'] + ref + '.clean.align.fasta -n ref.tre -f d -p 12345 -w ' + variables['ref_dir'], shell = True, executable = executable)
    raxml1.communicate()
    
    ## Root the tree.
    
    raxml2 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T 1 -m GTRGAMMA -f I -t ' + variables['ref_dir'] + 'RAxML_bestTree.ref.tre -n root.ref.tre', shell = True, executable = executable)  
    raxml2.communicate()
    
    ## Generate SH-like support values for the tree.
    
    raxml3 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + variables['cpus'] + ' -m GTRGAMMA -f J -p 12345 -t ' + variables['ref_dir'] + 'RAxML_bestTree.ref.tre -n conf.root.ref.tre -s ' + variables['ref_dir'] + ref + '.clean.align.fasta -w ' + variables['ref_dir'], shell = True, executable = executable)   
    raxml3.communicate()
     
    ## Generate the reference package using the tree with SH support values and a log file.
    
    rm = subprocess.call('rm -r ' + variables['ref_dir'] + ref + '.refpkg', shell = True, executable = executable)
    taxit = subprocess.Popen('taxit create -l 16S_rRNA -P ' + variables['ref_dir'] + ref + '.refpkg --aln-fasta ' + variables['ref_dir'] + ref + '.clean.align.fasta --tree-stats ' + variables['ref_dir'] + 'RAxML_info.ref.tre --tree-file ' + variables['ref_dir'] + 'RAxML_fastTreeSH_Support.conf.root.ref.tre', shell = True, executable = executable)
    taxit.communicate()
    
elif len(sys.argv) == 3:
    
    query = sys.argv[1]
    ref = sys.argv[2]
    
    clear_wspace = subprocess.call('rm ' + query + '.' + ref + '*', shell = True, executable = executable)
            
    combine = subprocess.Popen('cat ' + variables['ref_dir'] + ref + '.refpkg/' + ref + '.clean.align.fasta ' + query + '.fasta > ' + query + '.' + ref + '.fasta', shell = True, executable = executable)
    combine.communicate()
    
    clean_name(query + '.' + ref)
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + query + '.' + ref + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    infernal_commands = 'cmalign -o ' + query + '.' + ref + '.clean.align.sto --outformat Pfam ' + variables['cm'] + ' ' + query + '.' + ref + '.clean.fasta'      
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()      
  
    convert = subprocess.Popen('seqmagick convert ' + query + '.' + ref + '.clean.align.sto ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()
    
    transcribe = subprocess.Popen('seqmagick mogrify --transcribe rna2dna ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    transcribe.communicate() 
    
    pplacer = subprocess.Popen('pplacer --out-dir ' + os.getcwd() + ' -p --keep-at-most 20 -c ' + variables['ref_dir'] + ref + '.refpkg ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    pplacer.communicate()
    
    guppy1 = subprocess.Popen('guppy to_csv --point-mass --pp -o ' + query + '.' + ref + '.clean.align.csv ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy1.communicate()
    
    guppy2 = subprocess.Popen('guppy fat --node-numbers --point-mass --pp -o ' + query + '.' + ref + '.clean.align.phyloxml ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy2.communicate()    
    
else:
    print 'wrong number of positional arguments!'
    
    
    
