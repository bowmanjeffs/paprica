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
    python genome_finder_place_it.py -query [query] -ref [ref] -splits [splits] -n [nseqs] for analysis or
    python genome_finder_place_it.py -ref [ref] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (which must be .fasta).

"""
executable = '/bin/bash' # shell for executing commands

import re
import subprocess
import sys
import os
from joblib import Parallel, delayed
import datetime
import random

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
                
## Parse command line arguments.  Arguments that are unique to a run,
## such as the number of splits, should be specified in the command line and
## not in paprica_profile.txt
                
command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]
    
## Define function to clean record names.

from Bio import SeqIO

def clean_name(file_name):
    
    bad_character = re.compile('[=-@!%,;\(\):\'\"\s]') # probably need to add \.
    with open(file_name + '.clean.fasta', 'w') as fasta_out:
        for record in SeqIO.parse(file_name + '.fasta', 'fasta'): 
            record.name = re.sub(bad_character, '_', str(record.description))
            print >> fasta_out, '>' + record.name
            print >> fasta_out, record.seq
            
## Define function to split fasta file to run pplacer in parallel.  This greatly
## improves the speed of paprica_run.sh, but at the cost of memory overhead.
## Users need to be cautious of memory limits when considering the number of
## splits to make.
            
def split_fasta(file_in, nsplits):
    
    splits = []            
    tseqs = len(re.findall('>', open(file_in + '.fasta', 'r').read()))
    
    nseqs = tseqs / nsplits
    
    seq_i = 0
    file_n = 1
    
    file_out = open(file_in + '.temp' + str(file_n) + '.fasta', 'w')
    
    for record in SeqIO.parse(file_in + '.fasta', 'fasta'):
        seq_i = seq_i + 1
        
        if seq_i <= nseqs:
            SeqIO.write(record, file_out, 'fasta')
        elif seq_i > nseqs:
            splits.append(file_in + '.temp' + str(file_n))
            file_out.close()
            file_n = file_n + 1
            file_out = open(file_in + '.temp' + str(file_n) + '.fasta', 'w')
            SeqIO.write(record, file_out, 'fasta')
            seq_i = 0
        
    splits.append(file_in + '.temp' + str(file_n))
    file_out.close()
    
    return(splits)

## Define function to execute phylogenetic placement on a query fasta, or split query fasta
    
def place(query, ref, variables):
            
    clean_name(query)
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + query + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    ## Conduct alignment with Infernal (cmalign) against the bacteria profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    
    align = subprocess.Popen('cmalign --dna -o ' + query + '.clean.align.sto --outformat Pfam ' + variables['cm'] + ' ' + query + '.clean.fasta', shell = True, executable = executable)
    align.communicate()    
    
    combine = subprocess.Popen('esl-alimerge --outformat pfam --dna -o ' + query + '.' + ref + '.clean.align.sto ' + query + '.clean.align.sto ' + variables['ref_dir'] + ref + '.refpkg/' + ref + '.clean.align.sto', shell = True, executable = executable)
    combine.communicate()
      
    convert = subprocess.Popen('seqmagick convert ' + query + '.' + ref + '.clean.align.sto ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()
    
    pplacer = subprocess.Popen('pplacer -o ' + query + '.' + ref + '.clean.align.jplace --out-dir ' + os.getcwd() + ' -p --keep-at-most 20 -c ' + variables['ref_dir'] + ref + '.refpkg ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    pplacer.communicate()
    
## Define function to generate csv file of placements and fat tree    
    
def guppy(query, ref):
    
    guppy1 = subprocess.Popen('guppy to_csv --point-mass --pp -o ' + query + '.' + ref + '.clean.align.csv ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy1.communicate()
    
    guppy2 = subprocess.Popen('guppy fat --node-numbers --point-mass --pp -o ' + query + '.' + ref + '.clean.align.phyloxml ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy2.communicate()   
                
## Execute main program.
                
## First option is for testing purposes, allows for testing from with Python.
                
if len(sys.argv) == 1:
        
    query = 'test'
    ref = 'combined_16S.tax'
    
    clear_wspace = subprocess.call('rm ' + query + '.' + ref + '*', shell = True, executable = executable)        

    place(query, ref, variables)

elif 'query' not in command_args.keys():
    
    ## Build reference package
    
    ref = command_args['ref']        
    clean_name(variables['ref_dir'] + ref)
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + variables['ref_dir'] + ref + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()

    infernal_commands = 'cmalign --dna -o ' + variables['ref_dir'] + ref + '.clean.align.sto --outformat Pfam ' + variables['cm'] + ' ' + variables['ref_dir'] + ref + '.clean.fasta'      
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + variables['ref_dir'] + ref + '.clean.align.sto ' + variables['ref_dir'] + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()     
    
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
    
    ## Copy sto of the reference alignment to ref package.
    
    cp = subprocess.Popen('cp ' + variables['ref_dir'] + ref + '.clean.align.sto ' + variables['ref_dir'] + ref + '.refpkg/' + ref + '.clean.align.sto', shell = True, executable = executable)

    ## Create a file with the date/time of database creation.

    current_time = datetime.datetime.now().isoformat()
    n_aseqs = len(re.findall('>', open(variables['ref_dir'] + '/' + ref + '.fasta', 'r').read()))
    
    with open(variables['ref_dir'] + '/' + ref + '.database_info.txt', 'w') as database_info:
        print >> database_info, 'ref tree built at:', current_time
        print >> database_info, 'nseqs in reference alignment:', n_aseqs 
            
else:
    
    query = command_args['query']
    ref = command_args['ref']
    splits = int(command_args['splits'])
    
    clear_wspace = subprocess.call('rm ' + query + '.' + ref + '*', shell = True, executable = executable)
    
    ## Select a random subset of reads, if this option is specified.  This is useful for
    ## normalizing the number of sampled reads across different samples.    
    
    if 'n' in command_args.keys():
        
        nseqs = command_args['n']
        tseqs = len(re.findall('>', open(query + '.fasta', 'r').read()))
        nseqs_get = random.sample(range(1, tseqs), int(nseqs))
        
        seq_i = 0

        with open(query + '.sub.fasta', 'w') as fasta_sub:
            for record in SeqIO.parse(query + '.fasta', 'fasta'):
                seq_i = seq_i + 1
                if seq_i in nseqs_get:
                    SeqIO.write(record, fasta_sub, 'fasta')
        
        ## All downstream operations now need to take place on the subsamled
        ## query, easiest way to do this is to just change the query variable
            
        query = query + '.sub'
                        
    ## Create splits, if splits > 1
    
    if splits > 1:        
        split_list = split_fasta(query, splits)
            
        if __name__ == '__main__':  
            Parallel(n_jobs = splits, verbose = 5)(delayed(place)
            (split_query, ref, variables) for split_query in split_list)
            
        guppy_merge = subprocess.Popen('guppy merge ' + query + '*' + '.jplace -o ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
        guppy_merge.communicate()
        guppy(query, ref)
        
        cleanup = subprocess.Popen('rm ' + query + '.temp*', shell = True, executable = executable)
        cleanup.communicate()
        
    else:
        place(query, ref, variables)
        guppy(query, ref)
