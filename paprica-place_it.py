#!/usr/bin/env python
# -*- coding: utf-8 -*-

help_string = """
Created on Sat Jan 03 08:59:11 2015

@author: Jeff Bowman, bowmanjs@ldeo.columbia.edu

paprica is licensed under a Creative Commons Attribution-NonCommercial
4.0 International License.  IF you use any portion of paprica in your
work please cite:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described
by Metabolic Structure: A General Framework and Application to a Seasonally
Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic
Peninsula." PloS one 10.8 (2015): e0135868.

If your analysis makes specific use of pplacer, Infernal, or pathway-tools
please make sure that you also cite the relevant publications.

REQUIRES:
    Files:
        bacterial_ssu.cm
        
    Programs:
        infernal
        taxtastic
        seqmagick
        pplacer
        raxmlHPC-PTHREADS-AVX2
        
    Python modules:
        Bio
        joblib

CALL AS:
    python genome_finder_place_it.py -query [query] -ref [ref] -splits [splits] -n [nseqs] for analysis or
    python genome_finder_place_it.py -ref [ref] -cpus [ncpus] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (which must be .fasta).
    
    It is not recommended to use more than 8 cpus for tree building (i.e. -cpus 8).
    See the RAXML manual for guidance on this.
    
OPTIONS:
    -cpus:  The number of cpus for RAxML to use.
    -domain: The domain (bacteria or archaea) you are analyzing for.
    -n:  The number of reads to subsample your input query to.
    -query: The prefix of the query fasta (entire name without extension).
    -ref: The name of the reference package for pplacer.
    -ref_dir: The directory containing the paprica database.
    -splits: The number of files to split the query fasta into to facilitate
        parallel analysis with pplacer.

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
executable = '/bin/bash' # shell for executing commands

import re
import subprocess
import sys
import os
from joblib import Parallel, delayed
import datetime
import random

paprica_path = os.path.dirname(os.path.abspath(__file__)) + '/' # The location of the actual paprica scripts.
cwd = os.getcwd() + '/' # The current working directory.
                
## Parse command line arguments.  Arguments that are unique to a run,
## such as the number of splits, should be specified in the command line and
## not in paprica_profile.txt
                
command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in command_args.keys():
    print help_string
    quit()

## Parse command line for mandatory variables, providing defaults if they are
## not present.  Defaults are generally for testing.

try:
    ref_dir = paprica_path + command_args['ref_dir']  # The complete path to the reference directory being used for analysis.        
    domain = command_args['domain']  # The domain being used for analysis.
    ref = command_args['ref']  # The name of the reference package being used.
    
except KeyError:
    ref_dir = paprica_path + 'ref_genome_database/'
    domain = 'bacteria'
    ref = 'combined_16S.bacteria.tax'
        
## Parse command line for optional variables, providing defaults if they are not
## present.
    
try:
    cpus = str(command_args['cpus']) # The number of cpus for RAxML to use.
except KeyError:
    cpus = '1'
    
try:    
    splits = int(command_args['splits']) # The number of splits to make for pplacer parallel operation.
except KeyError:
    splits = 1
    
try:
    query = command_args['query'] # The prefix of the query file.
except KeyError:
    query = 'test.bacteria'
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
ref_dir_domain = ref_dir + domain + '/' # Complete path to the domain subdirectory within the reference directory.
cm = paprica_path + 'models/' + domain + '_ssu.cm' # Domain specific covariance model used by infernal.
    
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
    
def place(query, ref, ref_dir_domain, cm):
            
    clean_name(query)
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + query + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    ## Conduct alignment with Infernal (cmalign) against the bacteria profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    
    align = subprocess.Popen('cmalign --dna -o ' + query + '.clean.align.sto --outformat Pfam ' + cm + ' ' + query + '.clean.fasta', shell = True, executable = executable)
    align.communicate()    
    
    combine = subprocess.Popen('esl-alimerge --outformat pfam --dna -o ' + query + '.' + ref + '.clean.align.sto ' + query + '.clean.align.sto ' + ref_dir_domain + ref + '.refpkg/' + ref + '.clean.align.sto', shell = True, executable = executable)
    combine.communicate()
      
    convert = subprocess.Popen('seqmagick convert ' + query + '.' + ref + '.clean.align.sto ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()
    
    pplacer = subprocess.Popen('pplacer -o ' + query + '.' + ref + '.clean.align.jplace --out-dir ' + os.getcwd() + ' -p --keep-at-most 20 -c ' + ref_dir_domain + ref + '.refpkg ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
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
    
    clear_wspace = subprocess.call('rm -f ' + cwd + query + '.' + ref + '*', shell = True, executable = executable)        
    place(cwd + query, ref, ref_dir_domain, cm)

elif 'query' not in command_args.keys():
    
    ## Build reference package
    
    clean_name(ref_dir_domain + ref)
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + ref_dir_domain + ref + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()

    infernal_commands = 'cmalign --dna -o ' + ref_dir_domain + ref + '.clean.align.sto --outformat Pfam ' + cm + ' ' + ref_dir_domain + ref + '.clean.fasta'      
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + ref_dir_domain + ref + '.clean.align.sto ' + ref_dir_domain + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()     
    
    rm = subprocess.call('rm -f ' + ref_dir_domain + '*ref.tre', shell = True, executable = executable)
    raxml1 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -m GTRGAMMA -s ' + ref_dir_domain + ref + '.clean.align.fasta -n ref.tre -f d -p 12345 -w ' + ref_dir_domain, shell = True, executable = executable)
    raxml1.communicate()
    
    ## Root the tree.
    
    raxml2 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T 2 -m GTRGAMMA -f I -t ' + ref_dir_domain + 'RAxML_bestTree.ref.tre -n root.ref.tre -w ' + ref_dir_domain, shell = True, executable = executable)  
    raxml2.communicate()
    
    ## Generate SH-like support values for the tree.
    
    raxml3 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -m GTRGAMMA -f J -p 12345 -t ' + ref_dir_domain + 'RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s ' + ref_dir_domain + ref + '.clean.align.fasta -w ' + ref_dir_domain, shell = True, executable = executable)   
    raxml3.communicate()
     
    ## Generate the reference package using the tree with SH support values and a log file.
    
    rm = subprocess.call('rm -rf ' + ref_dir_domain + ref + '.refpkg', shell = True, executable = executable)
    taxit = subprocess.Popen('taxit create -l 16S_rRNA -P ' + ref_dir_domain + ref + '.refpkg --aln-fasta ' + ref_dir_domain + ref + '.clean.align.fasta --tree-stats ' + ref_dir_domain + 'RAxML_info.ref.tre --tree-file ' + ref_dir_domain + 'RAxML_fastTreeSH_Support.conf.root.ref.tre', shell = True, executable = executable)
    taxit.communicate()
    
    ## Copy sto of the reference alignment to ref package.
    
    cp = subprocess.Popen('cp ' + ref_dir_domain + ref + '.clean.align.sto ' + ref_dir_domain + ref + '.refpkg/' + ref + '.clean.align.sto', shell = True, executable = executable)

    ## Create a file with the date/time of database creation.

    current_time = datetime.datetime.now().isoformat()
    n_aseqs = len(re.findall('>', open(ref_dir_domain + '/' + ref + '.fasta', 'r').read()))
    
    with open(ref_dir_domain + '/' + ref + '.database_info.txt', 'w') as database_info:
        print >> database_info, 'ref tree built at:', current_time
        print >> database_info, 'nseqs in reference alignment:', n_aseqs 
            
else:
    
    splits = int(command_args['splits'])
    
    clear_wspace = subprocess.call('rm -f ' + cwd + query + '.' + ref + '*', shell = True, executable = executable)
    clear_wspace = subprocess.call('rm -f ' + cwd + query + '.sub*', shell = True, executable = executable)
    
    ## Select a random subset of reads, if this option is specified.  This is useful for
    ## normalizing the number of sampled reads across different samples.    
    
    if 'n' in command_args.keys():
        
        nseqs = command_args['n']
        tseqs = len(re.findall('>', open(cwd + query + '.fasta', 'r').read()))
        nseqs_get = random.sample(range(1, tseqs), int(nseqs))
        
        seq_i = 0

        with open(cwd + query + '.sub.fasta', 'w') as fasta_sub:
            for record in SeqIO.parse(cwd + query + '.fasta', 'fasta'):
                seq_i = seq_i + 1
                if seq_i in nseqs_get:
                    SeqIO.write(record, fasta_sub, 'fasta')
        
        ## All downstream operations now need to take place on the subsamled
        ## query, easiest way to do this is to just change the query variable
            
        query = query + '.sub'
                        
    ## Create splits, if splits > 1
    
    if splits > 1:        
        split_list = split_fasta(cwd + query, splits)
            
        if __name__ == '__main__':  
            Parallel(n_jobs = splits, verbose = 5)(delayed(place)
            (split_query, ref, ref_dir_domain, cm) for split_query in split_list)
            
        guppy_merge = subprocess.Popen('guppy merge ' + cwd + query + '*' + '.jplace -o ' + cwd + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
        guppy_merge.communicate()
        guppy(cwd + query, ref)
        
        cleanup = subprocess.Popen('rm -f ' + cwd + query + '.temp*', shell = True, executable = executable)
        cleanup.communicate()
        
    else:
        place(cwd + query, ref, ref_dir_domain, cm)
        guppy(cwd + query, ref)
