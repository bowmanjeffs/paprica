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
    python genome_finder_place_it.py -domain [domain] -query [query] -ref [ref] -splits [splits] -n [nseqs] for analysis or
    python genome_finder_place_it.py -domain [domain] -ref [ref] -cpus [ncpus] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (which must be .fasta).
    
    It is not recommended to use more than 8 cpus for tree building (i.e. -cpus 8).
    See the RAXML manual for guidance on this.
    
OPTIONS:
    -o: out directory
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
import pandas as pd
import numpy as np
from Bio import SeqIO

try:
    paprica_path = os.path.dirname(os.path.abspath(__file__)) + '/' # The location of the actual paprica scripts.
except NameError:
    paprica_path = os.path.dirname(os.path.abspath("__file__")) + '/'
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
except KeyError:
    ref_dir = paprica_path + 'ref_genome_database/'
try:     
    domain = command_args['domain']  # The domain being used for analysis.
except KeyError:
    domain = 'bacteria'
try:
    ref = command_args['ref']  # The name of the reference package being used.
except KeyError:
    ref = 'taxit'    
try:
    out = command_args['o']  # The name of the reference package being used.
except KeyError:
    out = 'outFolder/placeIt/'
    print 'Output will be stored in '+cwd+out
try :
    os.makedirs(out)
except OSError:
	pass
print "Place-it output being stored in = ", out

## If sys.argv == 1, you are probably running inside Python in testing mode.
## Provide some default values to make this possibe.  If > 1, parse command
## line for optional variables, providing reasonable defaults if not present.
## No default is currently provided for query.
    
if len(sys.argv) == 1:
    cpus = '8'
    #splits = 1
    #query = 'test.eukarya'

else:
    
    try:
        cpus = str(command_args['cpus']) # The number of cpus for RAxML to use.
    except KeyError:
        cpus = '1'
    
    try:    
        splits = int(command_args['splits']) # The number of splits to make for pplacer parallel operation.
    except KeyError:
        splits = 1
        
    try:
        query = command_args['query']
    except KeyError:
        pass
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'


if out.endswith('/') == False:
    out = out + '/'
    
ref_dir_domain = ref_dir + domain + '/' # Complete path to the domain subdirectory within the reference directory.
cm = paprica_path + 'models/' + domain + '_ssu.cm' # Domain specific covariance model used by infernal.

#%% Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]
    
#%% Define function to clean record names.  Not that bad_character is also used
## by make_tax, to insure that sequence names match between taxonomy database
## and tree.
    
bad_character = '[\|\\=-@!%,;\(\):\'\"\s]'



def clean_name(file_name, bad_character):

    bad_character = re.compile(bad_character)
    file_out=out+os.path.basename(file_name)+'.clean.fasta'

    fasta_out=open(file_out,'w+')
    for record in SeqIO.parse(file_name + '.fasta', 'fasta'): 
        record.name = re.sub(bad_character, '_', str(record.description))
        print >> fasta_out, '>' + record.name 
        print >> fasta_out, record.seq 
            
#%% Define function to split fasta file to run pplacer in parallel.  This greatly
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

#%% Define function to execute phylogenetic placement on a query fasta, or split query fasta
    
def place(query, ref, ref_dir_domain, cm):

    clean_name(query, bad_character)
    query=out+os.path.basename(query)

    degap = subprocess.Popen('seqmagick mogrify --ungap ' + query + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    ## Conduct alignment with Infernal (cmalign) against the bacteria profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    print query
    align = subprocess.Popen('cmalign --dna -o ' + query + '.clean.align.sto --outformat Pfam --mxsize 3000 --tfile ' +query+'.traceback.txt --sfile ' +query+'.table.txt '+ cm + ' ' + query + '.clean.fasta', shell = True, executable = executable)
    align.communicate()    
    
    combine = subprocess.Popen('esl-alimerge --outformat pfam --dna \
    -o ' + query + '.' + ref + '.clean.align.sto \
    ' + query + '.clean.align.sto \
    ' + ref_dir_domain + ref + '/'+domain+'.align.DNA.sto', shell = True, executable = executable)
    combine.communicate()
      
    convert = subprocess.Popen('seqmagick convert ' + query + '.' + ref + '.clean.align.sto ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()
    
    pplacer = subprocess.Popen('pplacer -o ' + query + '.' + ref + '.clean.align.jplace \
    --out-dir ' + out + ' \
    -p --keep-at-most 20 --map-identity\
    -c ' + ref_dir_domain + ref + ' ' + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
    pplacer.communicate()
    
#%% Define function to generate csv file of placements and fat tree    
    
def guppy(query, ref):
    query=out+os.path.basename(query)
    guppy1 = subprocess.Popen('guppy to_csv --point-mass --pp -o ' + query + '.' + ref + '.clean.align.csv ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy1.communicate()
    
    guppy2 = subprocess.Popen('guppy fat --node-numbers --point-mass --pp -o ' + query + '.' + ref + '.clean.align.phyloxml ' + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
    guppy2.communicate()
    
#%% Define a function to classify read placements.
    
def classify():
    
    ## Prep sqlite database for classification.  Would be great not to have to do
    ## this, but classification command requires sqlite output.
    
    prep_database = subprocess.Popen('rppr prep_db -c ' + ref_dir_domain + ref + '.refpkg \
    --sqlite ' + query + '.' + ref + '.clean.align.db',
    shell = True,
    executable = executable)
    
    prep_database.communicate()
    
    ## Now execute guppy classification command.
    
    guppy_classify = subprocess.Popen('guppy classify \
    -c ' + ref_dir_domain + ref + '.refpkg \
    --sqlite ' + query + '.' + ref + '.clean.align.db ' + query + '.' + ref + '.clean.align.jplace',
    shell = True,
    executable = executable)
    
    guppy_classify.communicate()    
    


#%% Define a function to tally the unique reads assigned to each edge.

def count_unique():

    query_seqs = pd.DataFrame()
    
    ## Add all sequences in aligned fasta to the dataframe, indexed by record id.
    
    for record in SeqIO.parse(out+os.path.basename(query) + '.' + ref + '.clean.align.fasta', 'fasta'):
        query_seqs.loc[str(record.id), 'sequence'] = str(record.seq)

    	
    ## Read in the guppy csv file, indexing by edge number.
    query_csv = pd.DataFrame.from_csv(out+os.path.basename(query) + '.' + ref + '.clean.align.csv', header = 0, index_col = 3)
    
    ## Iterating by edge number, get all the sequences for that edge, then count the number of unique sequences.
    
    unique_edge_counts = pd.DataFrame(columns = ['rep', 'abundance', 'edge_num'])
    for edge_num in query_csv.index.unique():
        print 'Finding unique sequences within edge', edge_num
        seqs=query_csv.loc[edge_num,'name']
        if np.size(seqs)==1:
            seqs=str(seqs)
        else:
            seqs= [str(x) for x in query_csv.loc[edge_num, 'name']]
        temp_df = query_seqs.loc[seqs]
        temp_seq_names = []
        try:
            temp_counts = temp_df.sequence.value_counts()
    		## Cludgy, but I can't find another way to capture unique sequence names!

            for seq in temp_counts.index:
                seq_name = temp_df.sequence[temp_df.sequence == seq].index.tolist()[0]
                temp_seq_names.append(seq_name)

        except AttributeError:
            temp_counts = pd.Series(1, name = temp_df.sequence)
            temp_counts.index = [temp_counts.name]
            temp_seq_names.append(temp_df.name)
        temp_df_out = pd.DataFrame(temp_counts)
        temp_df_out.columns = ['abundance']
        temp_df_out['edge_num'] = int(edge_num)
        temp_df_out['rep'] = temp_seq_names
    			   
        ## Hash the sequence to create a lighter-weight index that reflects the
        ## sequence.  This will be used for comparison across samples.
    	
        temp_df_out['seq'] = temp_counts.index                   
        temp_df_out.index = temp_df_out['seq'].apply(hash)
        temp_df_out.drop('seq', axis = 1, inplace = True)
        unique_edge_counts = pd.concat([unique_edge_counts, temp_df_out])

    unique_edge_counts.edge_num = unique_edge_counts.edge_num.astype(int)
    return(unique_edge_counts)
    
#%% Execute main program.
    
clear_wspace = subprocess.call('rm -f ' + query + '.' + ref + '*', shell = True, executable = executable)
clear_wspace = subprocess.call('rm -f ' + query + '.sub*', shell = True, executable = executable)
    
## Select a random subset of reads, if this option is specified.  This is useful for
## normalizing the number of sampled reads across different samples.    
    
if 'n' in command_args.keys():
        
    nseqs = command_args['n']
    tseqs = len(re.findall('>', open(cwd + query + '.fasta', 'r').read()))
    nseqs_get = random.sample(range(1, tseqs), int(nseqs))
        
    seq_i = 0

    with open(query + '.sub.fasta', 'w') as fasta_sub:
        for record in SeqIO.parse(cwd + query + '.fasta', 'fasta'):
            seq_i = seq_i + 1
            if seq_i in nseqs_get:
                SeqIO.write(record, fasta_sub, 'fasta')
        
        ## All downstream operations now need to take place on the subsampled
        ## query, easiest way to do this is to just change the query variable.
            
    query = query + '.sub'
                        
## Create splits, if splits > 1
    
if splits > 1:   
    query=out+os.path.basename(query)
    split_list = split_fasta( query, splits)
            
    if __name__ == '__main__':  
        Parallel(n_jobs = splits, verbose = 5)(delayed(place)
        (split_query, ref, ref_dir_domain, cm) for split_query in split_list)
            
        guppy_merge = subprocess.Popen('guppy merge ' + query + '*' + '.jplace -o ' + out + query + '.' + ref + '.clean.align.jplace', shell = True, executable = executable)
        guppy_merge.communicate()
        guppy(out + query, ref)
        
        merge = subprocess.Popen('cat ' + cwd + query + '.temp*.' + ref + '.clean.align.fasta > ' + out + query + '.' + ref + '.clean.align.fasta', shell = True, executable = executable)
        merge.communicate()
                
        cleanup = subprocess.Popen('rm -f ' + query + '.temp*', shell = True, executable = executable)
        cleanup.communicate()
        
else:

    place(query, ref, ref_dir_domain, cm)
    guppy(query, ref)
    unique_seqs = count_unique()
    unique_seqs.to_csv(out+os.path.basename(query) + '.unique.seqs.csv')
	#classify()

