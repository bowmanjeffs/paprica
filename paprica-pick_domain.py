#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Thu Feb 04 09:21:16 2016

@author: jeff

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
        models/all_domains.cm
        models/all_domains.cm.i1f
        models/all_domains.cm.i1i
        models/all_domains.cm.i1p
        models/all_domains.cm.i1m
        
    Programs:
        infernal
        
    Python modules:
        Bio
        pandas

CALL AS:
    python paprica-pick_domain.py -in [sample prefix] -min_score [min bit score]

    Note that [prefix] includes the entire file name without extension
    (which must be .fasta).
    
OPTIONS:
    -in:  The prefix of the query fasta (entire name without extension).
    -h:  Print this help message.
    -min_score: The minimum bit score to keep a read.

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""

import pandas as pd
from Bio import SeqIO
import sys
import os
import re
from joblib import Parallel, delayed
import multiprocessing

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]
        
if 'h' in list(command_args.keys()):
    print(help_string)
    quit()

try:        
    prefix = command_args['in']
except KeyError:
    prefix = 'test'
try:
    min_score = command_args['min_score']
except KeyError:
    min_score = 1
    
fasta_in = prefix + '.fasta'
    
try:
    paprica_path = os.path.dirname(os.path.realpath(__file__)) + '/' # The location of the actual paprica scripts.
except NameError:
    paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/'
cwd = os.getcwd() + '/' # The current working directory.

#%% Check to make sure the output files don't already exist.  Exit if they do.

file_set = set(os.listdir(cwd))
exit_test = 0

if prefix + '.bacteria.fasta' in file_set:
    exit_test = exit_test + 1
    
if prefix + '.archaea.fasta' in file_set:
    exit_test = exit_test + 1
    
if prefix + '.eukarya.fasta' in file_set:
    exit_test = exit_test + 1

if exit_test == 3:
    exit()
    
#%% Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])

#%% Define a function to search an input fasta against a set, and write
## reads with ids in that set to an output fasta.

def write_fasta(fasta_out, fasta_in, domain_set, domain):
    
    n = 0
    
    with open(fasta_out, 'w') as fasta_out:            
        for record in SeqIO.parse(fasta_in, 'fasta'):            
            if record.id in domain_set:
                
                temp_record = record
                temp_names = seq_names[str(record.seq)]
                for name in temp_names:
                    n = n + 1
                    temp_record.id = name
                    temp_record.description = ''
                    SeqIO.write(record, fasta_out, 'fasta')
                    
    print(domain, 'total:unique', str(n) + ':' + str(len(domain_set)))
                
#%% Define a function to split an input fasta file.
                
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
    
#%% Define a function to run cmscan

def run_cmscan(split_in):   
    os.system('cmscan --cpu 1 --tblout ' + split_in + '.txt ' + paprica_path + 'models/all_domains.cm ' + split_in + '.fasta > /dev/null')
    
#%% Define a function to create a unique fasta file from the input.
            
def make_unique(query):
    
    seq_count = {}
    seq_names = {}
    name_seq = {}
        
    for record in SeqIO.parse(query + '.fasta', 'fasta'):
        name = str(record.id)
        seq = str(record.seq)
        
        if seq in list(seq_count.keys()):
            seq_count[seq] = seq_count[seq] + 1
            temp_name_list = seq_names[seq]
            temp_name_list.append(name)
            seq_names[seq] = temp_name_list
        else:
            seq_count[seq] = 1
            seq_names[seq] = [name]
            
    with open(query + '.unique.fasta', 'w') as fasta_out, open(query + '.unique.count', 'w') as count_out:
        print('rep_name' + ',' + 'abundance', file = count_out)
        for seq in list(seq_names.keys()):
            print('>' + seq_names[seq][0], file = fasta_out)
            print(seq, file = fasta_out)
            
            print(seq_names[seq][0] + ',' + str(seq_count[seq]), file = count_out)
            name_seq[seq_names[seq][0]] = seq
            
    return seq_count, seq_names, name_seq
    
#%% Run program
    
## Make unique.  If program hangs here you probably didn't denoise your reads
## adequately!

seq_count, seq_names, name_seq = make_unique(cwd + prefix)

## Determine number of processors, but don't let the number of splits
## exceed the number of sequences.
    
nsplits = multiprocessing.cpu_count()

if nsplits > len(name_seq):
    nsplits = len(name_seq)
    
split_list = split_fasta(cwd + prefix + '.unique', nsplits)

## Search the input fasta with cmscan against the covariance model of all domains.
    
if __name__ == '__main__':  
    Parallel(n_jobs = nsplits, verbose = 5)(delayed(run_cmscan)
    (split_in) for split_in in split_list)
    
## Aggregate output.
    
colnames = ['target.name', 'target.accession', 'query.name', 'accession', 'mdl', 'mdl.from', 'mdl.to', 'seq.from', 'seq.to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description']
cmscan = pd.DataFrame(columns = colnames)

for split in split_list:
    try:
        cmscan_temp = pd.read_csv(split + '.txt', comment = '#', names = colnames, header = None, skiprows = [0,1], delim_whitespace = True, index_col = 2)
        cmscan = pd.concat([cmscan, cmscan_temp], sort = False)
    except IOError:
        continue
    
## Cleanup temp output.
        
os.system('rm -f ' + cwd + prefix + '.unique.temp*')
        
## Iterate across all lines, selecting for each read the domain with the lowest 
## E-value.  This is considered to the be true domain of the genome originating the read.
    
bacteria_set = []
archaea_set = []
eukarya_set = []

with open(cwd + prefix + '.bacterial16S.reads.txt', 'w') as bacteria_out, open(cwd + prefix + '.archaeal16S.reads.txt', 'w') as archaea_out, open(cwd + prefix + '.eukaryote18S.reads.txt', 'w') as eukarya_out:
    for index in set(cmscan.index):        
        temp = cmscan.loc[index]
        
        try:
            e_min = temp[temp['score'] == temp['score'].max()]
            domain = e_min.loc[index, 'target.name']
            
        ## If read returned a hit for only one domain you
        ## will have a series and not a dataframe.  For some reason different
        ## pandas versions return a different error, need to pass both here.

        except (AttributeError, KeyError):
            e_min = temp
            domain = temp['target.name']  
            
        ## If the bit score exceeds the specified minimum add to appropriate
        ## domain list.
            
        if float(e_min['score']) > min_score:
            if 'bacteria' in str(domain):
                print(index, file=bacteria_out)
                bacteria_set.append(index)
            elif 'archaea' in str(domain):
                print(index, file=archaea_out) 
                archaea_set.append(index)
            elif 'eukarya' in str(domain):
                print(index, file=eukarya_out)
                eukarya_set.append(index)
            
bacteria_set = set(bacteria_set)
archaea_set = set(archaea_set)
eukarya_set = set(eukarya_set)
            
## Write fasta files for each domain.
                
write_fasta(prefix + '.bacteria.fasta', prefix + '.unique.fasta', bacteria_set, 'bacteria')
write_fasta(prefix + '.archaea.fasta', prefix + '.unique.fasta', archaea_set, 'archaea')
write_fasta(prefix + '.eukarya.fasta', prefix + '.unique.fasta', eukarya_set, 'eukarya')
