#!/usr/bin/env python
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
    python pick_16S_domain.py -in [sample prefix]

    Note that [prefix] includes the entire file name without extension
    (which must be .fasta).
    
OPTIONS:
    -in:  The prefix of the query fasta (entire name without extension).
    -h:  Print this help message.

This script must be located in the 'paprica' directory as it makes use of relative
paths.


Call as:
    pick_16S_domain.py -in [sample prefix]

"""

import pandas as pd
from Bio import SeqIO
import sys
import os

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]
        
if 'h' in command_args.keys():
    print help_string
    quit()

try:        
    prefix = command_args['in']
except KeyError:
    prefix = 'test.bacteria'
    
fasta_in = prefix + '.fasta'
    
try:
    paprica_path = os.path.dirname(os.path.realpath(__file__)) + '/' # The location of the actual paprica scripts.
except NameError:
    paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/'
cwd = os.getcwd() + '/' # The current working directory.

## Define a function to search an input fasta against a set, and write
## reads with ids in that set to an output fasta.

def write_fasta(fasta_out, fasta_in, domain_set):
    
    with open(fasta_out, 'w') as fasta_out:            
        for record in SeqIO.parse(fasta_in, 'fasta'):            
            if record.id in domain_set:
                SeqIO.write(record, fasta_out, 'fasta')

## Search the input fasta with cmscan against the covariance model of all domains.
       
os.system('cmscan --tblout ' + cwd + prefix + '.txt ' + paprica_path + 'models/all_domains.cm ' + cwd + fasta_in + ' > /dev/null')

## Read in the output of cmscan, then iterate across all lines, selecting for each read the domain with the lowest 
## E-value.  This is considered to the be true domain of the genome originating the read.
    
colnames = ['target.name', 'target.accession', 'query.name', 'accession', 'mdl', 'mdl.from', 'mdl.to', 'seq.from', 'seq.to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description']

cmscan = pd.read_csv(cwd + prefix + '.txt', comment = '#', names = colnames, header = None, skiprows = [0,1], delim_whitespace = True, index_col = 2)

bacteria_set = []
archaea_set = []
eukarya_set = []

with open(cwd + prefix + '.bacterial16S.reads.txt', 'w') as bacteria_out, open(cwd + prefix + '.archaeal16S.reads.txt', 'w') as archaea_out, open(cwd + prefix + '.eukaryote18S.reads.txt', 'w') as eukarya_out:
    for index in set(cmscan.index):        
        temp = cmscan.loc[index]
        
        try:
            e_min = temp[temp['E-value'] == temp['E-value'].min()]
            domain = e_min.loc[index, 'target.name']

        except AttributeError:
            domain = temp['target.name']           
                    
        if 'bacteria' in str(domain):
            print >> bacteria_out, index
            bacteria_set.append(index)
        elif 'archaea' in str(domain):
            print >> archaea_out, index 
            archaea_set.append(index)
        elif 'eukarya' in str(domain):
            print >> eukarya_out, index
            eukarya_set.append(index)
            
bacteria_set = set(bacteria_set)
archaea_set = set(archaea_set)
eukarya_set = set(eukarya_set)
            
print prefix, len(bacteria_set), len(archaea_set), len(eukarya_set)

## Write fasta files for each domain.
                
write_fasta(prefix + '.bacteria.fasta', fasta_in, bacteria_set)
write_fasta(prefix + '.archaea.fasta', fasta_in, archaea_set)
write_fasta(prefix + '.eukarya.fasta', fasta_in, eukarya_set)