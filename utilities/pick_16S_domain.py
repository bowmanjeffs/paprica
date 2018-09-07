# -*- coding: utf-8 -*-
"""
Created on Thu Feb 04 09:21:16 2016

@author: jeff

Call as:
    python pick_16S_domain.py -a_in [archaea input prefix] -b_in [bacteria input prefix] -e_in [eukaryote input prefix] -out [output prefix]
    
Input prefix must be the same for input files from each domain, such as:

PAL_0007_bacterial_16S_hits.txt
PAL_0007_archaeal_16S_hits.txt
PAL_0007_eukaryote_16S_hits.txt

PAL_0007_bacterial_16S_hits.sto
PAL_0007_archaeal_16S_hits.sto
PAL_0007_eukaryote_16S_hits.sto

"""

import pandas as pd
from Bio import SeqIO
import sys

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]

try:        
    prefix = command_args['prefix']
    output = command_args['out']
except KeyError:
    prefix = 'SRR7156747.unique.fasta'
    output = 'SRR7156747.unique.fasta'

colnames = ['target.name', 'target.accession', 'query.name', 'accession', 'mdl', 'mdl.from', 'mdl.to', 'seq.from', 'seq.to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description']

cmscan = pd.read_csv(prefix + '.txt', comment = '#', names = colnames, header = None, skiprows = [0,1], delim_whitespace = True, index_col = 2)

bacteria_set = []
archaea_set = []
eukarya_set = []

## need to use groupby

with open(output + '.bacterial16S.reads.txt', 'w') as bacteria_out, open(output + '.archaeal16S.reads.txt', 'w') as archaea_out, open(output + '.eukaryote18S.reads.txt', 'w') as eukarya_out:
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
            
def write_fasta(fasta_out, prefix, domain_set):
    
    with open(fasta_out, 'w') as fasta_out:            
        for record in SeqIO.parse(prefix, 'fasta'):            
            if record.id in domain_set:
                SeqIO.write(record, fasta_out, 'fasta')
                
write_fasta(output + '.bacterial16S.fasta', prefix, bacteria_set)
write_fasta(output + '.archaeal16S.fasta', prefix, archaea_set)
write_fasta(output + '.eukaryote18S.fasta', prefix, eukarya_set)
