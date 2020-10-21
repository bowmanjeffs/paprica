#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Nov 28 15:45:43 2015

@author: jeff

This script aggregates information from multiple '.edge_data.csv' files
produces by running paprica on multiple samples.  It produces a matrix of edges
by sample, and a matrix of mean edge parameters, by sample.

For simple execution (works in most user cases) run as:
    ./combine_edge_results.py -domain [domain] -o [prefix for output] -d [analysis directory]
    
If your file names don't follow the default suffix pattern of (e.g.)
[domain].edge_data.csv, you will need to specify the suffix pattern like this:
    ./combine_edge_results.py -domain [domain]  -o [prefix for output] -edge_in [suffix pattern for edges] -path_in [suffix pattern for paths] -ec_in [suffix pattern for ec numbers] -unique_in [suffix pattern for unique sequences]

It will automatically loop through all files in the directory with the specified suffixes.

"""

import os
import pandas as pd
import re
import math
import sys

## Read in command line arguments.

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in list(command_args.keys()):
    pass
    #print help_string ## no help sting
    quit()
    
try:
    domain = command_args['domain']
except KeyError:
    domain = 'archaea'
    
try:
    prefix = command_args['o']
except KeyError:
    prefix = 'test.' + domain
    
try:
    unique_suffix = command_args['unique_in']
except KeyError:
    unique_suffix = domain + '.unique_seqs.csv'

try:
    edge_suffix = command_args['edge_in']
except KeyError:
    edge_suffix = domain + '.edge_data.csv'
    
if domain != 'eukarya':
    
    try:
        path_suffix = command_args['path_in']
    except KeyError:
        path_suffix = domain + '.sum_pathways.csv'
        
    try:
        ec_suffix = command_args['ec_in']
    except KeyError:
        ec_suffix = domain + '.sum_ec.csv'
        
cwd = os.getcwd() + '/'
prefix = cwd + prefix
            
## Delete old combined files, so that
## the script doesn't try to include them.

os.system('rm -f ' + prefix + '.' + domain + '.edge_tally.csv')
os.system('rm -f ' + prefix + '.' + domain + '.edge_data.csv')
os.system('rm -f ' + prefix + '.' + domain + '.path_tally.csv')
os.system('rm -f ' + prefix + '.' + domain + '.ec_tally.csv')
    
def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])

edge_tally = pd.DataFrame()
edge_data = pd.DataFrame()
ec_tally = pd.DataFrame()
path_tally = pd.DataFrame()

def fill_edge_data(param, name, df_in):
    temp = []
    for index in df_in.index:
        
        ## Exception is necessary for old version of build_core_genomes which
        ## did not provide a n16S value for draft genomes.
        
        try:
            n = list(range(int(math.ceil(df_in.loc[index, 'nedge_corrected']))))
        except ValueError:
            n = list(range(int(df_in.loc[index, 'nedge'])))
            
        for i in n:
            temp.append(df_in.loc[index, param])
            
    temp = pd.Series(temp)
            
    mean = temp.mean()
    sd = temp.std()
    
    print(name, param, mean, sd)
    return mean, sd

if domain == 'eukarya':
    ranks = ['kingdom', 'supergroup', 'division', 'class', 'order', 'family',
             'genus', 'species', 'taxon']
else:
    ranks = ['superkingdom', 'phylum', 'clade', 'class',
         'order', 'family', 'genus', 'species', 'strain', 'taxon']

taxon_map = pd.DataFrame(columns = ranks)
        
for f in os.listdir(cwd):
    if f.endswith(edge_suffix):
        
        temp_edge = pd.read_csv(cwd + f, index_col = 0)
        name = re.sub(edge_suffix, '', f)
        
        for edge in temp_edge.index:            
            for rank in ranks:
                try:
                    temp_rank = temp_edge.loc[edge, rank]
                    
                    ## Limit lowest taxonomic name to 'taxon'
                    
                    if rank == 'tax_name':
                        rank = 'taxon'
                        
                    taxon_map.loc[edge, rank] = temp_rank
                except KeyError:
                    continue
        
        if domain != 'eukarya':        
            for param in ['n16S', 'nge', 'ncds', 'genome_size', 'GC', 'phi', 'confidence']:
                edge_data.loc[name, param + '.mean'], edge_data.loc[name, param + '.sd'] = fill_edge_data(param, name, temp_edge)
                
        elif domain == 'eukarya':            
            print(f, 'please note that summary edge_data files are not created for domain eukarya')
            
        if len(pd.isnull(pd.DataFrame(temp_edge['nedge_corrected'])) > 0):
            temp_edge_abund = pd.DataFrame(temp_edge['nedge_corrected'])
        else:
            temp_edge_abund = pd.DataFrame(temp_edge['nedge'])
            
        temp_edge_abund.columns = [name]        
        edge_tally = pd.concat([edge_tally, temp_edge_abund], axis = 1)
        
if domain != 'eukarya':
    for f in os.listdir(cwd):
            
        if f.endswith(path_suffix):
            name = re.sub(path_suffix, '', f)
            temp_path = pd.read_csv(cwd + f, index_col = 0, names = [name])
            path_tally = pd.concat([path_tally, temp_path], axis = 1)
            
        elif f.endswith(ec_suffix):
            name = re.sub(ec_suffix, '', f)
            temp_ec = pd.read_csv(cwd + f, index_col = 0, names = [name])
            ec_tally = pd.concat([ec_tally, temp_ec], axis = 1)
        
## Write out datafiles
            
pd.DataFrame.to_csv(edge_tally.transpose(), prefix + '.' + domain + '.edge_tally.csv') 
pd.DataFrame.to_csv(taxon_map, prefix + '.' + domain + '.taxon_map.csv') 

if domain != 'eukarya':
    pd.DataFrame.to_csv(edge_data, prefix + '.' + domain + '.edge_data.csv')
    pd.DataFrame.to_csv(path_tally.transpose(), prefix + '.' + domain + '.path_tally.csv') 
    pd.DataFrame.to_csv(ec_tally.transpose(), prefix + '.' + domain + '.ec_tally.csv')
    
## Make combined unique file

unique_tally = pd.DataFrame()
unique_edge_num = {}

for f in os.listdir(cwd):
    if f.endswith(unique_suffix):
        name = re.sub(unique_suffix, '', f)
        temp_unique = pd.read_csv(cwd + f, index_col = None)
        temp_unique['seq'] = temp_unique['identifier'].str.split('|', expand = True).iloc[:,0]
        temp_unique.set_index('seq', inplace = True)
        
        for seq in temp_unique.index:
            try:
                temp_edge = unique_edge_num[seq]
                if temp_unique.loc[seq, 'global_edge_num'] not in temp_edge:
                    temp_edge.append(temp_unique.loc[seq, 'global_edge_num'])
            except KeyError:
                unique_edge_num[seq] = [temp_unique.loc[seq, 'global_edge_num']]
        
        unique_tally = pd.concat([unique_tally, temp_unique.abundance_corrected], axis = 1, sort = True)
        unique_tally.rename({'abundance_corrected':name})
    
pd.DataFrame.to_csv(unique_tally.transpose(), prefix + '.' + domain + '.unique_tally.csv')

## Write out a file mapping unique sequences to edges.  This is useful for identifying those sequences that
## placed to different edges in different samples.  These sequences should be manually curated for correct
## taxonomy.  Those sequences that placed to multiple edges will be listed first to make them easy to find.

with open(prefix + '.' + domain + '.seq_edge_map.csv', 'w') as seq_edge_out:
    for key in unique_edge_num.keys():
        if len(unique_edge_num[key]) > 1:
            temp_str = key + ',' + str(unique_edge_num[key]).strip('[]')
            print(temp_str, file = seq_edge_out)
    for key in unique_edge_num.keys():
        if len(unique_edge_num[key]) == 1:
            temp_str = key + ',' + str(unique_edge_num[key]).strip('[]')
            print(temp_str, file = seq_edge_out)                
        