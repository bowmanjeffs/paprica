#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Sun Oct 11 21:20:57 2015
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
    Python modules:
        pandas
        numpy
        Bio
        
CALL AS:
    python paprica_tally_pathways.py [options]
    
OPTIONS:
    -cutoff: The fraction of terminal daughters that need to have a pathway for it
        to be included in an internal node, between 0-1
    -domain: domain of analysis, either bacteria or archaea
    -i: input csv
    -o: prefix for output files
    -ref_dir: name of reference directory
    -override ["old|new,old|new"]: any known incorrect/correct edge pair
        replacements to be made, note that quotes are necessary.
        Example: -override ["527|13,1140:2139"] will substitute edge 13 for
        527, and 2139 for 1140.
    -omit [start:stop]: a range of edges (e.g., cyanobacteria if you need to 
        eliminate possible chloroplasts) that should be omitted.
        Example: -omit 5:9 to omit edges 5:9
        
This script must be located in the 'paprica' directory as it makes use of relative
paths.

Although you do not need to specify it as an input file, this script also requires
the stockholm-format combined reference and query alignment produced by paprica-place_it.py.

"""

import pandas as pd
import numpy as np

import sys
import os

try:
    paprica_path = os.path.dirname(os.path.realpath(__file__)) + '/' # The location of the actual paprica scripts.
except NameError:
    paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/'

cwd = os.getcwd() + '/'  # The current working directory
    
## Parse command line arguments.
                
command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in list(command_args.keys()):
    print(help_string)
    quit()
        
## If any command line options are specified all need to be specified except
## overrides and omit.

if len(sys.argv) > 2:               
    cutoff = float(command_args['cutoff'])  # The cutoff value used to determine pathways to include for internal nodes.
    domain = command_args['domain']  # The domain (bacteria or archaea) for analysis.
    ref_dir = paprica_path + command_args['ref_dir']  # The complete path to the reference directory being used for analysis.        
    name = command_args['o']
    query = command_args['i']
      
    try:
        overrides = command_args['override']
    except KeyError:
        overrides = ''
    try:
        omit = command_args['omit']
    except KeyError:
        omit = ''
    
else:
    query = 'test.bacteria.combined_16S.bacteria.tax.placements.csv'
    name = 'test_bacteria'
    cutoff = 0.5  # The cutoff value used to determine pathways to include for internal nodes.
    domain = 'bacteria'  # The domain (bacteria or archaea) for analysis.
    ref_dir = paprica_path + 'ref_genome_database'  # The complete path to the reference directory being used for analysis.        
    #omit = '674:818'
    #overrides = '5804|93,4619|4571'
    overrides = ''
    omit = ''
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
ref_dir_domain = ref_dir + domain + '/'  # The complete path the the domain subdirectory of the reference directory.

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])
    
## Create a dictionary of any edges that need replacement.

override_dic = {}

if len(overrides) > 0:
    overrides = overrides.split(',')
    for pair in overrides:
        pair = pair.split('|')
        override_dic[int(pair[0])] = int(pair[1])
        
## Read in the query csv file generated by paprica-place_it.

query_csv = pd.read_csv(cwd + query, header = 0, index_col = 'Pquery')

## Import csv files generated by paprica_build_core_genomes.

genome_data = pd.read_csv(ref_dir_domain + 'genome_data.final.csv.gz', header = 0, index_col = 0)
lineages = pd.read_csv(ref_dir_domain + 'node_lineages.csv.gz', header = 0, index_col = 0)

if domain != 'eukarya':
    internal_data = pd.read_csv(ref_dir_domain + 'internal_data.csv.gz', header = 0, index_col = 0)
    
    terminal_paths = pd.read_csv(ref_dir_domain + 'terminal_paths.csv.gz', header = 0, index_col = 0)
    internal_probs = pd.read_csv(ref_dir_domain + 'internal_probs.csv.gz', header = 0, index_col = 0)
    
    internal_ec_probs = pd.read_csv(ref_dir_domain + 'internal_ec_probs.csv.gz', header = 0, index_col = 0)
    internal_ec_n = pd.read_csv(ref_dir_domain + 'internal_ec_n.csv.gz', header = 0, index_col = 0)
    terminal_ec = pd.read_csv(ref_dir_domain + 'terminal_ec.csv.gz', header = 0, index_col = 0)
    
    internal_probs = internal_probs.fillna(0)
    internal_ec_probs = internal_ec_probs.fillna(0)
    internal_ec_n = internal_ec_n.fillna(0)
    terminal_ec = terminal_ec.fillna(0)

## Inflate the dataframe based on the values in the abundance column.

unique_csv = query_csv.copy()
query_csv = query_csv.loc[query_csv.index.repeat(query_csv.abundance)]
query_csv.drop(columns = 'abundance', inplace = True)

## Override bad edges.

for edge in list(override_dic.keys()):
    query_csv.loc[query_csv['global_edge_num'] == edge, 'global_edge_num'] = override_dic[edge]
    query_csv.loc[query_csv['global_edge_num'] == edge, 'like_weight_ratio'] = np.NaN
    query_csv.loc[query_csv['global_edge_num'] == edge, 'likelihood'] = np.NaN
    query_csv.loc[query_csv['global_edge_num'] == edge, 'distal_length'] = np.NaN
    query_csv.loc[query_csv['global_edge_num'] == edge, 'pendant_length'] = np.NaN

## Tally the number of occurences of each edge in the sample and
## get the mean posterior probability, overlap, and map ratio for each edge.

## If all placements are to internal edges there will be no map_id and map_ratio.
    
try:
    edge_map_id = query_csv.groupby('global_edge_num').map_id.mean()
    edge_map_ratio = query_csv.groupby('global_edge_num').map_ratio.mean()
except AttributeError:
    edge_map_id = np.nan
    edge_map_ratio = np.nan

edge_tally = query_csv.groupby('global_edge_num').size()
edge_lwr = query_csv.groupby('global_edge_num').like_weight_ratio.mean()
edge_l = query_csv.groupby('global_edge_num').likelihood.mean()
edge_edpl = query_csv.groupby('global_edge_num').EDPL.mean()

## Omit undesired edges.

if len(omit) > 0:
    omit = omit.split(':')
    drop_edges = list(range(int(omit[0]), int(omit[1]) + 1))
    edge_tally = edge_tally.drop(drop_edges, errors = 'ignore')

## Add the edge tally and mean pp to a new data frame that will hold other
## sample information.

edge_data = pd.DataFrame(index = edge_tally.index)
edge_data['nedge'] = edge_tally
edge_data['like_weight_ratio'] = edge_lwr
edge_data['likelihood'] = edge_l
edge_data['map_ratio'] = edge_map_ratio
edge_data['map_id'] = edge_map_id
edge_data['EDPL'] = edge_edpl

## This makes sure that there is an nedge_corrected entry for eukarya.
## Will be overwritten for bacteria and archaea.

edge_data['nedge_corrected'] = edge_data['nedge']

if domain == 'eukarya':
    edge_data['taxon'] = lineages.reindex(edge_data.index).consensus

## Dataframe to hold the number of occurences of pathway in sample, by edge.

if domain != 'eukarya':
    sample_pathways = pd.DataFrame(index = sorted(terminal_paths.columns))
    sample_ec = pd.DataFrame(index = sorted(terminal_ec.columns))

    for edge in list(edge_tally.index):
        print('generating data for edge', edge)
        
        ## If edge is an internal node...
        
        if edge in internal_probs.index:
                    
            ## Collect other information that you might want later. 
            
            edge_data.loc[edge, 'genome_size'] = internal_data.loc[edge, 'genome_size']
            edge_data.loc[edge, 'clade_size'] = internal_data.loc[edge, 'clade_size']
            edge_data.loc[edge, 'npaths_terminal'] = internal_data.loc[edge, 'npaths_terminal']
            edge_data.loc[edge, 'nec_terminal'] = internal_data.loc[edge, 'nec_terminal']
            edge_data.loc[edge, 'branch_length'] = internal_data.loc[edge, 'branch_length']
            edge_data.loc[edge, 'taxon'] = lineages.loc[edge, 'consensus']
            
            edge_data.loc[edge, 'n16S'] = internal_data.loc[edge, 'n16S']
            edge_data.loc[edge, 'GC'] = internal_data.loc[edge, 'GC']
            edge_data.loc[edge, 'phi'] = internal_data.loc[edge, 'phi']
            edge_data.loc[edge, 'ncds'] = internal_data.loc[edge, 'ncds']
            edge_data.loc[edge, 'nge'] = internal_data.loc[edge, 'nge']
            edge_data.loc[edge, 'nedge_corrected'] = float(edge_data.loc[edge, 'nedge']) / float(internal_data.loc[edge, 'n16S'])
        
            ## Get the pathways associated with the edge.  Report the abundance of
            ## pathways as the 16S copy number corrected abundance of edge.
            
            edge_pathways = internal_probs.loc[edge, internal_probs.loc[edge, :] >= cutoff]
            edge_pathways.loc[:] = edge_data.loc[edge, 'nedge_corrected']
            sample_pathways.loc[:, edge] = edge_pathways
            edge_data.loc[edge, 'npaths_actual'] = edge_pathways.count() # How many pathways are present in terminal daughters above the cutoff?
            
            ## Get the enzymes associated with the edge.  For this to work the
            ## columns for internal_ec_n and internal_ec_probs MUST be in the
            ## same order.
            
            edge_ec_n = internal_ec_n.loc[edge, internal_ec_probs.loc[edge, :] >= cutoff]
            edge_data.loc[edge, 'nec_actual'] = edge_ec_n.sum()
            edge_ec_n = edge_ec_n.mul(edge_data.loc[edge, 'nedge_corrected'])
            sample_ec.loc[:, edge] = edge_ec_n
    
            ## Calculate the confidence score.  This differs from PAPRICA_v0.11 in that the number
            ## of pathways in the edge relative to the terminal clade members is used in place of
            ## the number of CDS.
    
            npaths_actual = edge_data.loc[edge, 'npaths_actual']
            npaths_terminal = edge_data.loc[edge, 'npaths_terminal']
        
            phi = edge_data.loc[edge, 'phi']
            confidence = (npaths_actual / npaths_terminal) * (1 - phi)
            edge_data.loc[edge, 'confidence'] = confidence 
    
        ## If edge is a terminal node...
            
        else:
            
            ## Now get some useful data for the edge.  If the edge belongs to
            ## phylum_refs, meaning that it could not be placed to any of the
            ## nodes in a subtree, it will not have an entry except for clade_size.
            
            edge_data.loc[edge, 'clade_size'] = 1
            
            try:
            
                edge_data.loc[edge, 'taxon'] = lineages.loc[edge, 'consensus']
                edge_data.loc[edge, 'branch_length'] = genome_data.loc[genome_data['clade'] == edge, 'branch_length'][0]
                edge_data.loc[edge, 'nedge_corrected'] = edge_data.loc[edge, 'nedge']
                edge_data.loc[edge, 'GC'] = genome_data.loc[genome_data['clade'] == edge, 'GC'][0]
                edge_data.loc[edge, 'phi'] = genome_data.loc[genome_data['clade'] == edge, 'phi'][0]
                edge_data.loc[edge, 'genome_size'] = genome_data.loc[genome_data['clade'] == edge, 'genome_size'][0]
                edge_data.loc[edge, 'ncds'] = genome_data.loc[genome_data['clade'] == edge, 'ncds'][0]
                edge_data.loc[edge, 'nge'] = genome_data.loc[genome_data['clade'] == edge, 'nge'][0]
                edge_data.loc[edge, 'n16S'] = genome_data.loc[genome_data['clade'] == edge, 'n16S'][0]
                edge_data.loc[edge, 'confidence'] = genome_data.loc[genome_data['clade'] == edge, 'phi'][0] # Phi for terminal nodes
                edge_data.loc[edge, 'nedge_corrected'] = float(edge_data.loc[edge, 'nedge']) / float(genome_data.loc[genome_data['clade'] == edge, 'n16S'])

            except KeyError:
                continue
        
            ## Get the pathways associated with the edge.  The pathways are indexed by assembly not edge number.
            
            assembly = genome_data[genome_data['clade'] == edge].index.tolist()[0]
            edge_pathways = terminal_paths.loc[assembly, terminal_paths.loc[assembly, :] == 1]
            
            ## For bacteria and archaea, correct for multiple 16S rRNA gene copies.
            ## The assumption here is that each genome has only one copy of a pathway.
            
            edge_pathways.loc[:] = edge_data.loc[edge, 'nedge_corrected'] 
                    
            edge_data.loc[edge, 'npaths_terminal'] = np.nan
            edge_data.loc[edge, 'npaths_actual'] = genome_data.loc[genome_data['clade'] == edge, 'npaths_actual'][0]
    
            sample_pathways.loc[:, edge] = edge_pathways
            
            ## Get the EC numbers associated with the edge.
            
            edge_ec_n = terminal_ec.loc[assembly, terminal_ec.loc[assembly, :] >= 1]
            edge_data.loc[edge, 'nec_actual'] = edge_ec_n.sum()
            
            ## For bacteria and archaea, correct for multiple 16S rRNA gene copies.
            
            edge_ec_n = edge_ec_n.mul(edge_data.loc[edge, 'nedge_corrected'])
    
            edge_data.loc[edge, 'nec_actual'] = edge_ec_n.sum()
            edge_data.loc[edge, 'nec_terminal'] = np.nan
    
            sample_ec.loc[:, edge] = edge_ec_n      
        
    # if 'taxon' not in edge_data.columns:
        
    #     ## This means that none of the edges had valid taxon data, something
    #     ## to address later. Quick fix here to make sure that combine_edge_results.py
    #     ## doesn't fail.
        
    #     edge_data['taxon'] = ''
    
    ## Calculate the confidence score for the sample. If your input file contains
    ## only a single read and it does not place to a subtree there will be no
    ## confidence value.
    
    try:
        sample_confidence = sum((edge_data['confidence'] * edge_data['nedge_corrected'])) / edge_data['nedge_corrected'].sum() 
    except KeyError:
        sample_confidence = np.nan
    
## Add lineage data for each edge

edge_data = pd.concat([edge_data, lineages], axis = 1, join = 'inner')
    
#%% Prepare unique read file, annotating with corrected read number and taxonomy.

print('Normalizing abundance for unique sequences...')

## Try clause is necessary for cases where there are only internal placements,
## in which case there will be no map_ratio or map_id in the placements.csv
## file.

try:
    unique_csv = unique_csv[['seq', 'abundance', 'global_edge_num', 'map_ratio', 'map_id', 'EDPL']]
except KeyError:
    unique_csv = unique_csv[['seq', 'abundance', 'global_edge_num', 'EDPL']]
    unique_csv['map_ratio'] = np.nan
    unique_csv['map_id'] = np.nan
    unique_csv = unique_csv[['seq', 'abundance', 'global_edge_num', 'map_ratio', 'map_id', 'EDPL']]

unique_csv.loc[unique_csv.index, 'name'] = unique_csv.index
unique_csv.index = unique_csv.seq

## then you need abundance_corrected and taxon from edge_data

for unique in unique_csv.index:
    edge = unique_csv.loc[unique, 'global_edge_num']
    
    if domain != 'eukarya':
        try:
            n16S = edge_data.loc[edge, 'n16S']
        except KeyError:
            n16S = 1
    else:
        n16S = 1
    
    unique_csv.loc[unique, 'abundance_corrected'] = unique_csv.loc[unique, 'abundance'] / float(n16S)  
    unique_csv.loc[unique, 'identifier'] = str(unique) + '|' + edge
        
    try:
        unique_csv.loc[unique, 'taxon'] = edge_data.loc[edge, 'taxon']  
    except KeyError:
        pass

unique_csv.index = unique_csv.identifier
unique_csv.drop('identifier', axis = 1, inplace = True)
unique_csv.drop('seq', axis = 1, inplace = True)    
unique_csv.to_csv(cwd + name + '.unique_seqs.csv')

edge_data.to_csv(cwd + name + '.edge_data.csv')

#%%

npathways = 'NA'
ppathways = 'NA'
nreads = edge_data['nedge'].sum()

## Generate a single column table of the total (corrected) abundance for each
## pathway and enzyme.  Absent pathways/enzymes are included as 0 to make it
## easier to compare between samples.

if domain != 'eukarya':
    sample_pathways_sum = sample_pathways.sum(1)
    npathways = len(sample_pathways_sum[sample_pathways_sum != 0])
    ppathways = len(sample_pathways_sum)

    ## Generate a single column table of the total (corrected) abundance for each
    ## enzyme.  Absent enzymes are included as 0, to make it easier to compare
    ## between samples.
    
    sample_ec_sum = sample_ec.sum(1)
    nec = len(sample_ec_sum[sample_ec_sum != 0])
    pec = len(sample_ec_sum)

    ## Write out all the tables.
    
    sample_pathways = sample_pathways.fillna(0)
    sample_ec = sample_ec.fillna(0)

    sample_pathways_sum.to_csv(cwd + name + '.sum_pathways.csv', header = False)
    sample_pathways.to_csv(cwd + name + '.pathways.csv')
    
    sample_ec_sum.to_csv(cwd + name + '.sum_ec.csv', header = False)
    sample_ec.to_csv(cwd + name + '.ec.csv')

## Get the database creation time, this serves as a version.

for f in os.listdir(os.path.expanduser(ref_dir_domain)):
    if f.endswith('.database_info.txt'):
        with open(os.path.expanduser(ref_dir_domain) + f, 'r') as database_info:
            for line in database_info:
                if 'ref tree built at:' in line:
                    line = line.rstrip()
                    line = line.split(': ')
                    database_time = line[1]
                    database_time = database_time.strip()

## And a simple tab-delim for the sample data.

with open(cwd + name + '.sample_data.txt', 'w') as sample_data:
    print('name' + '\t' + name, file=sample_data)
    
    if domain != 'eukarya':
        print('sample_confidence' + '\t' + str(sample_confidence), file=sample_data)

    print('npathways' + '\t' + str(npathways), file=sample_data)
    print('ppathways' + '\t' + str(ppathways), file=sample_data)
    print('nreads' + '\t' + str(nreads), file=sample_data)
    print('database_created_at' + '\t' + database_time, file=sample_data)

