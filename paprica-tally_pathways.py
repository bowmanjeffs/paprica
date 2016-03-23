#!/usr/bin/env python
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
        
CALL AS:
    python paprica_tally_pathways.py [options]
    
OPTIONS:
    -cutoff: The fraction of terminal daughters that need to have a pathway for it
        to be included in an internal node, between 0-1
    -domain: domain of analysis, either bacteria or archaea
    -i: input csv
    -o: prefix for output files
    -ref_dir: name of reference directory

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""

import pandas as pd
import numpy as np

import sys
import os

paprica_path = os.path.dirname(os.path.abspath(__file__)) + '/' # The location of the actual paprica scripts.
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
        
if 'h' in command_args.keys():
    print help_string
    quit()
        
## If any command line options are specified all need to be specified.

if len(sys.argv) > 2:               
    cutoff = float(command_args['cutoff'])  # The cutoff value used to determine pathways to include for internal nodes.
    domain = command_args['domain']  # The domain (bacteria or archaea) for analysis.
    ref_dir = paprica_path + command_args['ref_dir']  # The complete path to the reference directory being used for analysis.        
    query = command_args['i']
    name = command_args['o']
    
else:
    query = 'test.archaea.combined_16S.archaea.tax.clean.align.csv'
    name = 'test.archaea'
    cutoff = 0.5  # The cutoff value used to determine pathways to include for internal nodes.
    domain = 'archaea'  # The domain (bacteria or archaea) for analysis.
    ref_dir = paprica_path + 'ref_genome_database'  # The complete path to the reference directory being used for analysis.        

## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
ref_dir_domain = ref_dir + domain + '/'  # The complete path the the domain subdirectory of the reference directory.

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]

## Import csv files generated by paprica_build_core_genomes.

genome_data = pd.DataFrame.from_csv(ref_dir_domain + 'genome_data.final.csv', header = 0, index_col = 0)
internal_data = pd.DataFrame.from_csv(ref_dir_domain + 'internal_data.csv', header = 0, index_col = 0)

terminal_paths = pd.DataFrame.from_csv(ref_dir_domain + 'terminal_paths.csv', header = 0, index_col = 0)
internal_probs = pd.DataFrame.from_csv(ref_dir_domain + 'internal_probs.csv', header = 0, index_col = 0)

internal_ec_probs = pd.DataFrame.from_csv(ref_dir_domain + 'internal_ec_probs.csv', header = 0, index_col = 0)
internal_ec_n = pd.DataFrame.from_csv(ref_dir_domain + 'internal_ec_n.csv', header = 0, index_col = 0)
terminal_ec = pd.DataFrame.from_csv(ref_dir_domain + 'terminal_ec.csv', header = 0, index_col = 0)

internal_probs = internal_probs.fillna(0)
internal_ec_probs = internal_ec_probs.fillna(0)
internal_ec_n = internal_ec_n.fillna(0)
terminal_ec = terminal_ec.fillna(0)

## Read in the query csv file generated by paprica_place_it.

query_csv = pd.DataFrame.from_csv(cwd + query, header = 0)

## Tally the number of occurences of each edge in the sample and
## get the mean posterior probability for each edge.

edge_tally = query_csv.groupby('edge_num').size()
edge_pp = edge_pp = query_csv.groupby('edge_num').post_prob.mean()

## Add the edge tally and mean pp to a new data frame that will hold other
##sample information.  Not necessary to define column names in advance.

edge_data = pd.DataFrame(index = edge_tally.index, columns = ['taxon', 'nedge', 'n16S', 'nedge_corrected', 'nge', 'ncds', 'genome_size', 'GC', 'phi', 'clade_size', 'branch_length', 'npaths_terminal', 'npaths_actual', 'confidence'])
edge_data['nedge'] = edge_tally
edge_data['post_prob'] = edge_pp

## Dataframe to hold the number of occurences of pathway in sample, by edge.

sample_pathways = pd.DataFrame(index = sorted(terminal_paths.columns))
sample_ec = pd.DataFrame(index = sorted(terminal_ec.columns))

for edge in list(edge_tally.index):
    print 'generating data for edge', edge
    
    ## If edge is an internal node...
    
    if edge in internal_probs.index:
                
        ## Collect other information that you might want later.
        
        edge_data.loc[edge, 'taxon'] = np.nan # Taxon is indeterminate for internal nodes.
        edge_data.loc[edge, 'n16S'] = internal_data.loc[edge, 'n16S']
        edge_data.loc[edge, 'nedge_corrected'] = float(edge_data.loc[edge, 'nedge']) / float(internal_data.loc[edge, 'n16S'])
        edge_data.loc[edge, 'nge'] = internal_data.loc[edge, 'nge']
        edge_data.loc[edge, 'ncds'] = internal_data.loc[edge, 'ncds']
        edge_data.loc[edge, 'genome_size'] = internal_data.loc[edge, 'genome_size']
        edge_data.loc[edge, 'phi'] = internal_data.loc[edge, 'phi']
        edge_data.loc[edge, 'clade_size'] = internal_data.loc[edge, 'clade_size']
        edge_data.loc[edge, 'npaths_terminal'] = internal_data.loc[edge, 'npaths_terminal']
        edge_data.loc[edge, 'nec_terminal'] = internal_data.loc[edge, 'nec_terminal']
        edge_data.loc[edge, 'branch_length'] = internal_data.loc[edge, 'branch_length']
        edge_data.loc[edge, 'GC'] = internal_data.loc[edge, 'GC']
        
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
        
        ## Now get some useful data for the edge.
        
        edge_data.loc[edge, 'taxon'] = genome_data.loc[genome_data['clade'] == edge, 'tax_name'][0]
        edge_data.loc[edge, 'n16S'] = genome_data.loc[genome_data['clade'] == edge, 'n16S'][0]
        edge_data.loc[edge, 'nedge_corrected'] = float(edge_data.loc[edge, 'nedge']) / float(genome_data.loc[genome_data['clade'] == edge, 'n16S'])
        edge_data.loc[edge, 'nge'] = genome_data.loc[genome_data['clade'] == edge, 'nge'][0]
        edge_data.loc[edge, 'ncds'] = genome_data.loc[genome_data['clade'] == edge, 'ncds'][0]
        edge_data.loc[edge, 'genome_size'] = genome_data.loc[genome_data['clade'] == edge, 'genome_size'][0]
        edge_data.loc[edge, 'phi'] = genome_data.loc[genome_data['clade'] == edge, 'phi'][0]
        edge_data.loc[edge, 'clade_size'] = 1
        edge_data.loc[edge, 'branch_length'] = genome_data.loc[genome_data['clade'] == edge, 'branch_length'][0]
        edge_data.loc[edge, 'GC'] = genome_data.loc[genome_data['clade'] == edge, 'GC'][0]
        
        ## Get the pathways associated with the edge.  The pathways are indexed by assembly not edge number.
        
        assembly = genome_data[genome_data['clade'] == edge].index.tolist()[0]
        edge_pathways = terminal_paths.loc[assembly, terminal_paths.loc[assembly, :] == 1]
        edge_pathways.loc[:] = edge_data.loc[edge, 'nedge_corrected']
        sample_pathways.loc[:, edge] = edge_pathways        
                
        edge_data.loc[edge, 'npaths_terminal'] = np.nan
        edge_data.loc[edge, 'npaths_actual'] = genome_data.loc[genome_data['clade'] == edge, 'npaths_actual'][0]
        edge_data.loc[edge, 'confidence'] = genome_data.loc[genome_data['clade'] == edge, 'phi'][0] # Phi for terminal nodes
        
        ## Get the EC numbers associated with the edge.
        
        edge_ec_n = terminal_ec.loc[assembly, terminal_ec.loc[assembly, :] >= 1]
        edge_data.loc[edge, 'nec_actual'] = edge_ec_n.sum()
        edge_ec_n = edge_ec_n.mul(edge_data.loc[edge, 'nedge_corrected'])
        sample_ec.loc[:, edge] = edge_ec_n        
        edge_data.loc[edge, 'nec_terminal'] = np.nan

## Calculate the confidence score for the sample.

sample_confidence = sum((edge_data['confidence'] * edge_data['nedge_corrected'])) / edge_data['nedge_corrected'].sum() 

## Generate a single column table of the total (corrected) abundance for each
## pathway.  Absent pathways are included as 0, to make it easier to compare
## between samples.

sample_pathways_sum = sample_pathways.sum(1)
npathways = len(sample_pathways_sum[sample_pathways_sum != 0])
ppathways = len(sample_pathways_sum)
nreads = edge_data['nedge'].sum()

## Generate a single column table of the total (corrected) abundance for each
## enzyme.  Absent enzymes are included as 0, to make it easier to compare
## between samples.

sample_ec_sum = sample_ec.sum(1)
nec = len(sample_ec_sum[sample_ec_sum != 0])
pec = len(sample_ec_sum)

## Write out all the tables.

sample_pathways = sample_pathways.fillna(0)
sample_ec = sample_ec.fillna(0)

edge_data.to_csv(cwd + name + '.edge_data.csv')

sample_pathways_sum.to_csv(cwd + name + '.sum_pathways.csv')
sample_pathways.to_csv(cwd + name + '.pathways.csv')

sample_ec_sum.to_csv(cwd + name + '.sum_ec.csv')
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
    print >> sample_data, 'name' + '\t' + name
    print >> sample_data, 'sample_confidence' + '\t' + str(sample_confidence)
    print >> sample_data, 'npathways' + '\t' + str(npathways)
    print >> sample_data, 'ppathways' + '\t' + str(ppathways)
    print >> sample_data, 'nreads' + '\t' + str(nreads)
    print >> sample_data, 'database_created_at' + '\t' + database_time