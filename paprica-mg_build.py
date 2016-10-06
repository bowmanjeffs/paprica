#!/usr/bin/env python
# -*- coding: utf-8 -*-

help_string = """
Created on Sun Jan 31 13:13:09 2016

@author: Jeff Bowman, bowmanjs@ldeo.columbia.edu

paprica is licensed under a Creative Commons Attribution-NonCommercial
4.0 International License.  IF you use any portion of paprica in your
work please cite:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described
by Metabolic Structure: A General Framework and Application to a Seasonally
Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic
Peninsula." PloS one 10.8 (2015): e0135868.

If your analysis makes specific use of pplacer, Infernal, DIAMOND, or pathway-tools
please make sure that you also cite the relevant publications.

This script identifies all the features in the Genbank files of completed
genomes that have an EC number and are thus useful for pathway prediction.  It
creates a nonredundant fasta of these sequences and a database of all the
feature qualifiers that are necessary to build Genbank files containing the
features "on the fly".

This script uses DIAMOND to create a database of the nonredundant fasta against
which query shotgun MG sequences reads can be searched.

CALL AS:
    python paprica-mg_build.py [options]
    
OPTIONS:
-ref_dir: The name of the directory containing the paprica database.  Not necessary
if your database is named ref_genome_database (the default).

"""

executable = '/bin/bash' # shell for executing commands

import os
import subprocess
import sys

from Bio import SeqIO

import pandas as pd
import numpy as np

## Define a stop function for troubleshooting.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]

## Read in profile.  Required variables are ref_dir. ###

## Read in command line arguments.

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

try:
    ref_dir = command_args['ref_dir']
except KeyError:
    ref_dir = 'ref_genome_database'
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'

paprica_path = os.path.dirname(os.path.abspath(__file__)) + '/' # The location of the actual paprica scripts.  
ref_dir = paprica_path + ref_dir

## Read in genome_data so that you can iterate by genomes that are actually
## used by paprica.

genome_data_bacteria = pd.read_csv(ref_dir + 'bacteria/genome_data.final.csv', index_col = 0, header = 0)
genome_data_bacteria = genome_data_bacteria.dropna(subset = ['clade'])
genome_data_bacteria['domain'] = 'bacteria'

genome_data_archaea = pd.read_csv(ref_dir + 'archaea/genome_data.final.csv', index_col = 0, header = 0)
genome_data_archaea = genome_data_archaea.dropna(subset = ['clade'])
genome_data_archaea['domain'] = 'archaea'

genome_data = pd.concat([genome_data_bacteria, genome_data_archaea])

## Iterate through all the files in refseq and find the gbk files.  First pass
## just counts the number of features with EC_number so that we can create a
## Numpy array the right size.

eci = 0

for d in genome_data.index:
    domain = genome_data.loc[d, 'domain']
    for f in os.listdir(ref_dir + domain + '/refseq/' + d):
        if f.endswith('gbff'):
            for record in SeqIO.parse(ref_dir + domain + '/refseq/' + d + '/' + f, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        if 'EC_number' in feature.qualifiers.keys():
                            ec = feature.qualifiers['EC_number']
                            for ec_number in ec:                                
                                eci = eci + 1
                                print 'counting features...', eci
                            
## Create numpy array for data and a 1D array that will become dataframe index.
                            
prot_array = np.empty((eci,7), dtype = 'object')
prot_array_index = np.empty(eci, dtype = 'object')

## Iterate through all the files in refseq and find the gbk files again.  Store
## the information necessary to create a Genbank record of each feature in the
## array.

i = 0

for d in genome_data.index:
    domain = genome_data.loc[d, 'domain']
    for f in os.listdir(ref_dir + domain + '/refseq/' + d):
        if f.endswith('gbff'):
            for record in SeqIO.parse(ref_dir + domain + '/refseq/' + d + '/' + f, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        if 'EC_number' in feature.qualifiers.keys():
                            
                            protein_id = feature.qualifiers['protein_id'][0]
                            trans = feature.qualifiers['translation'][0]
                            ec = feature.qualifiers['EC_number']
                            prod = feature.qualifiers['product'][0]
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            
                            for ec_number in ec:
                            
                                prot_array_index[i] = protein_id
                                prot_array[i,0] = d
                                prot_array[i,1] = domain
                                prot_array[i,2] = ec_number
                                prot_array[i,3] = trans
                                prot_array[i,4] = prod
                                prot_array[i,5] = start
                                prot_array[i,6] = end
                                
                                i = i + 1
                                print d, i, 'out of', eci, protein_id

## Convert array to pandas dataframe

columns = ['genome', 'domain', 'EC_number', 'translation', 'product', 'start', 'end']                                    
prot_df = pd.DataFrame(prot_array, index = prot_array_index, columns = columns)

## Determine how often each translation appears.  CDS with a translation
## that is not unique should not be used for taxonomic profiling.  Currently
## not using taxonomic profiling anyway.

prot_counts = pd.DataFrame(prot_df['translation'].value_counts())
prot_counts.columns = ['n_occurrences'] # The number of times that the sequence appears across all genomes.

## Add this information to prot_df.

prot_df = pd.merge(prot_df, prot_counts, left_on = 'translation', right_index = True)

## Now that you know which are duplicates, make nonredundant and print out to
## csv.

prot_unique_df = prot_df.drop_duplicates(subset = ['translation'])
prot_unique_df.to_csv(ref_dir + 'paprica-mg.ec.csv')

## Make a nonredundant fasta.

with open(ref_dir + 'paprica-mg.fasta', 'w') as fasta_out:
    for protein_id in prot_unique_df.index:
        print 'making fasta file', protein_id
        print >> fasta_out, '>' + protein_id
        print >> fasta_out, prot_unique_df.loc[protein_id, 'translation']

makedb = subprocess.Popen('diamond makedb --in ' + ref_dir + 'paprica-mg.fasta -d ' + ref_dir + 'paprica-mg', shell = True, executable = executable)
makedb.communicate()    
                        
                        
        