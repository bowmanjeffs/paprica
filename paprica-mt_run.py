#!/usr/bin/env python
# -*- coding: utf-8 -*-

help_string = """
Created on Mon Feb 01 16:50:01 2016

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

CALL AS:
    python paprica-mt_run.py [OPTIONS]
    
OPTIONS:
    -i: Input fasta or fasta.gz.  If you have PE reads you would enter two
        files, like paprica-mt_run.py -i file1.fasta.gz file2.fasta.gz...
    -o: Prefix for output files.
    -ref_dir: Name of the directory containing the paprica database.
    -pgdb_dir: The location where pathway-tools stores PGDBs.  Only necessary
        if pathways are predicted.
    -pathways: T or F, whether pathways should be predicted.  This can take a
        long time.
    -t: The number of threads for bwa to use.

"""

executable = '/bin/bash' # shell for executing commands

import sys
import subprocess
import os
import gzip

import pandas as pd

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]

## Identify directory locations.
    
cwd = os.getcwd() + '/' # The current working directory.

if len(sys.argv) == 1:
    paprica_path = '/volumes/hd1/paprica/' # The location of the paprica scripts during testing.
else:
    paprica_path = os.path.dirname(os.path.abspath(__file__)) + '/' # The location of the actual paprica scripts.
    
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

## Provide input switches for testing.

if 'i' not in command_args.keys():
    query = ['test.fasta.gz']
else:
    query = command_args['i']
    query = query.split()
    
if 'o' not in command_args.keys():
    name = 'test_mt'
else:
    name = command_args['o']
    
if 'ref_dir' not in command_args.keys():
    ref_dir = 'ref_genome_database'
else:
    ref_dir = command_args['ref_dir']
    
if 'pathways' not in command_args.keys():
    pathways = 'F'
else:
    pathways = command_args['pathways']
    
if 'pgdb_dir' not in command_args.keys():
    pgdb_dir = '/volumes/hd2/ptools-local/user/'
else:
    pgdb_dir = command_args['pgdb_dir']
    
if 't' not in command_args.keys():
    threads = 72
else:
    threads = command_args['t']
    
## Define path to reference directory.
    
ref_dir_path = paprica_path + ref_dir

if ref_dir_path.endswith('/') == False:
    ref_dir_path = ref_dir_path + '/'
    
#%% Use BWA to search the query against the paprica-mt database

if len(query) == 1:
    bwa_aln = subprocess.Popen('bwa mem ' \
    + '-t ' + str(threads) + ' ' \
    + paprica_path + ref_dir + '/paprica-mt.fasta ' \
    + cwd + query[0] + ' > ' \
    + cwd + name + '.sam', shell = True, executable = executable)
    bwa_aln.communicate() 

## Option for PE reads.

if len(query) > 1:
    
    bwa_aln = subprocess.Popen('bwa mem ' \
    + '-t ' + str(threads) + ' ' \
    + paprica_path + ref_dir + '/paprica-mt.fasta ' \
    + cwd + query[0] + \
    + ' ' + cwd + query[1] + ' > ' \
    + cwd + name + '.sam', shell = True, executable = executable)
    bwa_aln.communicate() 
    
## gzip the sam file to save space.
    
gz = subprocess.Popen('gzip -f ' + cwd + name + '.sam', shell = True, executable = executable)
gz.communicate()
    
#%% Iterate across sam file, tallying the number of hits to each reference that appears in the results.
    
prot_counts = pd.Series()
prot_counts.name = 'n_hits'
    
i = 0
f = 0
    
with gzip.open(cwd + name + '.sam.gz', 'rb') as sam:
    for line in sam:
        
        if line.startswith('@') == False:
            i = i + 1
            line = line.split('\t')
            
            ## For mapped reads, add to tally for the reference sequence.
            
            if line[1] in ['0', '2', '16']:            
                rname = line[2]
                f = f + 1
                
                try:
                    prot_counts[rname] = prot_counts[rname] + 1
                except KeyError:
                    prot_counts[rname] = 1
    
                print 'tallying hits for', name + ':', 'found', f, 'out of', i

## Add information from paprica-mt.ec.csv.

prot_unique_cds_df = pd.read_csv(paprica_path + ref_dir + '/paprica-mt.ec.csv', header = 0, index_col = 0)
prot_unique_cds_df = pd.concat([prot_unique_cds_df, prot_counts], axis = 1, join_axes = [prot_unique_cds_df.index])
prot_unique_cds_df.dropna(subset = ['n_hits'], inplace = True)
prot_unique_cds_df['length_cds'] = prot_unique_cds_df.translation.str.len() # CDS length, would be nice if this was precomputed

## Write out the final csv file.

print 'writing output csv:', cwd + name + '.tally_ec.csv...'
columns_out = ['genome', 'domain', 'EC_number', 'product', 'length_cds', 'n_hits']
prot_unique_cds_df.to_csv(cwd + name + '.tally_ec.csv', columns = columns_out)

## Write out report file.

with open(cwd + name + '.paprica-mt_report.txt', 'w') as report:
    print >> report, 'file', query
    print >> report, 'n_reads', i
    print >> report, 'n_hits', f
    print >> report, 'f_hits', float(f)/i
    print >> report, 'n_genomes', len(prot_unique_cds_df.genome.value_counts())

print 'done'
