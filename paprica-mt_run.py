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
    query = ['test.fastq.gz']
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
    threads = 4
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
    
#%% Parse SAM format output
    
## First need to determine the maximum number of fields.  This is variable because
## of optional field codes.  It may always be 17, though not clear why pandas doesn't
## recognize this on its own.
    
line_length = 12
i = 0
    
with open(cwd + name + '.sam') as sam:
    for line in sam:
        i = i + 1
        if line.startswith('@') == False:
            line = line.split('\t')
            
            if len(line) > line_length:
                line_length = len(line)
                
            print 'determing the number of fields:', i, line_length
            
## Then generate column names.

sam_columns = ['opt'] * line_length
sam_columns = [m + str(n) for m,n in zip(sam_columns, range(1, len(sam_columns)))]    
sam_columns[0:10] = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'mrnm', 'mpos', 'isize', 'seq', 'qual']

## Read the SAM file as a pandas df.

print 'reading ' + cwd + name + '.sam...'
sam = pd.read_csv(cwd + name + '.sam', header = None, comment = '@', sep = '\t', engine = 'python', names = sam_columns)

## Count the number of hits to each reference protein.

print 'tallying mapped reads...'
prot_counts = sam.rname.value_counts()
prot_counts.name = 'n_hits'

## Add information from paprica-mt.ec.csv.

prot_unique_cds_df = pd.read_csv(paprica_path + ref_dir + '/paprica-mt.ec.csv', header = 0, index_col = 0)
prot_unique_cds_df = pd.concat([prot_unique_cds_df, prot_counts], axis = 1, join_axes = [prot_unique_cds_df.index])
prot_unique_cds_df.dropna(subset = ['n_hits'], inplace = True)

## Write out the final csv file.

print 'writing output csv:', cwd + name + '.tally_ec.csv...'
columns_out = ['genome', 'domain', 'EC_number', 'product', 'n_hits']
prot_unique_cds_df.to_csv(cwd + name + '.tally_ec.csv', columns = columns_out)

print 'done'
