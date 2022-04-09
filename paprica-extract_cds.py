#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Sun March 20 12:36:39 2022

@author: Jeff Bowman, jsbowman@ucsd.edu

paprica is licensed under a Creative Commons Attribution-NonCommercial
4.0 International License.  IF you use any portion of paprica in your
work please cite:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described
by Metabolic Structure: A General Framework and Application to a Seasonally
Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic
Peninsula." PloS one 10.8 (2015): e0135868.

Please all cite core dependencies as indicated here:
    https://github.com/bowmanjeffs/paprica#citations.

REQUIRES:
    Files:
        [domain]/refseq/[assembly]/[assembly_identifier.gbff]

    Programs:
        none

    Python modules:
        Bio
        os
        sys
        
RUN AS:
    paprica-extract_cds.py [options]
    
OPTIONS:
    -domain: Which domain are you analyzing?  Either bacteria or archaea.
    -ref_dir: The name for the database you are building.  The default is "ref_genome_database".

This script must be located in the 'paprica' directory as it makes use of relative
paths.

This script loops across all genome directories in refseq for a given domain
and extracts all CDS to a new fasta file. The CDS fasta files are used by
paprica-run_gRodon.r to predict growth rates.
    
"""


from Bio import SeqIO
import os
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
    print(help_string)
    quit()
    
## set defaults
    
try:        
    domain = command_args['domain']
except KeyError:
    domain = 'archaea'
try:
    ref_dir = command_args['ref_dir']
except KeyError:
    ref_dir = 'ref_genome_database'

ref_dir_domain = ref_dir + '/' + domain + '/'

## get genomes that have already been done as set

try:
    done = set(os.listdir(ref_dir_domain + 'cds'))
except FileNotFoundError:
    os.mkdir(ref_dir_domain + 'cds')
    done = set()
    
## extract CDS
    
#!!! This should be functionalized and implemented in parallel.  You should
## also loop across genome_data.csv and not /refseq as that contains very many
## genomes that aren't used.

with open(domain + '_extract_cds_error_file.txt', 'w') as error_file:
    for d in os.listdir(ref_dir_domain + 'refseq'):
        with open(ref_dir_domain + 'cds/' + d + '.cds.fasta', 'w') as file_out:
            for f in os.listdir(ref_dir_domain + 'refseq/' + d):
                if f.endswith('.gbff'):
                    
                    basename = f.split('.gbff')
                    basename = basename[0:-1]
                    basename = '.'.join(basename)
                    
                    if basename + '.cds.fasta' not in done:
                        with open(ref_dir_domain + 'refseq/' + d + '/' + f, 'r') as file_in:
                            i = 0
                            try:
                                for record in SeqIO.parse(file_in, 'genbank'):
                                    for feature in record.features:
                                        if feature.type=="CDS":
                                            i += 1
                                            try:
                                                product = str(i) + '_' + feature.qualifiers['product'][0]
                                            except KeyError:
                                                product = str(i) + '_noproduct'
                                            seq = feature.extract(record.seq)
                                            print('>' + product + '\n' + seq, file = file_out)
                            except ValueError:
                                print(f, file = error_file)
        print('Extracting CDS for', d)

             

                    