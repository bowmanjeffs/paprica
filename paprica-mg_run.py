# -*- coding: utf-8 -*-
help_string = """
Created on Mon Feb 01 16:50:01 2016

@author: jeff

CALL AS:
    python paprica-mg_run.py [OPTIONS]
    
OPTIONS:
-i: Input fasta or fasta.gz.
-o: Prefix for output files.
-ref_dir: Name of the directory containing the paprica database.
-pgdb_dir: The location where pathway-tools stores PGDBs.

"""

executable = '/bin/bash' # shell for executing commands

import sys
import subprocess
import os
import shutil

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqFeature

## Identify the location of this script.
    
paprica_path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'

## Identify the current working directory for output.
    
cwd = os.getcwd() + '/'

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

ref_dir = paprica_path + command_args['ref_dir']
pgdb_dir = command_args['pgdb_dir']

## Read in command line arguments.

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]

## Provide input switches for testing.

if 'i' not in command_args.keys():
    query = 'test_mg.fasta.gz'
else:
    query = command_args['i']
    
if 'o' not in command_args.keys():
    name = 'test_mg'
else:
    name = command_args['o']
    
## Use DIAMOND blastx to search the query against the paprica-mg database.  Note
## that the index values in the resulting DF aren't unique, which works great.
## Note use of --top 0, this means only top tying hits are reported.
    
prot_df = pd.DataFrame.from_csv(ref_dir + 'paprica-mg.prot.csv')
ngenomes = len(prot_df['genome'].unique())
    
print 'executing DIAMOND blastx, this might take a while...'
    
diamond = subprocess.Popen('diamond blastx -k ' + str(ngenomes) + ' --min-score 35 -d ' + ref_dir + 'paprica-mg -q ' + cwd + query + ' -a ' + cwd + name + '.daa', shell = True, executable = executable)
diamond.communicate()

diamond_view = subprocess.Popen('diamond view -a ' + cwd + name + '.daa -o ' + cwd + name + '.txt', shell = True, executable = executable)
diamond_view.communicate()

diamond_df = pd.DataFrame.from_csv(cwd + name + '.txt', sep = '\t', index_col = 0, header = None)
diamond_df.columns = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

hit_tally = pd.DataFrame(diamond_df['sseqid'].value_counts())
hit_tally.columns = ['n_hits']

## Merge hit_tally with prot_df to create a new dataframe, which also eliminates
## all rows that don't have a tally. Write this file out.

annotation_df = pd.merge(prot_df, hit_tally, left_index = True, right_index = True)
annotation_df.to_csv(cwd + name + '.annotation.csv')

## Create a directory to hold the files that are required by pathologic, a
## Genbank file of all proteins with hits for each genome with hits, and an
## entry in genetic-elements.dat for each genome with hits.  This approach
## presents each genome to pathologic as a seperate chromosome.

## Additionally write out a file that is a quantitative list of all the
## EC_numbers.  This allows you to overlay the enzyme abundance data on the
## metabolic overview map in the pathway-tools GUI.

pathos_output_dir = cwd + name + '.pathologic/'
shutil.rmtree(pathos_output_dir, ignore_errors = True)
os.makedirs(pathos_output_dir)

with open(pathos_output_dir + 'genetic-elements.dat', 'w') as all_genetic_elements, open(cwd + name + '.ec_numbers.txt', 'w') as ec_out:
    for genome in annotation_df['genome'].unique():
     
        features = []
        temp = annotation_df[annotation_df['genome'] == genome]
        
        for hit in temp.index:
            
            print hit
            
            qualifiers = {}
            qualifiers['gene'] = hit
            qualifiers['product'] = temp.loc[hit, 'product']
            
            ec_number = temp.loc[hit, 'EC_number']
            ec_number = ec_number.split('|')

            for number in ec_number:
                print >> ec_out, number
                
            qualifiers['EC_number'] = ec_number
        
            start = int(temp.loc[hit, 'start'])
            end = int(temp.loc[hit, 'end'])
        
            location = SeqFeature.FeatureLocation(start, end)
        
            new_feature = SeqFeature.SeqFeature(type = 'CDS', qualifiers = qualifiers)
            new_feature.location = location
            features.append(new_feature)
        
        new_record = SeqRecord(Seq('nnnn', alphabet = IUPAC.ambiguous_dna), id = genome, name = genome, features = features)
        SeqIO.write(new_record, open(pathos_output_dir + name + '.' + genome + '.gbk', 'w'), 'genbank')
        
        print >> all_genetic_elements, 'ID' + '\t' + genome
        print >> all_genetic_elements, 'NAME' + '\t' + genome
        print >> all_genetic_elements, 'TYPE' + '\t' + ':CHRSM'
        print >> all_genetic_elements, 'CIRCULAR?' + '\t' + 'Y'
        print >> all_genetic_elements, 'ANNOT-FILE' + '\t' + name + '.' + genome + '.gbk'
        print >> all_genetic_elements, '//' 

## Now create the organism-params.dat file. 

with open(pathos_output_dir + 'organism-params.dat', 'w') as all_organism_params:
    print >> all_organism_params, 'ID' + '\t' + name
    print >> all_organism_params, 'Storage' + '\t' + 'File'
    print >> all_organism_params, 'Name' + '\t' + name
    print >> all_organism_params, 'Rank' + '\t' + 'Strain'
    print >> all_organism_params, 'Domain' + '\t' + 'TAX-2'
    print >> all_organism_params, 'Create?' + '\t' + 't'
    
## Attempt to remove the old PGDB if it exists, then run pathologic.

shutil.rmtree(pgdb_dir + name + 'cyc', ignore_errors = True)
pathos = subprocess.Popen('pathway-tools -lisp -no-cel-overview -patho ' + pathos_output_dir + ' -disable-metadata-saving &> ' + cwd + 'pathos.' + name + '.log', shell = True, executable = executable) 
pathos.communicate()

## Parse the pathways-report file.   

try:
    with open(pgdb_dir + name + 'cyc/1.0/reports/pathways-report.txt', 'r') as pathways_report, open(cwd + name + '.pathways.txt', 'w') as pathways_out:
        for line in pathways_report:
            if line.startswith('#') == False:
                if line.startswith('Pathway Name') == False:
                    if 'PWY-WAS-NOT-DELETED' not in line:
                        line = line.rstrip()
                        if line != '':
                            line = line.split('|')                                
                            if len(line) > 0:
                                path = line[0]
    
                                print >> pathways_out, path
                                print name, path

except IOError:
    print 'no report found, check for', pgdb_dir + name + 'cyc/1.0/reports/pathways-report.txt'
