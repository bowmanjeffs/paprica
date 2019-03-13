#!/usr/bin/env python2
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
    python paprica-mg_run.py [OPTIONS]
    
OPTIONS:
    -i: Input fasta or fasta.gz.
    -o: Prefix for output files.
    -ref_dir: Name of the directory containing the paprica database.
    -pgdb_dir: The location where pathway-tools stores PGDBs.  Only necessary
        if pathways are predicted.
    -pathways: T or F, whether pathways should be predicted.  This can take a
        long time.

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
    paprica_path = os.path.dirname(os.path.realpath(__file__)) + '/' # The location of the actual paprica scripts.

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
    query = 'ERR318619_1.qc.fasta.gz'
else:
    query = command_args['i']
    
if 'o' not in command_args.keys():
    name = 'test_mg'
else:
    name = command_args['o']
    
if 'ref_dir' not in command_args.keys():
    ref_dir = 'paprica-mgt.database/ref_genome_database'
else:
    ref_dir = 'paprica-mgt.database/' + command_args['ref_dir']
    
if 'pathways' not in command_args.keys():
    pathways = 'F'
else:
    pathways = command_args['pathways']
    
if 'pgdb_dir' not in command_args.keys():
    pgdb_dir = '/volumes/hd2/ptools-local/user/'
else:
    pgdb_dir = command_args['pgdb_dir']
    
## Define path to reference directory.
    
ref_dir_path = paprica_path + ref_dir

if ref_dir_path.endswith('/') == False:
    ref_dir_path = ref_dir_path + '/'
    
## Use DIAMOND blastx to search the query against the paprica-mg database.  Note
## that the index values in the resulting DF aren't unique, which works great.
## Note use of --top 0, this means only top tying hits are reported.
    
ec_df = pd.DataFrame.from_csv(ref_dir_path + 'paprica-mg.ec.csv')
ngenomes = len(ec_df['genome'].unique())
    
print 'executing DIAMOND blastx, this might take a while...'
    
## If pathways aren't being predicted than we only want a single hit for each
## query.  This isn't exactly correct, as each hit could have multiple
## correct EC numbers, but is a good compromise for speed and efficiency.
## Instead of trying to make the massive redundant annotation file used for
## pathway prediction nonredundant it is much simpler to just run diamond
## saving only a single hit.
    
diamond = subprocess.Popen('diamond blastx -k 1 --min-score 35 -d ' + ref_dir_path + 'paprica-mg -q ' + cwd + query + ' -a ' + cwd + name + '.paprica-mg.nr.daa', shell = True, executable = executable)
diamond.communicate()

diamond_view = subprocess.Popen('diamond view -a ' + cwd + name + '.paprica-mg.nr.daa -o ' + cwd + name + '.paprica-mg.nr.txt', shell = True, executable = executable)
diamond_view.communicate()

diamond_df = pd.read_csv(cwd + name + '.paprica-mg.nr.txt', sep = '\t', header = None, index_col = 0)
diamond_df.columns = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
hit_tally_nr = pd.DataFrame(diamond_df.sseqid.value_counts())

hit_tally_nr.columns = ['nr_hits'] # The nr_annotation column is the number of reads whose best match is the enzyme

###### Ideally at this point you would have a tally of EC numbers, not just seqids...

## Merge hit_tally with ec_df to create a new dataframe, which also eliminates
## all rows that don't have a tally.

annotation_df = pd.merge(ec_df, hit_tally_nr, left_index = True, right_index = True)

## Make a new dataframe with the number of occurrences of each EC number in the
## metagenome.  In this dataframe the protein id and genome of origin are only
## examples, that is, the hits may have been to many different proteins of this
## same EC number.  The nr_hits and EC_number columns together define the
## metabolic structure of the metagenome.

ec_tally = annotation_df.groupby('EC_number')['nr_hits'].sum()

if pathways == 'T':
    
    warning = """WARNING: You've elected to predict pathways on the metagenome
    with the -pathways T option.  That's great, but beware that this will
    require considerable time, memory, and hard drive space.  Please monitor
    your resources carefully.  In addition please make sure that you maintain
    the connection with your server so that pathway-tools can send you
    graphical messages"""
    
    print warning
    
    print 'running DIAMOND round 2 for redundand output, this will take a little while...'
    
    diamond = subprocess.Popen('diamond blastx -k ' + str(ngenomes) + ' --min-score 35 -d ' + ref_dir_path + 'paprica-mg -q ' + cwd + query + ' -a ' + cwd + name + '.paprica-mg.daa', shell = True, executable = executable)
    diamond.communicate()
    
    diamond_view = subprocess.Popen('diamond view -a ' + cwd + name + '.paprica-mg.daa -o ' + cwd + name + '.paprica-mg.txt', shell = True, executable = executable)
    diamond_view.communicate()
    
    print 'parsing redundant DIAMOND output, this will probably also take a while...'
    
    ## The DIAMOND output can be huge, necessary to read in in chunks.
    
    count_lines = subprocess.Popen('grep -c \'.\' ' + cwd + name + '.paprica-mg.txt',
                                   shell = True,
                                   executable = executable, stdout = subprocess.PIPE)
        
    n_lines = int(count_lines.communicate()[0].rstrip())
    max_chunk = n_lines/100000
    
    diamond_it = pd.read_csv(cwd + name + '.paprica-mg.txt', sep = '\t',
                             header = None,
                             chunksize = 100000,
                             index_col = 0)
                                                
    ichunk = 0                         
    for chunk in diamond_it:
        chunk.columns = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        ichunk += 1
        
        if ichunk == 1:
            temp = pd.DataFrame(chunk['sseqid'].value_counts())
            hit_tally = temp
        else:        
            hit_tally = hit_tally.add(pd.DataFrame(chunk['sseqid'].value_counts()), fill_value = 0)
        
        print ichunk, 'of', max_chunk
    
    hit_tally.columns = ['r_hits'] ## r hits reflects all significant hits to the database
    
    annotation_df = pd.merge(ec_df, hit_tally, left_index = True, right_index = True)
    annotation_df = pd.merge(annotation_df, hit_tally_nr, left_index = True, right_index = True)

    ## Create a directory to hold the files that are required by pathologic, a
    ## Genbank file of all enzymes with hits for each genome with hits, and an
    ## entry in genetic-elements.dat for each genome with hits.  This approach
    ## presents each genome to pathologic as a seperate chromosome.
    
    ## Additionally write out a file that is a quantitative list of all the
    ## EC_numbers.  This allows you to overlay the enzyme abundance data on the
    ## metabolic overview map in the pathway-tools GUI.
        
    pathos_output_dir = cwd + name + '.pathologic/'
    shutil.rmtree(pathos_output_dir, ignore_errors = True)
    os.makedirs(pathos_output_dir)
    
    print 'generating files and directory structure for pathway-tools'
    
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
                
                ## It isn't necessray to split by | anymore, because
                ## there is only one EC number allowed per line in the
                ## database csv.  However this shouldn't be causing any
                ## problems.
    
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
        
    print 'executing pathway-tools, this will take several hours'
    
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

## No need to include the translation column, which takes up excessive space.
        
annotation_df.drop('translation', 1).to_csv(cwd + name + '.annotation.csv')
ec_tally.to_csv(cwd + name + '.mg.sum_ec.csv')
    
