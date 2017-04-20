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
    paprica-mgt_build.py [options]
    
OPTIONS:
-ref_dir: The name of the directory containing the paprica database.  Not necessary
if your database is named ref_genome_database (the default).

"""

executable = '/bin/bash' # shell for executing commands

import os
import shutil
import subprocess
import sys
from joblib import Parallel, delayed

from Bio import SeqIO

import pandas as pd
import numpy as np

#%% Function definitions.

## Define a stop function for troubleshooting.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]
    
## Define a function to download viral genomes.
    
def download_assembly(ref_dir_domain, executable, assembly_accession):
    try:
        strain_ftp = genome_data_virus.loc[assembly_accession, 'ftp_path']
        
        mkdir = subprocess.Popen('mkdir ' + ref_dir_domain + 'refseq/' + assembly_accession, shell = True, executable = executable)
        mkdir.communicate()
        
        wget0 = subprocess.Popen('cd ' + ref_dir_domain + 'refseq/' + assembly_accession + ';wget --tries=10 -T30 -q -A "genomic.fna.gz","genomic.gbff.gz","protein.faa.gz" ' + strain_ftp + '/*', shell = True, executable = executable)
        wget0.communicate() 
        
        gunzip = subprocess.Popen('gunzip ' + ref_dir_domain + 'refseq/' + assembly_accession + '/*gz', shell = True, executable = executable)
        gunzip.communicate()
        
        print assembly_accession + ':' + strain_ftp
        
    except KeyError:
        print 'no', assembly_accession, 'path'

#%% Read in command line arguments.

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

paprica_path = os.path.dirname(os.path.abspath("__file__")) + '/' # The location of the actual paprica scripts.  
ref_dir = paprica_path + ref_dir

#%% Download virus sequences, since they aren't used anywhere else.  Since this
## isn't a particularly large database, just overwrite existing.

try:
    shutil.rmtree(ref_dir + 'virus')

except OSError:
    pass

os.mkdir(ref_dir + 'virus')
os.mkdir(ref_dir + 'virus' + '/refseq')

genome_data_virus = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt', header = 1, index_col = 0)
genome_data_virus = genome_data_virus[genome_data_virus.assembly_level == 'Complete Genome']

if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(download_assembly)
    (ref_dir + 'virus/', executable, assembly_accession) for assembly_accession in genome_data_virus.index)

## Check to make sure critical files downloaded.

incomplete_genome = []
    
for assembly_accession in genome_data_virus.index:

    ng_file_count = 0
    temp_faa = ''
    
    ## Genbank now puts some useless fna files in the directory, remove
    ## or they complicate things.
    
    for f in os.listdir(ref_dir + 'virus/' + 'refseq/' + assembly_accession):
        if f.endswith('from_genomic.fna'):
            os.remove(ref_dir + 'virus/' + 'refseq/' + assembly_accession + '/' + f)
            
    ## Now check to make sure that the files you want are in place.
    
    for f in os.listdir(ref_dir + 'virus/' + 'refseq/' + assembly_accession):
        if f.endswith('protein.faa'):
            temp_faa = f
            ng_file_count = ng_file_count + 1
        elif f.endswith('genomic.fna'):
            ng_file_count = ng_file_count + 1
        elif f.endswith('genomic.gbff'):
            ng_file_count = ng_file_count + 1
                    
    if ng_file_count != 3:
        print assembly_accession, 'is missing a Genbank file'
        incomplete_genome.append(assembly_accession)
                       
genome_data_virus = genome_data_virus.drop(incomplete_genome)
genome_data_virus['tax_name'] = genome_data_virus.organism_name.replace(' ', '_')
genome_data_virus.to_csv(ref_dir + 'virus/' + 'genome_data.final.csv')

#%%  Create database.

## Read in genome_data so that you can iterate by genomes that are actually
## used by paprica.

genome_data_bacteria = pd.read_csv(ref_dir + 'bacteria/genome_data.final.csv', index_col = 0, header = 0)
genome_data_bacteria = genome_data_bacteria.dropna(subset = ['clade'])
genome_data_bacteria['domain'] = 'bacteria'

genome_data_archaea = pd.read_csv(ref_dir + 'archaea/genome_data.final.csv', index_col = 0, header = 0)
genome_data_archaea = genome_data_archaea.dropna(subset = ['clade'])
genome_data_archaea['domain'] = 'archaea'

genome_data_eukarya = pd.read_csv(ref_dir + 'eukarya/genome_data.final.csv', index_col = 0, header = 0)
genome_data_eukarya = genome_data_eukarya.dropna(subset = ['clade'])
genome_data_eukarya['domain'] = 'eukarya'

genome_data_virus = pd.read_csv(ref_dir + 'virus/genome_data.final.csv', index_col = 0, header = 0)
genome_data_virus['domain'] = 'virus'

genome_data = pd.concat([genome_data_bacteria, genome_data_archaea, genome_data_eukarya, genome_data_virus])

## Iterate through all the files in refseq and find the gbff files.  First pass
## just counts the number of features with EC_number so that we can create a
## Numpy array the right size.

def count_ec(output, i):
    
    d = genome_data.index[i]
    eci = 0
    domain = genome_data.loc[d, 'domain']
    
    for f in os.listdir(ref_dir + domain + '/refseq/' + d):
        if f.endswith('gbff'):
            
            try:
                for record in SeqIO.parse(ref_dir + domain + '/refseq/' + d + '/' + f, 'genbank'):
                    for feature in record.features:
                        if feature.type == 'CDS':
                            if 'EC_number' in feature.qualifiers.keys():                              
                                eci = eci + 1
                
                output[i] = eci
                print d, 'has', eci, 'features'

            ## For some assemblies an error is raised on a second(?) record identified
            ## in the Genbank file.  It isn't clear why this is happening, pass the error
            ## here.
                                
            except AttributeError:
                pass    

sums = np.memmap(open('tmp.paprica.mmp', 'w+b'), shape = genome_data.index.shape[0], dtype = 'uint64')

Parallel(n_jobs = -1)(delayed(count_ec)(sums, i)
for i in range(0, len(genome_data.index)))
    
eci = int(sums.sum() * 2) # Count_ec is undercounting and not clear why.  Multiply by 2 to insure large enough array.

## Delete mmp

os.remove('tmp.paprica.mmp')
                            
## Create numpy array for data and a 1D array that will become dataframe index.
## You can probably parallelize this as above, but going to take some effort.
                            
prot_array = np.empty((eci,9), dtype = 'object')
prot_array_index = np.empty(eci, dtype = 'object')

## Iterate through all the files in refseq and find the gbk files again.  Store
## the information necessary to create a Genbank record of each feature in the
## array.

i = 0

for d in genome_data.index:
    domain = genome_data.loc[d, 'domain']    
    strain = genome_data.loc[d, 'tax_name']
    
    ## For the domain eukarya, the nuc records are not available as part of the
    ## Genbank format file, they need to be aquired separately and indexed so
    ## that they can be matched with the prot records.
        
    if domain == 'eukarya':
        eukarya_nt_dict = {}
        
        for record in SeqIO.parse(ref_dir + domain + '/refseq/' + d + '/' + d + '.nt.fa', 'fasta'):
            seq_name = str(record.description).split('NCGR_SEQ_ID=')[1]
            seq_name = seq_name.split(' /')[0]
            
            ## Dictionary maps peptide name to nucleotide sequence.
            
            eukarya_nt_dict[seq_name] = str(record.seq)
                
    for f in os.listdir(ref_dir + domain + '/refseq/' + d):
        if f.endswith('gbff'):
            
            try:
                for record in SeqIO.parse(ref_dir + domain + '/refseq/' + d + '/' + f, 'genbank'):
                    for feature in record.features:
                        if feature.type == 'CDS':
                            if 'EC_number' in feature.qualifiers.keys():
                                
                                ## These qualifiers must be present.
                            
                                try:
                                    protein_id = feature.qualifiers['protein_id'][0]
                                    trans = feature.qualifiers['translation'][0]
                                    ec = feature.qualifiers['EC_number']
                                    prod = feature.qualifiers['product'][0]
                                except KeyError:
                                    continue
                                
                                ## Get the nucleotide sequence.
                                
                                if domain != 'eukarya':
                                    
                                    ## Start and end are captured here, but these aren't
                                    ## always correct, which is why I'm now using the feature.extract method.                                    
                                    
                                    start = int(feature.location.start)
                                    end = int(feature.location.end)
                                    seq = str(feature.extract(record).seq)
                                else:
                                    start = 0
                                    end = 0
                                    
                                    try:
                                        seq = eukarya_nt_dict[protein_id.rstrip('_1')]
                                    except KeyError:
                                        seq = 'no_nucleotide_sequence_found'
                                    
                                ## Because enzymes can be identified by combinations of EC numbers,
                                ## need to keep multiple together, and ordered in a consistent way.
                                    
                                if len(ec) > 1:
                                    ec.sort()
                                    ec = ['|'.join(ec)]
                                
                                prot_array_index[i] = protein_id
                                prot_array[i,0] = d
                                prot_array[i,1] = domain
                                prot_array[i,2] = strain
                                prot_array[i,3] = ec[0]
                                prot_array[i,4] = seq
                                prot_array[i,5] = trans
                                prot_array[i,6] = prod
                                prot_array[i,7] = start
                                prot_array[i,8] = end
                                
                                i = i + 1
                                print d, i, 'out of around', int(eci/2), protein_id, ec[0]
                                
            ## For some assemblies an error is raised on a second(?) record identified
            ## in the Genbank file.  It isn't clear why this is happening, pass the error
            ## here.
                                
            except AttributeError:
                pass

## Convert array to pandas dataframe

columns = ['genome', 'domain', 'tax_name', 'EC_number', 'sequence', 'translation', 'product', 'start', 'end']                                    
prot_df = pd.DataFrame(prot_array, index = prot_array_index, columns = columns)

## Determine how often each translation appears.  CDS with a translation
## that is not unique should not be used for taxonomic profiling.  Currently
## not using taxonomic profiling anyway.

prot_counts = pd.DataFrame(prot_df['translation'].value_counts())
prot_counts.columns = ['trans_n_occurrences'] # The number of times that the sequence appears across all genomes.

## Add this information to prot_df.

prot_df = pd.merge(prot_df, prot_counts, left_on = 'translation', right_index = True)

## Determine how often each sequence appears.

nuc_counts = pd.DataFrame(prot_df['sequence'].value_counts())
nuc_counts.columns = ['cds_n_occurrences'] # The number of times that the sequence appears across all genomes.

## Add this information to prot_df.

prot_df = pd.merge(prot_df, nuc_counts, left_on = 'sequence', right_index = True)

## Now that you know which are duplicates, make nonredundant and print out to
## csv.  This is the basis of the paprica-mg database for metagenomic analysis.

prot_nr_trans_df = prot_df.drop_duplicates(subset = ['translation'])
prot_nr_trans_df.to_csv(ref_dir + 'paprica-mg.ec.csv')

## For metatranscriptome analysis we don't just want nonredunant, we want those
## coding sequences that are actually unique.

prot_unique_cds_df = prot_df[prot_df['cds_n_occurrences'] == 1]

## Identical products from genes with silent mutations will still have the same
## protein ID.  Need to eliminate.

prot_unique_cds_df = prot_unique_cds_df[prot_unique_cds_df.sequence != 'no_nucleotide_sequence_found']
prot_unique_cds_df.drop_duplicates(subset = ['translation'], inplace = True)
prot_unique_cds_df.to_csv(ref_dir + 'paprica-mt.ec.csv')

## Make a nonredundant fasta for the metagenomic analysis database.

with open(ref_dir + 'paprica-mg.fasta', 'w') as fasta_out:
    for row in prot_nr_trans_df.iterrows():
        protein_id = row[0]
        translation = row[1]['translation']
        print 'making fasta file', protein_id
        print >> fasta_out, '>' + protein_id
        print >> fasta_out, translation

makedb = subprocess.Popen('diamond makedb --in ' + ref_dir + 'paprica-mg.fasta -d ' + ref_dir + 'paprica-mg', shell = True, executable = executable)
makedb.communicate()  

## Make an unique fasta for the metatranscriptomic analysis database.

with open(ref_dir + 'paprica-mt.fasta', 'w') as fasta_out:
    for row in prot_unique_cds_df.iterrows():
        protein_id = row[0]
        sequence = row[1]['sequence']
        print 'making fasta file', protein_id
        print >> fasta_out, '>' + protein_id
        print >> fasta_out, sequence

makedb = subprocess.Popen('bwa index ' + ref_dir + 'paprica-mt.fasta', shell = True, executable = executable)
makedb.communicate()  
                        
                        
        