# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 17:06:39 2015

@author: jeff

REQUIRES:
    Files:
        kmer_top_1e5.txt
        bacterial_ssu.cm

    Programs:
        raxmlHPC-PTHREADS-AVX2
        infernal
        seqmagick

    Python modules:
        pandas
        Bio
        joblib
        
RUN AS:
    python paprica_make_ref_v0.20.py
    
"""

### User setable variables. ###

## Set to true to initiate fresh download of genomes.
download = True

## If there are assemblies that you would like to exclude from analysis
## put them in the list 'bad' below.  Bedellvibrio, for example, causes
## errors for the placement of some reads.

bad = ['GCF_000691605.1', \
'GCF_000348725.1', \
'GCF_000525675.1', \
'GCF_000317895.1']

### End user setable variables. Mucking with anything below this point might ###
### break the script.  Of course it might also improve it :) ###

import os
import subprocess
from joblib import Parallel, delayed
import gzip
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils

import pandas as pd
import numpy as np
from scipy import spatial
import math

executable = '/bin/bash'

## Read in profile.  Required variables are ref_dir and cutoff. ###

variables = {}

def get_variable(line, variables):
    
    line = line.rstrip()
    line = line.split('=')
    variable = line[0]
    value = line[1]
    
    variables[variable] = value
    return variables

with open('paprica_profile.txt', 'r') as profile:
    for line in profile:
        if line.startswith('#') == False:
            if line != '\n':
                get_variable(line, variables)

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]

#%% Download fna, gbff, faa files for each completed bacterial genome.

## Define a function so that the download of assemblies, if desired, can be parallelized.
## This should reduce a major bottleneck, as I don't think connection speed is the limiting
## factor but the process of establishing the connection.  This allows multiple connections
## to be established simultaneously.

def download_assembly(ref_dir, executable, assembly_accession):
    try:
        strain_ftp = summary_complete.loc[assembly_accession, 'ftp_path']
        
        subprocess.call('mkdir ' + variables['ref_dir'] + 'refseq/' + assembly_accession + ';cd ' + variables['ref_dir'] + 'refseq/' + assembly_accession + ';wget -A genomic.fna.gz ' + strain_ftp + '/*', shell = True, executable = executable)
        wget1 = subprocess.Popen('cd ' + variables['ref_dir'] + 'refseq/' + assembly_accession + ';wget --tries=3 -T30 -A genomic.gbff.gz ' + strain_ftp + '/*', shell = True, executable = executable)
        wget1.communicate()
        wget2 = subprocess.Popen('cd ' + variables['ref_dir'] + 'refseq/' + assembly_accession + ';wget --tries=3 -T30 -A protein.faa.gz ' + strain_ftp + '/*', shell = True, executable = executable)
        wget2.communicate()
        gunzip = subprocess.Popen('gunzip ' + variables['ref_dir'] + 'refseq/' + assembly_accession + '/*gz', shell = True, executable = executable)
        gunzip.communicate()
        
    except KeyError:
        print 'no', assembly_accession, 'path'

if download == True:
    
    ## Remove old reference directory
    
    subprocess.call('rm -r ' + variables['ref_dir'], shell = True, executable = executable)
    subprocess.call('mkdir ' + variables['ref_dir'], shell = True, executable = executable)
    subprocess.call('mkdir ' + variables['ref_dir'] + '/refseq', shell = True, executable = executable)
    
    ## Download all the completed genomes, starting with the Genbank assembly_summary.txt file.
    
    summary = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt', header = 0, index_col = 0)
    summary_complete = summary[summary.assembly_level == 'Complete Genome']
    summary_complete = summary_complete.drop(bad)
    
    if __name__ == '__main__':  
        Parallel(n_jobs = -1, verbose = 5)(delayed(download_assembly)
        (variables['ref_dir'], executable, assembly_accession) for assembly_accession in summary_complete.index)
    
    ## Add columns to dataframe that will be used later.
    
    summary_complete['n16S'] = np.nan
    summary_complete['nge'] = np.nan
    summary_complete['ncds'] = np.nan
    summary_complete['genome_size'] = np.nan
    summary_complete['tax_name'] = np.nan
    summary_complete['phi'] = np.nan
    summary_complete['GC'] = np.nan
    
    summary_complete.to_csv(variables['ref_dir'] + 'genome_data.csv')
    
else:
    summary_complete = pd.DataFrame.from_csv(variables['ref_dir'] + 'genome_data.csv', header = 0, index_col = 0)

#%% Get the 16S rRNA genes for each assembly and genome parameters.
       
## Find 16S rRNA genes in fna files. Get some paramenters on the genome; number of
## 16S genes, number of elements, size of genome, and add these to summary_complete.
## Generate two fasta files of the 16S rRNA genes.  One will be used later to build
## the reference tree and has sensible taxonomic names.  One is used to calculate
## the phi values and is named by assembly.

## Define a function to conduct the 16S rRNA gene search.

def search_16S(ref_dir, d):
    for f in os.listdir(ref_dir + 'refseq/' + d):
        if f.endswith('fna'):
            
            find_16S = subprocess.Popen('cmsearch --cpu 1 --tblout ' + ref_dir + 'refseq/' + d + '/' + d + '.16S.hits -A ' + ref_dir + 'refseq/' + d + '/' + d + '.16S.sto ' + variables['cm'] + ' ' + ref_dir + 'refseq/' + d + '/' + f, shell = True, executable = executable)
            find_16S.communicate()
                        
            convert = subprocess.Popen('seqmagick convert ' + ref_dir + 'refseq/' + d + '/' + d + '.16S.sto ' + ref_dir + 'refseq/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable)
            convert.communicate()   
            
if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(search_16S)
    (variables['ref_dir'], d) for d in os.listdir(variables['ref_dir'] + 'refseq'))
        
faa_paths = [] # This is used during compositional vector calculation.

with open(variables['ref_dir'] + 'combined_16S.fasta', 'w') as fasta_out:
    for d in os.listdir(variables['ref_dir'] + 'refseq'):
        for f in os.listdir(variables['ref_dir'] + 'refseq/' + d):
            if f.endswith('fna'):
                
                ## Count the number of 16S genes.
                
                count_16S = subprocess.Popen('grep -c \'>\' ' + variables['ref_dir'] + 'refseq/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable, stdout = subprocess.PIPE)
                n16S = count_16S.communicate()[0]
                n16S = n16S.rstrip()
                n16S = int(n16S)
                
                ## Print one of the 16S rRNA genes to combined_16.
                                
                keep = True
                for record in SeqIO.parse(variables['ref_dir'] + 'refseq/' + d + '/' + d + '.16S.fasta', 'fasta'):
                    if keep == True:
                        new_record = SeqRecord(record.seq)                        
                        new_record.id = d                        
                        new_record.description = ''
                        
                        SeqIO.write(new_record, fasta_out, 'fasta')                        
                        keep = False
                    else:
                        continue
                        
                ## Count the number of genetic elements.
            
                count_ge = subprocess.Popen('grep -c \'>\' ' + variables['ref_dir'] + 'refseq/' + d + '/' + f, shell = True, executable = executable, stdout = subprocess.PIPE)
                nge = count_ge.communicate()[0]
                nge = nge.rstrip()
                nge = int(nge)

                ## Get genome size.
                
                gs = ''                
                with open(variables['ref_dir'] + 'refseq/' + d + '/' + f, 'r') as fna_in:
                    for line in fna_in:
                        if line.startswith('>') == False:
                            line = line.rstrip()
                            gs = gs + line
                genome_size = len(gs)
                
                summary_complete.loc[d, 'n16S'] = n16S
                summary_complete.loc[d, 'nge'] = nge
                summary_complete.loc[d, 'genome_size'] = genome_size
                
                ## Get GC content
                
                temp_gc = []
                for record in SeqIO.parse(variables['ref_dir'] + 'refseq/' + d + '/' + f, 'fasta'):
                    temp_gc.append(SeqUtils.GC(record.seq))
                if len(temp_gc) > 1:
                    GC = float(sum(temp_gc) / len(temp_gc))
                    summary_complete.loc[d, 'GC'] = GC
                else:
                    summary_complete.loc[d, 'GC'] = temp_gc[0]
                                    
            elif f.endswith('faa'):
                
                ## Count the number of CDS in the genome.
                
                count_cds = subprocess.Popen('grep -c \'>\' ' + variables['ref_dir'] + 'refseq/' + d + '/' + f, shell = True, executable = executable, stdout = subprocess.PIPE)
                ncds = count_cds.communicate()[0]
                ncds = ncds.rstrip()
                ncds = int(ncds)
                
                summary_complete.loc[d, 'ncds'] = ncds

#%% Generate a distance matrix for the extracted 16S rRNA genes.
## Align the 16S rRNA genes, first remove existing gaps as the "aligned" sequences
## are not of the same length.

degap = subprocess.Popen('seqmagick mogrify --ungap ' + variables['ref_dir'] + 'combined_16S.fasta', shell = True, executable = executable)
degap.communicate()

unique = subprocess.Popen('seqmagick convert --deduplicate-sequences ' + variables['ref_dir'] + 'combined_16S.fasta ' + variables['ref_dir'] + 'combined_16S.unique.fasta', shell = True, executable = executable)
unique.communicate()

align_16S = subprocess.Popen('cmalign --outformat Pfam -o ' + variables['ref_dir'] + 'combined_16S.align.sto ' + variables['cm'] + ' ' + variables['ref_dir'] + 'combined_16S.unique.fasta', shell = True, executable = executable)        
align_16S.communicate() 

## Use RAxML to calculate the ML-based distance between taxon pairs.  I thought
## RAxML could handle the sto format but it doesn't seem to like it, so convert
## to fasta first. 

convert = subprocess.Popen('seqmagick convert ' + variables['ref_dir'] + 'combined_16S.align.sto ' + variables['ref_dir'] + 'combined_16S.align.fasta', shell = True, executable = executable)
convert.communicate()

remove_dist = subprocess.Popen('rm ' + variables['ref_dir'] + '*dist', shell = True, executable = executable)
remove_dist.communicate()

dist = subprocess.Popen('cd ' + variables['ref_dir'] + ';raxmlHPC-PTHREADS-AVX2 -T 2 -f x -p 12345 -s ' + variables['ref_dir'] + 'combined_16S.align.fasta -m GTRGAMMA -n dist', shell = True, executable = executable)
dist.communicate()

#%% Calculate the compositional vectors.
## Calculate the compositional vectors for the faa files.  Start by reading in
## the top 100,000 kmers that account for variability between genomes.

bins = set()

with open('kmer_top_1e5.txt', 'r') as top_kmers:
    for line in top_kmers:
        line = line.rstrip()
        bins.add(line)
        
## Define a function to tally the 5, 4, and 3mers in each faa.

def calc_vector(path, bins):
    
    k = 5
    k1 = k - 1
    k2 = k - 2

    k1_found_bins = {}
    #k1_used_bins = set()
    k2_found_bins = {}
    #k2_used_bins = set()
    found_bins = {}
    
    print 'working on', path
    path_split = path.split('/')
    assembly = path_split[-2]
    
    for record in SeqIO.parse(path, 'fasta'):
        query = str(record.seq)
        
        ## k1 and k2
        
        for i in range(0,len(query)):
            kmer = query[i:i + k1]
            try:
                k1_found_bins[kmer] = k1_found_bins[kmer] + 1
            except KeyError:
                k1_found_bins[kmer] = 1
            
        for i in range(0,len(query)):
            kmer = query[i:i+k2]
            try:
                k2_found_bins[kmer] = k2_found_bins[kmer] + 1
            except KeyError:
                k2_found_bins[kmer] = 1
                
        ## k
            
        for i in range(0,len(query)):
            kmer = query[i:i+k]
            try:
                found_bins[kmer] = found_bins[kmer] + 1
            except KeyError:
                found_bins[kmer] = 1
                
    ## k0 - calculate the normalized kmer abundance
        
    norm_bins = {}
    
    for kmer in found_bins.keys():
        if kmer in bins: # if not a kmer of interest, don't bother normalizing
            kmer_1 = kmer[0:-1]
            kmer_2 = kmer[1:]
            kmer_3 = kmer[1:-1]
            bigL = len(query)
            
            kmer_0 = ((k1_found_bins[kmer_1] * k1_found_bins[kmer_2])
            / float(k2_found_bins[kmer_3])) * (((bigL - k + 1) * (bigL - k + 3))
            / float((bigL - k + 2) ** 2))
            
            kmer_norm = (found_bins[kmer] - kmer_0) / kmer_0
            norm_bins[kmer] = kmer_norm
        
    with gzip.open(variables['ref_dir'] + 'refseq/' + assembly + '/' + assembly + '_' + str(k) + 'mer_bins.txt.gz', 'wb') as bins_out:
        for each in sorted(bins):
            try:
                print >> bins_out, each + '\t' + str(norm_bins[each])
            except KeyError:
                print >> bins_out, each + '\t' + str(0)
                
## Get a list of the genomes corresponding to the unique 16S rRNA genes.
                
unique_genomes = []
unique_assembly = []

for record in SeqIO.parse(variables['ref_dir'] + 'combined_16S.unique.fasta', 'fasta'):
    for f in os.listdir(variables['ref_dir'] + 'refseq/' + record.id):
        if f.endswith('faa'):
            unique_genomes.append(variables['ref_dir'] + 'refseq/' + record.id + '/' + f)
            unique_assembly.append(record.id)
    
## Run the composition vector function in parallel.
    
unique_genomes = sorted(unique_genomes)
unique_assembly = sorted(unique_assembly)
  
if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(calc_vector)
    (genome, bins) for genome in unique_genomes)

#%% Generate symmetrical 16S and cv matrices.
## All the output was saved as a seperate file to facilitate parallel computation.
## Concantenate into common numpy array, then write to a csv.gz file in case you want
## it later.  The current method is a little dangerous because it does not check
## to make sure that the bins are in the right order as they are added to the 
## array, although they should be.  This is done for reasons of speed as a pandas
## dataframe takes orders of magnitude too long to complete these task.

table_out = np.zeros(shape = (len(bins), len(unique_assembly)))

i = 0
l = len(unique_assembly)

for d in sorted(unique_assembly):
    print 'writing vector matrix', i + 1, 'of', l
    temp_bins = np.loadtxt(variables['ref_dir'] + 'refseq/' + d + '/' + d + '_5mer_bins.txt.gz', usecols = [1])
    table_out[:,i] = temp_bins
    i = i + 1

table_out_df = pd.DataFrame(table_out, index = sorted(bins), columns = unique_assembly)
table_out_df.to_csv(variables['ref_dir'] + '5mer_compositional_vectors.csv.gz')

#%% Calculate phi from the two distance measures
## Generate a distance matrix of the genomes according to compositional vectors.

dist_cv = spatial.distance.pdist(np.transpose(table_out), metric = 'braycurtis')
dist_cv = spatial.distance.squareform(dist_cv)
dist_cv = pd.DataFrame(dist_cv, columns = unique_assembly, index = unique_assembly)
dist_cv.to_csv(variables['ref_dir'] + '5mer_compositional_vectors.dist')

## Define a function to convert a column matrix to a square matrix.  The current
## function is pretty slow.

def col2square(column):   
    ## Input column matrix must be pandas dataframe in column order index, taxa, taxa, distance.
    column.columns = ['taxa1', 'taxa2', 'distance']
    col_row_names = pd.concat([column.taxa1, column.taxa2])
    col_row_names = col_row_names.unique()
    square = pd.DataFrame(index = col_row_names, columns = col_row_names)
    
    matrix_size = (len(col_row_names) * (len(col_row_names) + 1)) / 2 # n(n + 1) / 2
    
    for index in column.index:
        print 'generating square matrix', index, 'of', matrix_size
        taxa1 = column.taxa1[index]
        taxa2 = column.taxa2[index]
        dist = column.distance[index]
        
        square.loc[taxa1, taxa2] = dist
        square.loc[taxa2, taxa1] = dist
        
    return square

## Read in the 16S distance matrix and convert to a symmetrical square matrix.

dist_16S = pd.read_table(variables['ref_dir'] + 'RAxML_distances.dist', names = ['taxa1', 'taxa2', 'distance'], delim_whitespace = True)
new_dist_16S = col2square(dist_16S)
new_dist_16S = new_dist_16S.fillna(0)

## Normalize so that mean = 0 and variance = 1, then shift so that in the
## range of 0 to 1.  There might be a better way to do this than the notation
## used here (calling on original array, then output array).

dist_16S_norm = (new_dist_16S - new_dist_16S.mean().mean()) / new_dist_16S.std().mean()
dist_16S_norm = dist_16S_norm + abs(dist_16S_norm.min().min())
dist_16S_norm = dist_16S_norm / dist_16S_norm.max().max()

dist_cv_norm = (dist_cv - dist_cv.mean().mean()) / dist_cv.std().mean()
dist_cv_norm = dist_cv_norm + abs(dist_cv_norm.min().min())
dist_cv_norm = dist_cv_norm / dist_cv_norm.max().max()

## I'm not comfortable curve fitting in Python, so export a random sample of the
## data to fit a curve to dist_16S_norm as a function of dist_cv_norm.

rmax = 10000

with open(variables['ref_dir'] + 'distance_measure_subsample.txt', 'w') as dist_sample:
    for r in range(0, rmax):
        r1 = np.random.choice(dist_16S_norm.columns)
        r2 = np.random.choice(dist_16S_norm.index)
        
        rx = dist_cv_norm.loc[r1, r2]
        ry = dist_16S_norm.loc[r1, r2]

        print >> dist_sample, r1, r2, rx, ry
        
## I fit an exponential curve: 16S = 0.0005 * exp(7.2259 * cv) at R2 = 0.74.
## These values are used below, but if you're going through all the trouble
## to run this scrip you should fit your own curve.
        
## Calculate the residuals.  This loop is pretty slow.
        
residuals = pd.DataFrame(index = dist_16S_norm.index, columns = dist_16S_norm.columns)

i = 0
l = str(len(dist_16S_norm.index))
        
for index in dist_16S_norm.index:
    i = i + 1
    for column in dist_16S_norm.columns:
        real_16S = dist_16S_norm.loc[index, column]
        real_cv = dist_cv_norm.loc[index, column]
        pred_16S = 0.0005 * math.exp(real_cv * 7.2259)
        residual = pred_16S - real_16S
        
        residuals.loc[index, column] = residual
        residuals.loc[column, index] = residual
        
    print 'completed phi calculation for', index, str(i), 'of', l
        
phi = residuals.mean()
phi = phi + abs(phi.min())
phi = phi / phi.max()

summary_complete.loc[:, 'phi'] = phi
#summary_complete.loc[phi.idxmax, :] # Use this if you'd like to check the max.
        
## Generate a fasta file with meaningful taxonomic names consisting of only
## those 16S rRNA genes used in the final alignment.
        
with open(variables['ref_dir'] + 'combined_16S.tax.fasta', 'w') as tax_fasta_out:
    for record in SeqIO.parse(variables['ref_dir'] + 'combined_16S.unique.fasta', 'fasta'):
        try:
            tax_name = record.id + '_' + summary_complete.loc[record.id, 'organism_name'] + '_' + summary_complete.loc[record.id, 'infraspecific_name']
        except TypeError:
            tax_name = record.id + '_' + summary_complete.loc[record.id, 'organism_name']
            
        tax_name = re.sub('\[', '', tax_name)
        tax_name = re.sub('\]', '', tax_name)
        
        new_record_tax = SeqRecord(record.seq)        
        new_record_tax.id = tax_name
        new_record_tax.description = ''
        
        SeqIO.write(new_record_tax, tax_fasta_out, 'fasta')
        summary_complete.loc[record.id, 'tax_name'] = tax_name
        
## Write out the final summary_complete file.

summary_complete.to_csv(variables['ref_dir'] + 'genome_data.csv')
        
        