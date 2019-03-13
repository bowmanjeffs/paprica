#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Sun Jan 04 17:06:39 2015

@author: Jeff Bowman, jsbowman@ucsd.edu

paprica is licensed under a Creative Commons Attribution-NonCommercial
4.0 International License.  IF you use any portion of paprica in your
work please cite:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described
by Metabolic Structure: A General Framework and Application to a Seasonally
Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic
Peninsula." PloS one 10.8 (2015): e0135868.

If your analysis makes specific use of pplacer, Infernal, pathway-tools, or
any other software please make sure that you also cite the relevant publications.

REQUIRES:
    Files:
        kmer_top_1e5.txt
        bacterial_ssu.cm
        archaea_ssu.cm
        eukarya_ssu.cm

    Programs:
        raxmlHPC-PTHREADS-AVX2
        infernal
        seqmagick

    Python modules:
        pandas
        Bio
        joblib
        Scipy
        Numpy
        
RUN AS:
    python paprica-make_ref.py [options]
    
OPTIONS:
    -cpus: The number of cpus for RAxML to use.
    -domain: Which domain are you analyzing?  Either bacteria, archaea, or eukarya.
    -download: Initiate a fresh download from Genbank?  Either T, F, or test.  Test
    allows you to use the small test set of genomes provided here: http://www.polarmicrobes.org/extras/ref_genome_database.tgz.
    -ref_dir: The name for the database you are building.  The default is "ref_genome_database".

This script must be located in the 'paprica' directory as it makes use of relative
paths.
    
"""

### User setable variables. ###

## If there are assemblies that you would like to exclude from analysis
## put them in the list 'bad' below.  Bedellvibrio, for example, causes
## errors for the placement of some reads.  For the eukarya, several of the
## 18S genes fall well outside the scope of the covariance model, presumably
## due to limitations in the training set.  This all need to be excluded or
## the tree is useless.

bad_bacteria = ['GCF_000691605.1', \
'GCF_000348725.1', \
'GCF_000525675.1', \
'GCF_000317895.1']

bad_archaea = []

bad_eukarya = ['MMETSP0017', \
'MMETSP0027', \
'MMETSP0103', \
'MMETSP0151', \
'MMETSP0200', \
'MMETSP0267', \
'MMETSP0403', \
'MMETSP0409', \
'MMETSP0414', \
'MMETSP0434', \
'MMETSP0448', \
'MMETSP0468', \
'MMETSP0562', \
'MMETSP0595', \
'MMETSP0689', \
'MMETSP0898', \
'MMETSP0918', \
'MMETSP0945', \
'MMETSP1015', \
'MMETSP1074', \
'MMETSP1317', \
'MMETSP1392', \
'MMETSP1446']

### End user setable variables. Mucking with anything below this point might ###
### break the script.  Of course it might also improve it :) ###

import os
import subprocess
from joblib import Parallel, delayed
import gzip
import re
import sys
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio.Seq import Seq

import pandas as pd
import numpy as np
from scipy import spatial

executable = '/bin/bash'
                
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
        
## Define some variables based on these arguments.  If any arguments are
## missing replace with a default value.
        
try:        
    domain = command_args['domain']
except KeyError:
    domain = 'archaea'
try:
    cpus = str(command_args['cpus'])
except KeyError:
    cpus = str(8)
try:
    download = command_args['download']
except KeyError:
    download = 'T'
try:
    ref_dir = command_args['ref_dir']
except KeyError:
    ref_dir = 'ref_genome_database'
        
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/' # The location of the actual paprica scripts.    

ref_dir_domain = paprica_path + ref_dir + domain + '/'

## Select model based on domain specified.
    
if domain == 'bacteria':
    cm = paprica_path + 'models/bacteria_ssu.cm'
elif domain == 'archaea':
    cm = paprica_path + 'models/archaea_ssu.cm'
elif domain == 'eukarya':
    cm = paprica_path + 'models/eukarya_ssu.cm'
else:
    print('Error, you must specify either -domain bacteria, -domain archaea, or -domain eukarya!')
    quit()

#%% Define functions

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])

## Define a function to tally the 5, 4, and 3mers in each faa (i.e. to construct
## compositional vectors).

def calc_vector(path, bins):
    
    k = 5
    k1 = k - 1
    k2 = k - 2

    k1_found_bins = {}
    #k1_used_bins = set()
    k2_found_bins = {}
    #k2_used_bins = set()
    found_bins = {}
    
    print('working on', path)
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
    
    for kmer in list(found_bins.keys()):
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
        
    with gzip.open(ref_dir_domain + 'refseq/' + assembly + '/' + assembly + '_' + str(k) + 'mer_bins.txt.gz', 'wt') as bins_out:
        for each in sorted(bins):
            try:
                print(each + '\t' + str(norm_bins[each]), file = bins_out)
            except KeyError:
                print(each + '\t' + str(0), file = bins_out)

## Define a function to find the 16S rRNA gene in a fna genome file.

def search_16S(directory, d):
    for f in os.listdir(directory + '/' + d):
        if f.endswith('fna'):
            print(f) ## testing
            
            find_16S = subprocess.Popen('cmsearch --cpu 1 --tblout ' + directory + '/' + d + '/' + d + '.16S.hits -A ' + directory + '/' + d + '/' + d + '.16S.sto ' + cm + ' ' + directory + '/' + d + '/' + f, shell = True, executable = executable)
            find_16S.communicate()
                        
            convert = subprocess.Popen('seqmagick convert ' + directory + '/' + d + '/' + d + '.16S.sto ' + directory + '/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable)
            convert.communicate()

## Define a function so that the download of assemblies, if desired, can be parallelized.
## This should reduce a major bottleneck, as I don't think connection speed is the limiting
## factor but the process of establishing the connection.  This allows multiple connections
## to be established simultaneously (note: this is now parallelized)

def download_assembly(ref_dir_domain, executable, assembly_accession):
    try:
        strain_ftp = summary_complete.loc[assembly_accession, 'ftp_path']
        
        ## Required to use proxy server on SIO network.  Bafflingly there is bug in wget that
        ## disallows the use of a wildcard with ftp when a proxy server is used.  So ftp
        ## path must be converted to http.  Note that this requires the use of the -nd flag
        ## because the http protocol will try to create the entire directory structure (grr!)
        
        base_name = strain_ftp.split('/')[-1]

        mkdir = subprocess.Popen('mkdir ' + ref_dir_domain + 'refseq/' + assembly_accession, shell = True, executable = executable)
        mkdir.communicate()
        
        for extension in ['_genomic.fna.gz', '_genomic.gbff.gz', '_protein.faa.gz']:
            wget0 = subprocess.Popen('cd ' + ref_dir_domain + 'refseq/' + assembly_accession + ';wget \
                                     --tries=10 -q -r -nd -T30 -e robots=off ' \
                                     + strain_ftp + '/' + base_name + extension, \
                                     shell = True, executable = executable)
            wget0.communicate() 
        
        gunzip = subprocess.Popen('gunzip ' + ref_dir_domain + 'refseq/' + assembly_accession + '/*', shell = True, executable = executable)
        gunzip.communicate()
        
        print(assembly_accession + ':' + strain_ftp)
        
    except KeyError:
        print('no', assembly_accession, 'path')
        
def download_euks(online_directory):
    try:
        assembly_accession = summary_complete.loc[online_directory, 'sample_name']
        #strain_ftp = 'ftp://ftp.imicrobe.us/projects/104/samples/' + online_directory
        strain_ftp = 'https://de.cyverse.org/anon-files//iplant/home/shared/imicrobe/projects/104/samples/' + online_directory
        
        strain_nt = assembly_accession + '.nt.fa'
        strain_fa = assembly_accession + '.pep.fa'
        
        mkdir = subprocess.Popen('mkdir ' + ref_dir_domain + 'refseq/' + assembly_accession, shell = True, executable = executable)
        mkdir.communicate()
        
        for f in [strain_nt, strain_fa]:
            wget0 = subprocess.Popen('cd ' + ref_dir_domain + 'refseq/' + assembly_accession + ';wget --tries=10 -T30 -q ' + strain_ftp + '/' + f, shell = True, executable = executable)
            wget0.communicate() 
        
        wget1 = subprocess.Popen('cd ' + ref_dir_domain + 'refseq/' + assembly_accession + ';wget --tries=10 -T30 -q ' + strain_ftp + '/annot/swissprot.gff3', shell = True, executable = executable)
        wget1.communicate() 
        
        ## gunzip commands no longer necessary since cyverse does not compress
    
        #gunzip = subprocess.Popen('gunzip ' + ref_dir_domain + 'refseq/' + assembly_accession + '/*', shell = True, executable = executable)
        #gunzip.communicate()
        
        print(assembly_accession + ':' + ref_dir_domain + 'refseq/' + assembly_accession)
        
    except KeyError:
        print('no', online_directory, 'online directory') 
        
## Define function to deterime if all the necessary files are present in existing genome directories
## for bacteria and archaea.

def check_directory(ref_dir_domain, genome):
    file_count = 0
    
    print('Checking files for accession', ref_dir_domain + 'refseq/' + genome)
    
    for f in os.listdir(ref_dir_domain + 'refseq/' + genome):
        if f.endswith('protein.faa'):
            file_count = file_count + 1
        elif f.endswith('genomic.fna'):
            file_count = file_count + 1
        elif f.endswith('genomic.gbff'):
            file_count = file_count + 1
        elif f.endswith('16S.fasta'):
            file_count = file_count + 1
        elif f.endswith('bins.txt.gz'):
            file_count = file_count + 1
        
            ## Sometime empty compositional vectors break the script downstream.  Not sure why this is happening,
            ## test to see if vectors are valid and add to new_genome_faa if not so that they can be recalculated.
        
            test_bins = np.loadtxt(ref_dir_domain + 'refseq/' + genome + '/' + f, usecols = [1])
            if len(test_bins) != 1e5:
                file_count = file_count - 1
            
    if file_count != 5:
        shutil.rmtree(ref_dir_domain + 'refseq/' + genome)

## The download and parsing of eukaryotic genomes is very different from bacterial and archaeal
## genomes and needs to be handled separately.

def get_eukaryotes():

    ## Get eukaryote sample data from MMETSP.
    
    if 'sample-attr.tab.gz' not in os.listdir(ref_dir_domain):
        wget0 = subprocess.Popen('cd ' + ref_dir_domain + ';wget --tries=10 -T30 -q https://de.cyverse.org/anon-files//iplant/home/shared/imicrobe/projects/104/sample-attr.tab', shell = True, executable = executable)
        wget0.communicate()
    
    ## Parse this file into a dataframe.
    
    summary_complete = pd.DataFrame()
    l = 0    
    
    with open(ref_dir_domain + 'sample-attr.tab', 'rb') as sample_attr:
        for line in sample_attr:
            line = line.rstrip()
            line = line.split('\t')
            l = l + 1
            if l != 1:
                summary_complete.loc[line[0], 'sample_name'] = line[1]
                summary_complete.loc[line[0], line[2]] = line[3]   
                        
    ## Get the 18S sequences.

    if 'combined_18S.fasta' not in os.listdir(ref_dir_domain):
        wget_18S = subprocess.Popen('cd ' + ref_dir_domain + ';wget https://de.cyverse.org/anon-files//iplant/home/shared/imicrobe/projects/104/18s/18s.fa;mv 18s.fa combined_18S.fasta', shell = True, executable = executable)
        wget_18S.communicate()
    
    return(summary_complete)    
    
#%% Download fna, gbff, faa files for each completed bacterial or archaeal genome, or execute eukaryote function.

## !! The logic in this section need to be revised.  Overly complex right now.

if download in ['T', 'test']:
    
    ## If the necessary directory structure isn't present, add it.
    
    print('Checking for reference database directories, will create if necessary...')

    try:
        os.listdir(ref_dir)
    except OSError:
        os.mkdir(ref_dir)
        os.mkdir(ref_dir + '/user')
        os.mkdir(ref_dir + '/user/bacteria')
        os.mkdir(ref_dir + '/user/archaea')
        
    try:
        os.listdir(ref_dir + domain)
    except OSError:
        os.mkdir(ref_dir + domain)
        
    if domain in ['bacteria', 'archaea']:
        
        try:
            os.listdir(ref_dir + domain + '/refseq')
        except OSError:
            os.mkdir(ref_dir + domain + '/refseq')
            
    ## Execute some tasks specific to the domain Eukarya.
        
    if domain == 'eukarya':
        
        ## For eukarya only you need to start with an empty refseq directory every
        ## time you want to download. There is no option to just update the database
        
        shutil.rmtree(ref_dir_domain)
        os.mkdir(ref_dir + domain)
        os.mkdir(ref_dir + domain + '/refseq')
        
        summary_complete = get_eukaryotes()
        
        ## Execute the download function.  You can't do this from inside the function
        ## because parallel requires input variables to be global.
        
        if __name__ == '__main__':  
            Parallel(n_jobs = -1, verbose = 5)(delayed(download_euks)
            (online_directory) for online_directory in summary_complete.index)
            
        ## Check to make sure that each downloaded directory has a .fa, .nt, and .gff3 file
        ## extension.  Remove if it does not, and add to bad_eukarya.  The .nt
        ## file is not required for regular paprica, but this makes the database
        ## compatible with paprica-mgt.
            
        for genome in summary_complete.sample_name:
            file_count = 0
            
            try:
                for f in os.listdir(ref_dir_domain + 'refseq/' + genome):
                    if f.endswith('pep.fa'):
                        file_count = file_count + 1
                    if f.endswith('nt.fa'):
                        file_count = file_count + 1
                    elif f.endswith('swissprot.gff3'):
                        file_count = file_count + 1 
            except OSError:
                pass
                    
            if file_count != 3:
                try:
                    shutil.rmtree(ref_dir_domain + 'refseq/' + genome)
                except OSError:
                    pass
                bad_eukarya.append(genome)
            
        ## Remove incomplete downloads from summary_complete.
            
        summary_complete = summary_complete[~summary_complete.sample_name.isin(bad_eukarya)]
        
        ## For Eukarya, dataframe is currently indexed by online directory number.
        ## Downstream scripts need it to be indexed by acccession.  

        summary_complete = summary_complete.set_index('sample_name')

    ## If bacteria or archaea continue here.

    else:
            
        ## Identify which genomes are already contained in the database.  Confirm that these
        ## genome directories have the necessary files and eliminate if they do not.
            
        if __name__ == '__main__':  
            Parallel(n_jobs = -1, verbose = 5)(delayed(check_directory)
            (ref_dir_domain, genome) for genome in os.listdir(ref_dir_domain + 'refseq/'))
                
        old_genomes = pd.Series(os.listdir(ref_dir_domain + 'refseq/'))
        
        ## If download == T, download all the completed genomes, starting with
        ## the Genbank assembly_summary.txt file.
        
        if download == 'T':
            summary = pd.read_table('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/' + domain + '/assembly_summary.txt', header = 1, index_col = 0)
            summary_complete = summary[summary.assembly_level == 'Complete Genome']
        else:
            summary_complete = pd.read_csv(ref_dir_domain + 'genome_data.csv', header = 0, index_col = 0)
        
        ## Drop the bad genomes if new download.
        
        if download == 'T':
            if domain == 'bacteria':
                summary_complete = summary_complete.drop(bad_bacteria)
            elif domain == 'archaea':
                summary_complete = summary_complete.drop(bad_archaea)
            
        ## Determine which genomes need to be downloaded.
            
        new_genomes = summary_complete.index[summary_complete.index.isin(old_genomes) == False]
            
        ## Download the good genomes.
        
        if download == 'T':
        
            if __name__ == '__main__':  
                Parallel(n_jobs = -1, verbose = 5)(delayed(download_assembly)
                (ref_dir_domain, executable, assembly_accession) for assembly_accession in new_genomes)
            
        ## Sometime wget will fail to download a valid file.  This causes problems
        ## downstream.  Check that each directory has an faa and fna file, and remove
        ## from summary_complete if it does not.
    
        new_genome_faa = [] # This will be used for compositional vector creation, and will only hold new genomes with valid faa
        incomplete_genomes = []
            
        for assembly_accession in new_genomes:
    
            ng_file_count = 0
            temp_faa = ''
            
            ## Genbank now puts some useless fna files in the directory, remove
            ## or they complicate things.
            
            for f in os.listdir(ref_dir_domain + 'refseq/' + assembly_accession):
                if f.endswith('from_genomic.fna'):
                    os.remove(ref_dir_domain + 'refseq/' + assembly_accession + '/' + f)
                    
            ## Now check to make sure that the files you want are in place.
            
            for f in os.listdir(ref_dir_domain + 'refseq/' + assembly_accession):
                if f.endswith('protein.faa'):
                    temp_faa = f
                    ng_file_count = ng_file_count + 1
                elif f.endswith('genomic.fna'):
                    ng_file_count = ng_file_count + 1
                elif f.endswith('genomic.gbff'):
                    ng_file_count = ng_file_count + 1
                            
            if ng_file_count != 3:
                print(assembly_accession, 'is missing a Genbank file')
                incomplete_genomes.append(assembly_accession)
            else:
                new_genome_faa.append(ref_dir_domain + 'refseq/' + assembly_accession + '/' + temp_faa)
                               
        summary_complete = summary_complete.drop(incomplete_genomes)
        new_genomes = new_genomes.drop(incomplete_genomes)   
            
## If download == F, assume that directories have already been checked for files and incomplete
## directories removed.  This might be a bad assumption...
    
if download == 'F':
    
    new_genome_faa = []
    
    summary_complete = pd.read_csv(ref_dir_domain + 'genome_data.csv', header = 0, index_col = 0)
    
    if domain != 'eukarya':
        summary_complete = summary_complete[summary_complete.assembly_level != 'Draft']
    
    new_genomes = summary_complete.index
    
    for assembly_accession in new_genomes:
        for f in os.listdir(ref_dir_domain + 'refseq/' + assembly_accession):
            if f.endswith('protein.faa'):
                temp_faa = f
                new_genome_faa.append(ref_dir_domain + 'refseq/' + assembly_accession + '/' + temp_faa)
        
## Add columns to dataframe that will be used later.

summary_complete['n16S'] = np.nan
summary_complete['nge'] = np.nan
summary_complete['ncds'] = np.nan
summary_complete['genome_size'] = np.nan
summary_complete['tax_name'] = np.nan
summary_complete['phi'] = np.nan
summary_complete['GC'] = np.nan

## Write out summary_complete.

summary_complete.to_csv(ref_dir_domain + 'genome_data.csv')

#%% Get the 16S/18S rRNA genes for each assembly and genome parameters.  Eukarya
## are a special case and must be handled separate from bacteria and archaea.

## Eliminate duplicate 18S, then remove 18S associated with incomplete downloads.

if domain == 'eukarya':

    ## It's possible that additional genome were added to bad_eukarya since the
    ## last download.  To avoid needing to download everything just to remove,
    ## check to make sure that summary_complete doesn't contain any bad genomes.

    summary_complete = summary_complete[~summary_complete.index.isin(bad_eukarya)]
    
    ## Duplicate sequences are not allowed as they are incompatible with
    ## phylogenetic placement, and shouldn't really be used for building trees
    ## anyway.  Remove these.
    
    unique = subprocess.Popen('seqmagick convert --deduplicate-sequences ' + ref_dir_domain + 'combined_18S.fasta ' + ref_dir_domain + 'combined_18S.unique.fasta', shell = True, executable = executable)
    unique.communicate()
    
    kept_genomes = [] # This will hold the names of genomes with unique 18S and correct data files (identified earlier).
    
    ## Check to make sure that each sequence in combined_18S.unique.fasta is
    ## present in summary_complete.  If it is not it is bad and should not
    ## be used.  Also add a reasonable taxonomic name for each reference sequence
    ## to summary_complete as tax_id.

    with open(ref_dir_domain + 'combined_18S.unique.fasta', 'r') as fasta_18S, open(ref_dir_domain + 'combined_18S.' + domain + '.tax.fasta', 'w') as good_fasta_18S:
        for line in fasta_18S:
            if line.startswith('>'):
                tax_name = line.strip('>')
                tax_name = tax_name.rstrip()
                genome = tax_name.split('|')[0]
                
                if genome in summary_complete.index:
                    summary_complete.loc[genome, 'tax_name'] = tax_name
                    keep = True
                    kept_genomes.append(genome)
                    print(line, end=' ', file=good_fasta_18S)
                else:
                    keep = False
            else:
                if keep == True:
                    print(line, end=' ', file=good_fasta_18S)
                    
    summary_complete = summary_complete[summary_complete.index.isin(kept_genomes)]
    
    ## Write out summary_complete and exit.
    
    summary_complete.to_csv(ref_dir_domain + 'genome_data.csv') 
    quit()
       
## For bacteria and archaea, find 16S rRNA genes in fna files. Get some paramenters on the genome; number of
## 16S genes, number of elements, size of genome, and add these to summary_complete.
## Generate two fasta files of the 16S rRNA genes.  One will be used later to build
## the reference tree and has sensible taxonomic names.  One is used to calculate
## the phi values and is named by assembly.

## Execute the 16S search function.
            
if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(search_16S)
    (ref_dir_domain + 'refseq', d) for d in new_genomes)

## Iterating by summary_complete.index should eliminate issues with adding 16S
## rRNA sequences for genomes that don't have an faa file.

with open(ref_dir_domain + 'combined_16S.fasta', 'w') as fasta_out:
    for d in summary_complete.index:
        for f in os.listdir(ref_dir_domain + 'refseq/' + d):
            if f.endswith('fna'):
                
                ## Count the number of 16S genes.
                
                count_16S = subprocess.Popen('grep -c \'>\' ' + ref_dir_domain + 'refseq/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable, stdout = subprocess.PIPE)
                n16S = count_16S.communicate()[0]
                n16S = n16S.rstrip()
                n16S = int(n16S)
                
                ## Print one of the 16S rRNA genes to combined_16.
                                
                keep = True
                for record in SeqIO.parse(ref_dir_domain + 'refseq/' + d + '/' + d + '.16S.fasta', 'fasta'):
                    if keep == True:
                        new_record = SeqRecord(record.seq)                        
                        new_record.id = d                        
                        new_record.description = ''
                        
                        SeqIO.write(new_record, fasta_out, 'fasta')                        
                        keep = False
                    else:
                        continue
                        
                ## Count the number of genetic elements.
            
                count_ge = subprocess.Popen('grep -c \'>\' ' + ref_dir_domain + 'refseq/' + d + '/' + f, shell = True, executable = executable, stdout = subprocess.PIPE)
                nge = count_ge.communicate()[0]
                nge = nge.rstrip()
                nge = int(nge)

                ## Get genome size.
                
                genome_size = 0                
                with open(ref_dir_domain + 'refseq/' + d + '/' + f, 'r') as fna_in:
                    for line in fna_in:
                        if line.startswith('>') == False:
                            line = line.rstrip()
                            genome_size = genome_size + len(line)
                            
                ## Add genome data to summary_complete.
                
                summary_complete.loc[d, 'n16S'] = n16S
                summary_complete.loc[d, 'nge'] = nge
                summary_complete.loc[d, 'genome_size'] = genome_size
                
                ## Get GC content.
                
                temp_seq_str = ''
                for record in SeqIO.parse(ref_dir_domain + 'refseq/' + d + '/' + f, 'fasta'):
                    temp_seq_str = temp_seq_str + str(record.seq)
                
                temp_seq = Seq(temp_seq_str)                            
                summary_complete.loc[d, 'GC'] = SeqUtils.GC(temp_seq)
                                   
            elif f.endswith('faa'):
                
                ## Count the number of CDS in the genome.
                
                count_cds = subprocess.Popen('grep -c \'>\' ' + ref_dir_domain + 'refseq/' + d + '/' + f, shell = True, executable = executable, stdout = subprocess.PIPE)
                ncds = count_cds.communicate()[0]
                ncds = ncds.rstrip()
                ncds = int(ncds)
                
                summary_complete.loc[d, 'ncds'] = ncds

#%% Generate a distance matrix for the extracted 16S rRNA genes.
## Align the 16S rRNA genes, first remove existing gaps as the "aligned" sequences
## are not of the same length.

degap = subprocess.Popen('seqmagick mogrify --ungap ' + ref_dir_domain + 'combined_16S.fasta', shell = True, executable = executable)
degap.communicate()

## Duplicate sequences are removed before alignment, may make more sense to remove after.

unique = subprocess.Popen('seqmagick convert --deduplicate-sequences ' + ref_dir_domain + 'combined_16S.fasta ' + ref_dir_domain + 'combined_16S.unique.fasta', shell = True, executable = executable)
unique.communicate()

align_16S = subprocess.Popen('cmalign --outformat Pfam -o ' + ref_dir_domain + 'combined_16S.align.sto ' + cm + ' ' + ref_dir_domain + 'combined_16S.unique.fasta', shell = True, executable = executable)        
align_16S.communicate() 

## Use RAxML to calculate the ML-based distance between taxon pairs.  I thought
## RAxML could handle the sto format but it doesn't seem to like it, so convert
## to fasta first. 

convert = subprocess.Popen('seqmagick convert ' + ref_dir_domain + 'combined_16S.align.sto ' + ref_dir_domain + 'combined_16S.align.fasta', shell = True, executable = executable)
convert.communicate()

remove_dist = subprocess.Popen('rm ' + ref_dir_domain + '*dist', shell = True, executable = executable)
remove_dist.communicate()

dist = subprocess.Popen('cd ' + ref_dir_domain + ';raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -f x -p 12345 -s ' + ref_dir_domain + 'combined_16S.align.fasta -m GTRGAMMA -n dist', shell = True, executable = executable)
dist.communicate()

#%% Calculate the compositional vectors.
## Calculate the compositional vectors for the faa files.  Start by reading in
## the top 100,000 kmers that account for variability between genomes.

bins = set()

with open(paprica_path + 'models/kmer_top_1e5.txt', 'r') as top_kmers:
    for line in top_kmers:
        line = line.rstrip()
        bins.add(line)
    
## Run the composition vector function in parallel.
      
if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(calc_vector)
    (genome_path, bins) for genome_path in new_genome_faa)

#%% Generate symmetrical 16S and cv matrices.
## All the output was saved as a seperate file to facilitate parallel computation.
## Concantenate into common numpy array, then write to a csv.gz file in case you want
## it later.  The current method is a little dangerous because it does not check
## to make sure that the bins are in the right order as they are added to the 
## array, although they should be.  This is done for reasons of speed as a pandas
## dataframe takes orders of magnitude too long to complete these task.

## Get a list of the genomes corresponding to the unique 16S rRNA genes.
                
unique_assembly = []

for record in SeqIO.parse(ref_dir_domain + 'combined_16S.unique.fasta', 'fasta'):
    for f in os.listdir(ref_dir_domain + 'refseq/' + record.id):
        if f.endswith('faa'):
            unique_assembly.append(record.id)
            
unique_assembly = sorted(unique_assembly)

## Read in the compositional vectors belonging to unique assemblies.

table_out = np.zeros(shape = (len(bins), len(unique_assembly)))

i = 0
l = len(unique_assembly)

for d in sorted(unique_assembly):
    print('reading vector matrix,', i + 1, 'of', l)
    temp_bins = np.loadtxt(ref_dir_domain + 'refseq/' + d + '/' + d + '_5mer_bins.txt.gz', usecols = [1])
    table_out[:,i] = temp_bins
    i = i + 1
    
## Currently this matrix is being written out but is not being used anywhere.

table_out_df = pd.DataFrame(table_out, index = sorted(bins), columns = unique_assembly)
table_out_df.to_csv(ref_dir_domain + '5mer_compositional_vectors.csv.gz', compression = 'gzip')

#%% Calculate phi from the two distance measures
## Generate a distance matrix of the genomes according to compositional vectors.

## Note: the distance function will work on a pandas dataframe which will avoid the need
## to transpose if you want to simply execute on table_out_df

dist_cv = spatial.distance.pdist(np.transpose(table_out), metric = 'braycurtis')
dist_cv = spatial.distance.squareform(dist_cv)
dist_cv = pd.DataFrame(dist_cv, columns = unique_assembly, index = unique_assembly)
dist_cv.to_csv(ref_dir_domain + '5mer_compositional_vectors.dist')

## Read in the 16S distance matrix and convert to a symmetrical square matrix.

dist_16S = pd.read_csv(ref_dir_domain + 'RAxML_distances.dist', names = ['taxa1', 'taxa2', 'distance'], delim_whitespace = True)
col_row_names = pd.concat([dist_16S.taxa1, dist_16S.taxa2])
col_row_names = pd.Series(col_row_names.unique())
new_dist_16S = pd.DataFrame(spatial.distance.squareform(dist_16S['distance']), columns = col_row_names, index = col_row_names)

## Normalize so that mean = 0 and variance = 1, then shift so that in the
## range of 0 to 1.  There might be a better way to do this than the notation
## used here (calling on original array, then output array).

print('normalizing distance matrices...')

dist_16S_norm = (new_dist_16S - new_dist_16S.mean().mean()) / new_dist_16S.std().mean()
dist_16S_norm = dist_16S_norm + abs(dist_16S_norm.min().min())
dist_16S_norm = dist_16S_norm / dist_16S_norm.max().max()

dist_cv_norm = (dist_cv - dist_cv.mean().mean()) / dist_cv.std().mean()
dist_cv_norm = dist_cv_norm + abs(dist_cv_norm.min().min())
dist_cv_norm = dist_cv_norm / dist_cv_norm.max().max()

## I'm not comfortable curve fitting in Python, so export a random sample of the
## data to fit a curve to dist_16S_norm as a function of dist_cv_norm.

rmax = 10000

with open(ref_dir_domain + 'distance_measure_subsample.txt', 'w') as dist_sample:
    for r in range(0, rmax):
        r1 = np.random.choice(dist_16S_norm.columns)
        r2 = np.random.choice(dist_16S_norm.index)
        
        rx = dist_cv_norm.loc[r1, r2]
        ry = dist_16S_norm.loc[r1, r2]

        print(r1, r2, rx, ry, file=dist_sample)
        
## I fit an exponential curve: 16S = 0.0005 * exp(7.2259 * cv) at R2 = 0.74.
## These values are used below, but if you're going through all the trouble
## to read this script you should fit your own curve.
        
print('calculating phi parameter...')
        
pred_16S = np.exp(dist_cv_norm.mul(7.2259)).mul(0.0005) # predicted 16S values based on above function
residuals = pred_16S.sub(dist_16S_norm) # difference between predicted and observed
        
phi = residuals.mean()
phi = phi + abs(phi.min())
phi = phi / phi.max()

summary_complete.loc[:, 'phi'] = phi
#summary_complete.loc[phi.idxmax, :] # Use this if you'd like to check the max.

## If there are manually downloaded, user specified draft genomes add these to
## refseq, add entry to summary_complete, and find 16S genes.  If no 16S gene
## found report this and don't include.

## Make sure that each draft genome directory has at a minimum a gbff and fna
## file.

good_drafts = {}

try:
    os.remove(ref_dir + 'user/' + domain + '/' + 'draft.combined_16S.fasta')
except OSError:
    pass

if len(os.listdir(ref_dir + 'user/' + domain)) > 0:
    for assembly_accession in os.listdir(ref_dir + 'user/' + domain):
        
        ## Exception clause is necessary because of .gitignore file, and
        ## possible other old files that might lurk here.        
        
        try:

            dir_contents = os.listdir(ref_dir + 'user/' + domain + '/' + assembly_accession)
            gbff = False
            fna = False
            
            for item in dir_contents:
                if item.endswith('gbff'):
                    gbff = True
                    
                    ## Get a name while you're at it.
    
                    for record in SeqIO.parse(ref_dir + 'user/' + domain + '/' + assembly_accession + '/' + item, 'genbank'):
                        name = record.annotations['organism']
                        
                if item.endswith('fna'):
                    fna = True
                    
            if gbff == False or fna == False:
                print('draft', assembly_accession, 'is missing either fna or gbff')
            else:
                good_drafts[assembly_accession] = re.sub('\s', '_', name)

        except OSError:
            continue
            
## Find the 16S rRNA genes in each good draft, and add an entry to
## summary_complete.  All the entries should be nan because we don't want
## to use the draft values for any predictions.
            
if len(list(good_drafts.keys())) > 0:
        
    if __name__ == '__main__':  
        Parallel(n_jobs = -1, verbose = 5)(delayed(search_16S)
        (ref_dir + 'user/' + domain, d) for d in list(good_drafts.keys()))

with open(ref_dir + 'user/' + domain + '/' + 'draft.combined_16S.fasta', 'w') as draft_fasta_out:        
        for d in list(good_drafts.keys()):
            
            count_16S = subprocess.Popen('grep -c \'>\' ' + ref_dir + 'user/' + domain + '/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable, stdout = subprocess.PIPE)
            n16S = count_16S.communicate()[0]
            n16S = n16S.rstrip()
            n16S = int(n16S)
            
            ## Were any 16S rRNA genes found in the draft assembly?
            
            if n16S > 0:
                
                ## Make an entry for the draft in summary_complete
                
                summary_complete.loc[d, 'organism_name'] = d + '_' + 'DRAFT' + '_' + good_drafts[d]
                summary_complete.loc[d, 'assembly_level'] = 'Draft'
                summary_complete.loc[d, 'tax_name'] = d + '_' + 'DRAFT' + '_' + good_drafts[d]
                
                ## Current behavior is to assign just 1 16S gene copy to draft genomes.  In
                ## the future it will be desirable to estimate from LCA with a completed genome.
                
                summary_complete.loc[d, 'n16S'] = 1
                
                ## Get GC content.
                
                temp_seq_str = ''
                
                for f in os.listdir(ref_dir + 'user/' + domain + '/' + d):
                    if f.endswith('.fna'):
                        for record in SeqIO.parse(ref_dir + 'user/' + domain + '/' + d + '/' + f, 'fasta'):
                            temp_seq_str = temp_seq_str + str(record.seq)
                            
                temp_seq = Seq(temp_seq_str)                            
                summary_complete.loc[d, 'GC'] = SeqUtils.GC(temp_seq)
                
                ## Copy the directory to refseq so that pathway-tools can find it.
                
                cp = subprocess.Popen('cp -r ' + ref_dir + 'user/' + domain + '/' + d + ' ' + ref_dir_domain + 'refseq/', shell = True, executable = executable)
                
                ## Print one of the 16S rRNA genes to combined_16.
                
                keep = True
                                    
                for record in SeqIO.parse(ref_dir + 'user/' + domain + '/' + d + '/' + d + '.16S.fasta', 'fasta'):
                    if keep == True:
                        
                        new_record = SeqRecord(record.seq)                        
                        new_record.id = summary_complete.loc[d, 'organism_name']                      
                        new_record.description = ''
                        
                        SeqIO.write(new_record, draft_fasta_out, 'fasta')                        
                        keep = False
                        
                    else:
                        continue
                    
            else:
                print('sorry, no 16S rRNA genes found in draft assembly', d)
         
## Generate a fasta file with meaningful taxonomic names consisting of only
## those 16S rRNA genes used in the final alignment.
         
unique_16S = set()
        
print('writing out data files')
        
with open(ref_dir_domain + 'combined_16S.' + domain + '.tax.fasta', 'w') as tax_fasta_out:
    for record in SeqIO.parse(ref_dir_domain + 'combined_16S.unique.fasta', 'fasta'):
        
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
        
        unique_16S.add(str(record.seq))
        
    for record in SeqIO.parse(ref_dir + 'user/' + domain + '/' + 'draft.combined_16S.fasta', 'fasta'):        
        if str(record.seq) not in unique_16S:
            SeqIO.write(record, tax_fasta_out, 'fasta')
            unique_16S.add(str(record.seq))
        else:
            print('sorry, a 16S rRNA gene sequence identical to', record.id, 'is already in use, not including')
            summary_complete = summary_complete[summary_complete.organism_name != record.id]
        
## Write out the final summary_complete file.

summary_complete.to_csv(ref_dir_domain + 'genome_data.csv')
        
        
