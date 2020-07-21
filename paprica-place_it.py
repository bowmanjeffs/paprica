#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Sat Jan 03 08:59:11 2015

@author: Jeff Bowman, jsbowman@ucsd.edu

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
    Files:
        bacterial_ssu.cm
        
    Programs:
        infernal
        seqmagick
        pplacer (guppy)
        raxml-ng
        
    Python modules:
        Bio
        joblib
        ete3
        pandas

CALL AS:
    paprica-place_it.py -domain [domain] -query [query] -ref [ref] -splits [splits] -n [nseqs] for analysis or
    paprica-place_it.py -domain [domain] -ref [ref] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (which must be .fasta).
    
OPTIONS:
    -domain: The domain (bacteria or archaea) you are analyzing for.
    -large: If included will increase the --mxsize flag in cmalign
        to 4000 Mb.  Note that this requires 4000 Mb per thread, so you
        will need to have that much memory!
    -n:  The number of reads to subsample your input query to.
    -query: The prefix of the query fasta (entire name without extension).
    -ref: The name of the reference package for pplacer.
    -ref_dir: The directory containing the paprica database.
    -splits: The number of files to split the query fasta into to facilitate
        parallel analysis with pplacer.

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
executable = '/bin/bash' # shell for executing commands

import re
import subprocess
import sys
import os
from joblib import Parallel, delayed
import multiprocessing
import datetime
import random
import pandas as pd
import json

## Assume the system is using hyperthreading.

hyperthreading = True

## Determine path to paprica scripts.

try:
    paprica_path = os.path.dirname(os.path.realpath(__file__)) + '/' # The location of the actual paprica scripts.
except NameError:
    paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/'
cwd = os.getcwd() + '/' # The current working directory.
                
## Parse command line arguments.  Arguments that are unique to a run,
## such as the number of splits, should be specified in the command line and
## not in paprica_profile.txt
                
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

## Parse command line for mandatory variables, providing defaults if they are
## not present.  Defaults are generally for testing.

try:
    ref_dir = paprica_path + command_args['ref_dir']  # The complete path to the reference directory being used for analysis.        
except KeyError:
    ref_dir = paprica_path + 'ref_genome_database/'
try:    
    domain = command_args['domain']  # The domain being used for analysis.
except KeyError:
    domain = 'archaea'
try:
    ref = command_args['ref']  # The name of the reference package being used.
except KeyError:
    ref = 'combined_16S.' + domain +'.tax'

## If sys.argv == 1, you are probably running inside Python in testing mode.
## Provide some default values to make this possibe.  If > 1, parse command
## line for optional variables, providing reasonable defaults if not present.
## No default is currently provided for query.
    
if len(sys.argv) == 1:
    splits = 1
#    query = 'test.' + domain
#    command_args['query'] = query

else:
    
    try:    
        splits = int(command_args['splits']) # The number of splits to make for pplacer parallel operation.
    except KeyError:
        splits = 1
        
    try:
        query = command_args['query']
        query = query.split('.')
        
        ## Finally add some handling of extensions!
        
        if query[-1] in ['fasta', 'fna', 'fa']:
            query = query[0:-1]
        query = '.'.join(query)
        
    except KeyError:
        pass
    
## Count the system cpus, and divide by number of splits to determine cpus
## that should be used by cmalign

system_cpus = multiprocessing.cpu_count()
cmalign_cpus = str(int(system_cpus / splits))

## Figure out an appropriate number of cores for building trees.

if hyperthreading == True:
    physical_cpus = system_cpus/2
    
if physical_cpus <= 24:
    raxml_cpus = 1
else:
    raxml_cpus = int(physical_cpus / 24)
    
## Regardless of the number of cores available, should never use more than
## 6 per tree.

if raxml_cpus > 6:
    raxml_cpus = 6
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
ref_dir_domain = ref_dir + domain + '/' # Complete path to the domain subdirectory within the reference directory.
cm = paprica_path + 'models/' + domain + '_ssu.cm' # Domain specific covariance model used by infernal.

#%% Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])
    
#%% Define function to clean record names.  Note that bad_character is also used
## by make_tax, to insure that sequence names match between taxonomy database
## and tree.
    
from Bio import SeqIO

def clean_name(prefix, bad_character):
    
    bad_character = re.compile('[\[\]\|\\=-@!%,;\(\):\'\"\s]')
        
    with open(prefix + '.clean.fasta', 'w') as fasta_out:
        for record in SeqIO.parse(prefix + '.fasta', 'fasta'): 
            record.name = re.sub(bad_character, '_', str(record.description))
            record.id = record.name
            record.description = ''
            SeqIO.write(record, fasta_out, 'fasta')
            
#%% Define a function to create a unique fasta file from the input.
            
def make_unique(query):
    
    seq_count = {}
    seq_names = {}
    name_seq = {}
    
    clean_name(query)
    
    for record in SeqIO.parse(query + '.clean.fasta', 'fasta'):
        name = str(record.id)
        seq = str(record.seq)
        
        if seq in list(seq_count.keys()):
            seq_count[seq] = seq_count[seq] + 1
            temp_name_list = seq_names[seq]
            temp_name_list.append(name)
            seq_names[seq] = temp_name_list
        else:
            seq_count[seq] = 1
            seq_names[seq] = [name]
            
    with open(query + '.clean.unique.fasta', 'w') as fasta_out, open(query + '.clean.unique.count', 'w') as count_out:
        print('rep_name' + ',' + 'abundance', file = count_out)
        for seq in list(seq_names.keys()):
            print('>' + seq_names[seq][0], file=fasta_out)
            print(seq, file=fasta_out)
            
            print(seq_names[seq][0] + ',' + str(seq_count[seq]), file=count_out)
            name_seq[seq_names[seq][0]] = seq
            
    return seq_count, seq_names, name_seq
            
#%% Define function to split fasta file to run pplacer in parallel.  This greatly
## improves the speed of paprica_run.sh, but at the cost of memory overhead.
## Users need to be cautious of memory limits when considering the number of
## splits to make.
            
def split_fasta(file_in, nsplits):
    
    splits = []            
    tseqs = len(re.findall('>', open(file_in + '.fasta', 'r').read()))
    
    nseqs = tseqs / nsplits
    
    seq_i = 0
    file_n = 1
    
    file_out = open(file_in + '.temp' + str(file_n) + '.fasta', 'w')
    
    for record in SeqIO.parse(file_in + '.fasta', 'fasta'):
        seq_i = seq_i + 1
        
        if seq_i <= nseqs:
            SeqIO.write(record, file_out, 'fasta')
        elif seq_i > nseqs:
            splits.append(file_in + '.temp' + str(file_n))
            file_out.close()
            file_n = file_n + 1
            file_out = open(file_in + '.temp' + str(file_n) + '.fasta', 'w')
            SeqIO.write(record, file_out, 'fasta')
            seq_i = 0
        
    splits.append(file_in + '.temp' + str(file_n))
    file_out.close()
    
    return(splits)
    
#%% Define a function to split a combined query and reference (such as one
## combined by a call to esl-alimerge).
    
def split_query_ref(query, ref, combined):
    
    ref_out_name = ref.split('/')[-1]
    ref_out_name = '.'.join(ref_out_name.split('.')[0:-1])
    ref_out_name = ref_out_name + '.newlength.fasta'
    
    query_out_name = '.'.join(query.split('.')[0:-1])
    query_out_name = query_out_name + '.newlength.fasta'
    
    query_names = []
    ref_names = []
    
    for record in SeqIO.parse(query, 'stockholm'):
        query_names.append(record.id)
        
    for record in SeqIO.parse(ref, 'stockholm'):
        ref_names.append(record.id)
        
    query_names = set(query_names)
    ref_names = set(ref_names)
    
    with open(ref_out_name, 'w') as ref_out, open(query_out_name, 'w') as query_out:
        for record in SeqIO.parse(combined, 'stockholm'):
            if record.id in query_names:
                SeqIO.write(record, query_out, 'fasta')
            elif record.id in ref_names:
                SeqIO.write(record, ref_out, 'fasta')
            else:
                print('Error, record.id not in query or ref!')
                break

#%% Define function to execute phylogenetic placement on a query fasta, or split query fasta
    
def place(query, ref, ref_dir_domain, cm):
            
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + query + '.clean.unique.fasta', shell = True, executable = executable)
    degap.communicate()
    
    ## Conduct alignment with Infernal (cmalign) against the domain profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    
    if 'large' in list(command_args.keys()):
        align = subprocess.Popen('cmalign --cpu ' + cmalign_cpus + ' --mxsize 12000 --dna -o ' + query + '.clean.unique.align.sto --outformat Pfam ' + cm + ' ' + query + '.clean.unique.fasta', shell = True, executable = executable)
        align.communicate()
    else:
        align = subprocess.Popen('cmalign --cpu ' + cmalign_cpus + ' --dna -o ' + query + '.clean.unique.align.sto --outformat Pfam ' + cm + ' ' + query + '.clean.unique.fasta', shell = True, executable = executable)
        align.communicate()    
    
    combine = subprocess.Popen('esl-alimerge --outformat Pfam --dna \
    -o ' + query + '.' + ref + '.clean.unique.align.sto \
    ' + query + '.clean.unique.align.sto \
    ' + ref_dir_domain + ref + '.clean.align.sto', shell = True, executable = executable)
    combine.communicate()
    
    split_query_ref(query + '.clean.unique.align.sto', ref_dir_domain + ref + '.clean.align.sto', query + '.' + ref + '.clean.unique.align.sto')
    
    ## Create a specific working directory in case multiple instances are
    ## being run at same time.
    
    os.system('mkdir ' + query)
    
    epa_ng = subprocess.Popen('epa-ng -q ' + query + '.clean.unique.align.newlength.fasta \
    --model ' + ref_dir_domain + ref + '.final.bestModel \
    -s ' + ref + '.clean.align.newlength.fasta \
    -t ' + ref_dir_domain + ref + '.final.bestTree \
    -w ' + query, shell = True, executable = executable)
    epa_ng.communicate()
    
    ## Move final files to primary working directory, and delete subdirectory.
    
    os.system('mv ' + query + '/epa_result.jplace ' + query + '.' + ref + '.clean.unique.align.jplace')
    os.system('rm -r ' + query)
    
#%% Define function to generate csv file of placements and fat tree 
    
def json_to_csv(query, ref):
    with open(query + '.' + ref + '.clean.unique.align.jplace', 'r') as jfile:
        data = json.load(jfile)
        colnames = data['fields']
  
    placements = pd.DataFrame(columns = colnames)

    for placement in data['placements']:
        placements.loc[placement['n'][0]] = placement['p'][0]
        
    placements.to_csv(query + '.' + ref + '.clean.unique.align.csv')
    
def guppy(query, ref):
    
    guppy1 = subprocess.Popen('guppy to_csv --point-mass -o ' + query + '.' + ref + '.clean.unique.align.csv ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy1.communicate()
    
    guppy2 = subprocess.Popen('guppy fat --node-numbers --point-mass -o ' + query + '.' + ref + '.clean.unique.align.phyloxml ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy2.communicate()
    
    guppy3 = subprocess.Popen('guppy edpl --csv -o ' + query + '.' + ref + '.clean.unique.align.edpl.csv ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy3.communicate()
    
#%% Define function to execute a single raxml-ng run
    
def run_raxml(i, prefix, msa, raxml_cpus):
    
    ## Create trees with parsimony starts.
    
    subprocess.run('raxml-ng \
                   --redo \
                   --search \
                   --msa ' + msa + ' \
                   --tree pars{1} \
                   --prefix ' + prefix + str(i) + ' \
                   --seed ' + str(i) + ' \
                   --threads ' + str(raxml_cpus), shell = True, executable = '/bin/bash')
    
#%% Replacement for make_tax, based on ete3.  Build a table of lineages for each
    ## reference.
                   
def classify_ref():

    from ete3 import NCBITaxa
    
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database() 
    summary_complete = pd.read_csv(ref_dir_domain + 'genome_data.csv.gz', header = 0, index_col = 0)
    
    ref_lineage = pd.DataFrame()
    
    for strain in summary_complete['taxid']:
        lineage = ncbi.get_lineage(strain)
        lineage_ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        
        for taxid in lineage_ranks.keys():
            rank = lineage_ranks[taxid]
            ref_lineage.loc[strain, rank] = names[taxid]
            
    ref_lineage.to_csv(ref_dir_domain + 'edge_lineages.csv')
    
#%% Define function to run infernal and return a fasta format file.
    
def infernal_align(prefix, cm):
    
    ## Make sure input is unaligned.  Input file must end with ".clean.fasta".
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + prefix + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    infernal_commands = 'cmalign --dna -o ' + prefix + '.clean.align.sto --outformat Pfam ' + cm + ' ' + prefix + '.clean.fasta'    
    
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + prefix + '.clean.align.sto ' + prefix + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()  

#%% Execute main program.

if 'query' not in list(command_args.keys()):
    
    ## If the query flag is not given this is taken as instruction to build
    ## the reference package. 
    
    classify_ref()
    
    clean_name(ref_dir_domain + ref)
    
    ## If cmalign returns an error complaining about the size of the DP matrix, there are probably
    ## reference 16S or 18S sequences that fall outside the bounds of the covariance model.
    ## These need to be removed (add to the appropriate bad_genomes list in paprica-make_ref.py)
    ## as they will result in a malformed tree.
    
    infernal_align(ref_dir_domain + ref, cm)
    
    ## Check the input alignment, this identifies the optimal number of threads,
    ## eliminates duplicates, and converts to binary format. Alignment file with
    ## no duplicates is identified by *.reduced.phy.
    
    rm = subprocess.call('rm -f ' + ref_dir_domain + '*raxml*', shell = True, executable = executable)
    
    raxml0 = subprocess.call('raxml-ng \
    --parse \
    --msa ' + ref_dir_domain + ref + '.clean.align.fasta \
    --model GTR+G \
    --prefix ' + ref_dir_domain + ref, shell = True, executable = executable)
    
    ## Use the coarse-grained parallelization scheme outlined in the RAxML
    ## manual to create 24 trees from parsimony starts.  This takes a long
    ## time for the domain Bacteria.
    
    if __name__ == '__main__':  
        Parallel(n_jobs = 24, verbose = 5)(delayed(run_raxml)
        (i, ref_dir_domain + ref, ref_dir_domain + ref + '.raxml.rba', raxml_cpus) for i in range(1, 25))
        
    ## Parse the log files and find the tree with the highest maximum likelihood.
        
    final_lls = pd.Series(dtype = 'float64')
        
    for f in os.listdir(ref_dir_domain):
        if f.endswith('raxml.log'):
            with open(ref_dir_domain + f, 'r') as raxml_log:
                for line in raxml_log:
                    if line.startswith('Final LogLikelihood'):
                        line = line.rstrip()
                        line = line.split()
                        final_ll = float(line[-1])
                        final_lls[f] = final_ll
                        
    besttree = final_lls.idxmax()
    besttree_base = '.'.join(besttree.split('.')[0:-1])
    besttree_name = besttree_base + '.bestTree'
    besttree_log = besttree_base + '.log'
    besttree_model = besttree_base + '.bestModel'
    
    ## now delete the other trees and rename the winning
    
    os.rename(ref_dir_domain + besttree_name, ref_dir_domain + ref + '.final.bestTree')
    os.rename(ref_dir_domain + besttree_log, ref_dir_domain + ref + '.final.log')
    os.rename(ref_dir_domain + besttree_model, ref_dir_domain + ref + '.final.bestModel')
    
    subprocess.run('rm ' + ref_dir_domain + '*raxml*', shell = True, executable = '/bin/bash')
    
    ## Create a file with the date/time of database creation.

    current_time = datetime.datetime.now().isoformat()
    n_aseqs = len(re.findall('>', open(ref_dir_domain + '/' + ref + '.fasta', 'r').read()))
    
    with open(ref_dir_domain + '/' + ref + '.database_info.txt', 'w') as database_info:
        print('ref tree built at:', current_time, file = database_info)
        print('nseqs in reference alignment:', n_aseqs, file = database_info) 
            
else:
        
    clear_wspace = subprocess.call('rm -f ' + cwd + query + '.' + ref + '*', shell = True, executable = executable)
    clear_wspace = subprocess.call('rm -f ' + cwd + query + '.sub*', shell = True, executable = executable)
    
    ## Select a random subset of reads, if this option is specified.  This is useful for
    ## normalizing the number of sampled reads across different samples.    
    
    if 'n' in list(command_args.keys()):
        
        nseqs = command_args['n']
        tseqs = len(re.findall('>', open(cwd + query + '.fasta', 'r').read()))
        nseqs_get = random.sample(list(range(1, tseqs)), int(nseqs))
        
        seq_i = 0

        with open(cwd + query + '.sub.fasta', 'w') as fasta_sub:
            for record in SeqIO.parse(cwd + query + '.fasta', 'fasta'):
                seq_i = seq_i + 1
                if seq_i in nseqs_get:
                    SeqIO.write(record, fasta_sub, 'fasta')
        
        ## All downstream operations now need to take place on the subsampled
        ## query, easiest way to do this is to just change the query variable.
            
        query = query + '.sub'
        
    ## Make a unique fasta file.
    
    unique_count, unique_names, unique_seq_names = make_unique(cwd + query)
                            
    ## Create splits, if splits > 1
    
    if splits > 1:        
        split_list = split_fasta(cwd + query, splits)
            
        if __name__ == '__main__':  
            Parallel(n_jobs = splits, verbose = 5)(delayed(place)
            (split_query, ref, ref_dir_domain, cm) for split_query in split_list)
            
        guppy_merge = subprocess.Popen('guppy merge ' + cwd + query + '*' + domain + '*' + '.jplace -o ' + cwd + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
        guppy_merge.communicate()
        guppy(cwd + query, ref)
        
        ## Note that the below merge command sloppily adds the reference sequences
        ## for each of the split files, those this doesn't have any negative impact
        ## other than space.
        
        merge = subprocess.Popen('cat ' + cwd + query + '.temp*.' + ref + '.clean.unique.align.fasta > ' + cwd + query + '.' + ref + '.clean.unique.align.fasta', shell = True, executable = executable)
        merge.communicate()
                
        cleanup = subprocess.Popen('rm -f ' + cwd + query + '.temp*', shell = True, executable = executable)
        cleanup.communicate()
        
    else:
        place(cwd + query, ref, ref_dir_domain, cm)
        guppy(cwd + query, ref)
        
    ## Add the abundance of each unique read to the guppy csv file and the edpl value.
    
    guppy_csv = pd.read_csv(cwd + query + '.' + ref + '.clean.unique.align.csv', header = 0, index_col = 'name')
    edpl_csv = pd.read_csv(cwd + query + '.' + ref + '.clean.unique.align.edpl.csv', header = None, names = ['edpl'])
    unique_abund = pd.read_csv(cwd + query + '.clean.unique.count', header = 0, index_col = 'rep_name')
    
    guppy_csv.loc[edpl_csv.index, 'edpl'] = edpl_csv.edpl
    guppy_csv.loc[unique_abund.index, 'abund'] = unique_abund.abundance
    
    ## Add the seq of each unique read to the guppy csv file
    
    for name in guppy_csv.index:
        row_seq = unique_seq_names[name]
        guppy_csv.loc[name, 'seq'] = str(row_seq)
        
    guppy_csv.to_csv(cwd + query + '.' + ref + '.clean.unique.align.csv')
