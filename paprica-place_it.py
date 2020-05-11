#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Sat Jan 03 08:59:11 2015

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

REQUIRES:
    Files:
        bacterial_ssu.cm
        
    Programs:
        infernal
        taxtastic
        seqmagick
        pplacer
        raxmlHPC-PTHREADS-AVX2
        
    Python modules:
        Bio
        joblib

CALL AS:
    paprica-place_it.py -domain [domain] -query [query] -ref [ref] -splits [splits] -n [nseqs] for analysis or
    paprica-place_it.py -domain [domain] -ref [ref] -cpus [ncpus] to generate a reference package.  

    Note that [ref] or [query] includes the entire file name without extension
    (which must be .fasta).
    
    It is not recommended to use more than 8 cpus for tree building (i.e. -cpus 8).
    See the RAXML manual for guidance on this.
    
OPTIONS:
    -cpus:  The number of cpus for RAxML to use.
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
    domain = 'bacteria'
try:
    ref = command_args['ref']  # The name of the reference package being used.
except KeyError:
    ref = 'combined_16S.bacteria.tax'

## If sys.argv == 1, you are probably running inside Python in testing mode.
## Provide some default values to make this possibe.  If > 1, parse command
## line for optional variables, providing reasonable defaults if not present.
## No default is currently provided for query.
    
if len(sys.argv) == 1:
    cpus = '8'
    splits = 1
    query = 'test.bacteria'
    command_args['query'] = query

else:
    
    try:
        cpus = str(command_args['cpus']) # The number of cpus for RAxML to use.
    except KeyError:
        cpus = '1'
    
    try:    
        splits = int(command_args['splits']) # The number of splits to make for pplacer parallel operation.
    except KeyError:
        splits = 1
        
    try:
        query = command_args['query']
        query = query.split('.')
        
        ## Finally add some handling of extenstions!
        
        if query[-1] in ['fasta', 'fna', 'fa']:
            query = query[0:-1]
        query = '.'.join(query)
        
    except KeyError:
        pass
    
## Count the system cpus, and divide by number of splits to determine cpus
## that should be used by cmalign

system_cpus = multiprocessing.cpu_count()
cmalign_cpus = str(int(system_cpus / splits))
    
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
    
#%% Define function to clean record names.  Not that bad_character is also used
## by make_tax, to insure that sequence names match between taxonomy database
## and tree.
    
bad_character = '[\[\]\|\\=-@!%,;\(\):\'\"\s]'

from Bio import SeqIO

def clean_name(file_name, bad_character):
    
    bad_character = re.compile(bad_character)
        
    with open(file_name + '.clean.fasta', 'w') as fasta_out:
        for record in SeqIO.parse(file_name + '.fasta', 'fasta'): 
            record.name = re.sub(bad_character, '_', str(record.description))
            record.id = record.name
            record.description = ''
            SeqIO.write(record, fasta_out, 'fasta')
            
#%% Define a function to create a unique fasta file from the input.
            
def make_unique(query):
    
    seq_count = {}
    seq_names = {}
    name_seq = {}
    
    clean_name(query, bad_character)
    
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
        align = subprocess.Popen('cmalign --cpu ' + cmalign_cpus + ' --mxsize 4000 --dna -o ' + query + '.clean.unique.align.sto --outformat Pfam ' + cm + ' ' + query + '.clean.unique.fasta', shell = True, executable = executable)
        align.communicate()
    else:
        align = subprocess.Popen('cmalign --cpu ' + cmalign_cpus + ' --dna -o ' + query + '.clean.unique.align.sto --outformat Pfam ' + cm + ' ' + query + '.clean.unique.fasta', shell = True, executable = executable)
        align.communicate()    
    
    combine = subprocess.Popen('esl-alimerge --outformat Pfam --dna \
    -o ' + query + '.' + ref + '.clean.unique.align.sto \
    ' + query + '.clean.unique.align.sto \
    ' + ref_dir_domain + ref + '.refpkg/' + ref + '.clean.align.sto', shell = True, executable = executable)
    combine.communicate()
    
    split_query_ref(query + '.clean.unique.align.sto', ref_dir_domain + ref + '.refpkg/' + ref + '.clean.align.sto', query + '.' + ref + '.clean.unique.align.sto')
    
    #os.system('rm -f epa_info.log') ## can remove this
    os.system('mkdir ' + query)
    
    epa_ng = subprocess.Popen('epa-ng -q ' + query + '.clean.unique.align.newlength.fasta \
    -s ' + ref + '.clean.align.newlength.fasta \
    -t ' + ref_dir_domain + ref + '.refpkg/RAxML_result.recalc.root.ref.tre \
    -w ' + query + ' \
    --model ' + ref_dir_domain + ref + '.refpkg/RAxML_info.recalc.root.ref.tre', shell = True, executable = executable)
    epa_ng.communicate()
    
    os.system('mv ' + query + '/epa_result.jplace ' + query + '.' + ref + '.clean.unique.align.jplace')
    os.system('rm -r ' + query)
    
#%% Define function to generate csv file of placements and fat tree    
    
def guppy(query, ref):
    
    guppy1 = subprocess.Popen('guppy to_csv --point-mass -o ' + query + '.' + ref + '.clean.unique.align.csv ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy1.communicate()
    
    guppy2 = subprocess.Popen('guppy fat --node-numbers --point-mass -o ' + query + '.' + ref + '.clean.unique.align.phyloxml ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy2.communicate()
    
    guppy3 = subprocess.Popen('guppy edpl --csv -o ' + query + '.' + ref + '.clean.unique.align.edpl.csv ' + query + '.' + ref + '.clean.unique.align.jplace', shell = True, executable = executable)
    guppy3.communicate()
    
#%% Define a function to classify read placements.
    
def classify():
    
    ## Prep sqlite database for classification.  Would be great not to have to do
    ## this, but classification command requires sqlite output.
    
    prep_database = subprocess.Popen('rppr prep_db -c ' + ref_dir_domain + ref + '.refpkg \
    --sqlite ' + query + '.' + ref + '.clean.unique.align.db',
    shell = True,
    executable = executable)
    
    prep_database.communicate()
    
    ## Now execute guppy classification command.
    
    guppy_classify = subprocess.Popen('guppy classify \
    -c ' + ref_dir_domain + ref + '.refpkg \
    --sqlite ' + query + '.' + ref + '.clean.unique.align.db '+ query + '.' + ref + '.clean.unique.align.jplace',
    shell = True,
    executable = executable)
    
    guppy_classify.communicate()    
    
#%% Define function to generate taxonomic information for ref package.

def make_tax(bad_character):
        
    ## Output taxon_id to tax_ids.txt, and other information to seq_info.csv for taxit.
    ## See http://fhcrc.github.io/taxtastic/quickstart.html for more information.
    ## Note that taxon_id or taxid are floating numbers, where they should be integers.
    ## Because there are NaN values, must remove before conversion to integers.
        
    ## Blank tax_id would be better than removing whole line, but not clear
    ## how to do this for now.
        
    ## Create csv file holding all information required for seq_info.csv.

    seq_info = pd.DataFrame(columns = ['seqname', 'accession', 'tax_id', 'species_name', 'is_type'])
    
    summary_complete = pd.read_csv(ref_dir_domain + 'genome_data.csv.gz', header = 0, index_col = 0)
    
    ## Need filler taxid for entries that don't have one.
    
    taxid = {'eukarya':2759, 'bacteria':2, 'archaea':2157}
    
    if domain == 'eukarya':
        summary_complete.taxon_id = summary_complete.taxon_id.fillna(value = taxid[domain])
        summary_complete.taxon_id = summary_complete.taxon_id.astype(dtype = 'uint64')
        seq_info['tax_id'] = summary_complete['taxon_id']

    else:
        summary_complete.taxid = summary_complete.taxid.fillna(value = taxid[domain])
        summary_complete.taxid = summary_complete.taxid.astype(dtype = 'uint64')
        seq_info['tax_id'] = summary_complete['taxid']
        
    ## Writing out genome_data.csv here allows the placeholder taxids
    ## for draft genomes, or other genomes without taxid, to be used
    ## downstream.  A better solution would be to find and add taxids for
    ## draft genomes.
        
    summary_complete.to_csv(ref_dir_domain + 'genome_data.csv.gz')
        
    ## Sequence names must be cleaned exactly as in clean_name.  Drop any entries
    ## that do not have a seqname.
    
    seq_info['seqname'] = summary_complete.tax_name.str.replace(bad_character, '_')
    seq_info.dropna(subset = ['seqname'], inplace = True)
        
    ## Write out the seq_info.csv file.
        
    seq_info.to_csv(ref_dir_domain + 'seq_info.csv', index = False)
        
    ## Download NCBI taxonomy database and load into sqlite database.  In case
    ## the taxonomy database has been updated, delete the old one.  Rebuilding
    ## the database takes a bit, however.
        
    taxtable_1 = subprocess.Popen('rm -f ' + ref_dir + 'taxonomy.db;\
    rm -f ' + ref_dir + 'taxdmp.zip;\
    taxit new_database ' + ref_dir + 'taxonomy.db', shell = True, executable = executable)
    taxtable_1.communicate()
    
    ## Probably some of your taxids are old.  Update them.
    
    taxtable_2 = subprocess.Popen('taxit update_taxids \
    -o ' + ref_dir_domain + 'seq_info.updated.csv \
    ' + ref_dir_domain + 'seq_info.csv \
    ' + ref_dir + 'taxonomy.db', shell = True, executable = executable)
    taxtable_2.communicate()
    
    ## Generate the file tax_ids.txt based on the newly generated seq_info.updates.csv.
    
    seq_info = pd.read_csv(ref_dir_domain + 'seq_info.updated.csv', header = 0)
    seq_info.to_csv(ref_dir_domain + 'tax_ids.txt', columns = ['tax_id'], header = False, index = False)
    
    taxtable_3 = subprocess.Popen('taxit taxtable \
    ' + ref_dir + 'taxonomy.db \
    -f ' + ref_dir_domain + 'tax_ids.txt \
    -o ' + ref_dir_domain + 'taxa.csv', shell = True, executable = executable)
    taxtable_3.communicate()

#%% Execute main program.

if 'query' not in list(command_args.keys()):
    
    ## If the query flag is not given this is taken as instruction to build
    ## the reference package.  
    
    clean_name(ref_dir_domain + ref, bad_character)
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + ref_dir_domain + ref + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
        
    ## If cmalign returns an error complaining about the size of the DP matrix, there are probably
    ## reference 16S or 18S sequences that fall outside the bounds of the covariance model.
    ## These need to be removed (add to the appropriate bad_genomes list in paprica-make_ref.py)
    ## as they will result in a malformed tree.
    
    infernal_commands = 'cmalign --dna -o ' + ref_dir_domain + ref + '.clean.align.sto --outformat Pfam ' + cm + ' ' + ref_dir_domain + ref + '.clean.fasta'    
    
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + ref_dir_domain + ref + '.clean.align.sto ' + ref_dir_domain + ref + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()   
    
    # Construct an initial tree using the GTRGAMMA model and the user provided
    # number of cpus.
    
    rm = subprocess.call('rm -f ' + ref_dir_domain + '*ref.tre', shell = True, executable = executable)
    raxml1 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -m GTRGAMMA -s ' + ref_dir_domain + ref + '.clean.align.fasta -n ref.tre -f d -p 12345 -w ' + ref_dir_domain, shell = True, executable = executable)
    raxml1.communicate()
    
    ## Root the tree.
    
    raxml2 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T 2 -m GTRGAMMA -f I -t ' + ref_dir_domain + 'RAxML_bestTree.ref.tre -n root.ref.tre -w ' + ref_dir_domain, shell = True, executable = executable)  
    raxml2.communicate()
    
    ## Generate SH-like support values for the tree.
    
#    raxml3 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -m GTRGAMMA -f J -p 12345 -t ' + ref_dir_domain + 'RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s ' + ref_dir_domain + ref + '.clean.align.fasta -w ' + ref_dir_domain, shell = True, executable = executable)   
#    raxml3.communicate()
    
    ## Generate info file according the epa-ng instructions.
    
    raxml4 = subprocess.Popen('raxmlHPC-PTHREADS-AVX2 -T ' + cpus + ' -m GTRGAMMA -f e -t ' + ref_dir_domain + 'RAxML_rootedTree.root.ref.tre -n recalc.root.ref.tre -s ' + ref_dir_domain + ref + '.clean.align.fasta -w ' + ref_dir_domain, shell = True, executable = executable)   
    raxml4.communicate()
    
    ## Generate taxonomy information for the reference package.
    
    make_tax(bad_character)
     
    ## Generate the reference package using the rooted tree with SH-like support values and a log file.
    ## Will not overwrite existing reference package, so delete if present.
    
    rm = subprocess.call('rm -rf ' + ref_dir_domain + ref + '.refpkg', shell = True, executable = executable)
    
    taxit = subprocess.Popen('taxit create -l 16S_rRNA -P ' + ref_dir_domain + ref + '.refpkg \
    --aln-fasta ' + ref_dir_domain + ref + '.clean.align.fasta \
    --tree-stats ' + ref_dir_domain + 'RAxML_info.recalc.root.ref.tre \
    --tree-file ' + ref_dir_domain + 'RAxML_result.recalc.root.ref.tre \
    --aln-sto ' + ref_dir_domain + ref + '.clean.align.sto \
    --seq-info ' + ref_dir_domain + 'seq_info.updated.csv \
    --taxonomy ' + ref_dir_domain + 'taxa.csv \
    --no-reroot \
    ', shell = True, executable = executable)
    
    taxit.communicate()
    
    ## Create a file with the date/time of database creation.

    current_time = datetime.datetime.now().isoformat()
    n_aseqs = len(re.findall('>', open(ref_dir_domain + '/' + ref + '.fasta', 'r').read()))
    
    with open(ref_dir_domain + '/' + ref + '.database_info.txt', 'w') as database_info:
        print('ref tree built at:', current_time, file=database_info)
        print('nseqs in reference alignment:', n_aseqs, file=database_info) 
            
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
    
    ## Executing the classify function to add information to the taxonomy db,
    ## but this information is not currently being used downstream.
        
    classify()
