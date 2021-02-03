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
    -ref_dir: The directory containing the paprica database.
    -splits: The number of files to split the query fasta into to facilitate
        parallel analysis with pplacer.

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
import re
import subprocess
import sys
import os
import multiprocessing
import datetime
import random
import pandas as pd
import json
from Bio import Phylo, SeqIO, Seq
from io import StringIO
import tempfile
import shutil
import numpy as np

executable = '/bin/bash' # shell for executing commands
random.seed(1)

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
    domain = 'bacteria'
try:
    ref = command_args['ref']  # The name of the reference package being used.
except KeyError:
    if domain == 'eukarya':
        ref = 'combined_18S.' + domain +'.tax'
    else:
        ref = 'combined_16S.' + domain +'.tax'
try:
    query = command_args['query']
except KeyError:
    pass

## If sys.argv == 1, you are probably running inside Python in testing mode.
## Provide some default values to make this possibe.  If > 1, parse command
## line for optional variables, providing reasonable defaults if not present.
## No default is currently provided for query.
    
if len(sys.argv) == 1:
    query = 'test.' + domain
    command_args['query'] = query

## Figure out an appropriate number of cores for building trees.  RAxML
## manual suggests using only physical cores.

system_cpus = multiprocessing.cpu_count()

if hyperthreading == True:
    physical_cpus = system_cpus/2
    
if physical_cpus <= 24:
    raxml_cpus = 1
else:
    raxml_cpus = int(physical_cpus / 24)
    
#!!! override for testing
raxml_cpus = 3
    
## Regardless of the number of cores available, should never use more than
## 6 per tree.

if raxml_cpus > 6:
    raxml_cpus = 6
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'
    
ref_dir_domain = ref_dir + domain + '/' # Complete path to the domain subdirectory within the reference directory.
cm16S = paprica_path + 'models/' + domain + '_ssu.cm' # Domain specific covariance model used by infernal.
cm23S = paprica_path + 'models/' + domain + '_lsu.cm'

#%% Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])
    
#%% Define function to clean record names.  Note that bad_character is also used
## by make_tax, to insure that sequence names match between taxonomy database
## and tree.
    
def clean_name(query):
    
    bad_character = re.compile('[\[\]\\=-@!%,;\(\):\'\"\s]')
    
    basename = query.split('.')[0:-1]
    basename = '.'.join(basename)
    name_out = basename + '.clean.fasta'
        
    with open(name_out, 'w') as fasta_out:
        for record in SeqIO.parse(query, 'fasta'): 
            record.name = re.sub(bad_character, '_', str(record.description))
            record.id = record.name
            record.description = ''
            SeqIO.write(record, fasta_out, 'fasta')
            
#%% Define a function to create a unique fasta file from the input.
            
def make_unique(query):
    
    seq_count = {}
    seq_names = {}
    name_seq = {}
    
    basename = query.split('.')[0:-1]
    basename = '.'.join(basename)
    fasta_out_name = basename + '.unique.fasta'
    count_out_name = basename + '.unique.count'
    
    for record in SeqIO.parse(query, 'fasta'):
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
            
    with open(fasta_out_name, 'w') as fasta_out, open(count_out_name, 'w') as count_out:
        print('rep_name' + ',' + 'abundance', file = count_out)
        for seq in list(seq_names.keys()):
            print('>' + seq_names[seq][0], file = fasta_out)
            print(seq, file=fasta_out)
            
            print(seq_names[seq][0] + ',' + str(seq_count[seq]), file=count_out)
            name_seq[seq_names[seq][0]] = seq
            
    return seq_count, seq_names, name_seq
            
#%% Define function to split fasta file to run pplacer in parallel.  This greatly
## improves the speed of paprica_run.sh, but at the cost of memory overhead.
## Users need to be cautious of memory limits when considering the number of
## splits to make.
    
## Currently this function is not used.
            
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
    
def split_query_ref(query_alignment, ref_alignment, combined, temp_dir, part_file = None):
    
    ## We assume that reads are placing to the first partition. Right now this
    ## only works for a 2-partition alignment, but could be easily expanded to
    ## accomodate additional partitions.

    if part_file != None:  
        i = 0
        with open(part_file) as part:
            for line in part:
                if i == 1:
                    line = line.rstrip()
                    line = line.split(',')[1]
                    needed_gaps = int(line.split('-')[1])
                i += 1
    else:
        needed_gaps = 0
        
    ## turn off the extensioning for now
        
    needed_gaps = 0
        
    query_basename = query_alignment.split('/')[-1]
    query_basename = query_basename.split('.')[0:-1]
    query_basename = '.'.join(query_basename)
    
    ref_basename = ref_alignment.split('/')[-1]
    ref_basename = ref_basename.split('.')[0:-1]
    ref_basename = '.'.join(ref_basename)
    query_out = query_basename + '.newlength.fasta'
    ref_out = ref_basename + '.newlength.fasta'
    
    query_names = []
    ref_names = []
    
    for record in SeqIO.parse(query_alignment, 'stockholm'):
        query_names.append(record.id)
        
    for record in SeqIO.parse(ref_alignment, 'stockholm'):
        ref_names.append(record.id)
        
    query_names = set(query_names)
    ref_names = set(ref_names)
    
    ## Outputs saved in query specific results directory.
    
    with open(temp_dir + ref_out, 'w') as ref_out, open( temp_dir + query_out, 'w') as query_out:
        for record in SeqIO.parse(combined, 'stockholm'):
            
            record.letter_annotations = {}
            
            if record.id in query_names:
                temp_seq = str(record.seq) + '-' * needed_gaps
                record.seq = Seq.Seq(temp_seq)
                SeqIO.write(record, query_out, 'fasta')
            elif record.id in ref_names:
                temp_seq = str(record.seq) + '-' * needed_gaps
                record.seq = Seq.Seq(temp_seq)
                SeqIO.write(record, ref_out, 'fasta')
            else:
                print('Error, record.id not in query or ref!')
                break    
    
#%% Define function to execute a single raxml-ng run
    
def run_raxml(prefix, msa):
    
    ## Create trees with parsimony starts.
    
    subprocess.run('raxml-ng \
                --redo \
                --search \
                --msa ' + msa + ' \
                --prefix ' + prefix + ' \
                --workers 18 \
                --tree pars{9},rand{9} \
                --threads 36', shell = True, executable = '/bin/bash')
    
#%% Replacement for make_tax, based on ete3.  Build a table of lineages for each
    ## reference.
                   
def classify_ref(ref_dir_domain):

    from ete3 import NCBITaxa
    import numpy as np
    
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database() 
    summary_complete = pd.read_csv(ref_dir_domain + 'genome_data.csv.gz', header = 0, index_col = 0)
    
    #!!! This statement will be unnecessary after checking that summary-complete
    ## and final_genes are reconciled at end of paprica-make_ref
    
    summary_complete.dropna(subset = ['tax_name'], inplace = True)
    
    ref_lineage = pd.DataFrame()
    
    for strain in summary_complete['taxid']:
        
        #!!! In theory a try clause shouldn't be necessary here, as there really
        ## isn't a reason this should fail.  Previously ete3 was failing on
        ## the most recent version of the NCBI taxonomy, so this was necessary.
        ## Solved by inserting NA values for taxids with missing lineage, however,
        ## this will probably break the determination of consensus lineages
        ## in build_core_genomes.
        
        try:
            lineage = ncbi.get_lineage(strain)
            lineage_ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)
            
            for taxid in lineage_ranks.keys():
                rank = lineage_ranks[taxid]
                ref_lineage.loc[strain, rank] = names[taxid]
        
        except ValueError:
            ref_lineage.loc[strain, 'no rank'] = np.nan
            
    ref_lineage = ref_lineage.fillna(np.nan)
    ref_lineage.to_csv(ref_dir_domain + 'edge_lineages.csv')
    
    return(ref_lineage, summary_complete)
    
#%% Define function to run infernal and return a fasta format file.
    
def reference_align(prefix, cm):
    
    ## Make sure input is unaligned.  Input file must end with ".clean.fasta".
    
    degap = subprocess.Popen('seqmagick mogrify --ungap ' + prefix + '.clean.fasta', shell = True, executable = executable)
    degap.communicate()
    
    infernal_commands = 'cmalign --dna -o ' + prefix + '.clean.align.sto --outformat Pfam ' + cm + ' ' + prefix + '.clean.fasta'    
    
    infernal = subprocess.Popen(infernal_commands, shell = True, executable = executable)
    infernal.communicate()
    
    convert = subprocess.Popen('seqmagick convert ' + prefix + '.clean.align.sto ' + prefix + '.clean.align.fasta', shell = True, executable = executable)
    convert.communicate()
    
def concatenate_alignments(align1, align2, partition_file, file_out):
    
    alignments = pd.DataFrame()
    
    for record in SeqIO.parse(align1, 'fasta'):
        alignments.loc[str(record.id), 'align1'] = str(record.seq)
    for record in SeqIO.parse(align2, 'fasta'):
        alignments.loc[str(record.id), 'align2'] = str(record.seq)
        
    alignments.dropna(inplace = True)
    part1_len = len(alignments.iloc[0]['align1'])
    part2_len = len(alignments.iloc[0]['align2'])
    
    with open(partition_file, 'w') as part_out:
        print('GTR+G+FO,', '16S=1-' + str(part1_len), file = part_out)
        print('GTR+G+FO,', '23S=' + str(part1_len + 1) + '-' + str(part1_len + part2_len), file = part_out)

    with open(file_out, 'w') as concat_out:                    
        for index, row in alignments.iterrows():
            print('>' + index, file = concat_out)            
            print(row['align1'] + row['align2'], file = concat_out) 
    
def trim_alignment(prefix):
    
    bases = re.compile('[AGCTagct]')
    first_positions = []
    last_positions = []
        
    for record in SeqIO.parse(prefix + '.clean.align.fasta', 'fasta'):
        
        seq = str(record.seq)
        i_front = re.search(bases, seq).span()[0]
        i_back = re.search(bases, seq[::-1]).span()[0]
        
        first_positions.append(i_front)
        last_positions.append(i_back)
        
    first_position = max(first_positions)
    last_position = max(last_positions)
    
    ## If all sequences spanned the full length of the alignment last_position
    ## will be nonsensical.  In this case the full sequence should be retained.
    
    if last_position == 0:
        last_position = 1
    
    i = 0
    
    with open(prefix + '.clean.align.trimmed.fasta', 'w') as trimmed_fasta:
        for record in SeqIO.parse(prefix + '.clean.align.fasta', 'fasta'):
            record.seq = record.seq[first_position:(-1 * last_position)]
            SeqIO.write(record, trimmed_fasta, 'fasta')
            i = i + 1
            
    return(i)

def parse_alignment(prefix, part):
    
    if part  != None:
        subprocess.call('raxml-ng \
        --parse \
        --msa ' + prefix + '.clean.align.trimmed.fasta \
        --model ' + part + ' \
        --prefix ' + prefix, shell = True, executable = executable)
        
    else:
        subprocess.call('raxml-ng \
        --parse \
        --msa ' + prefix + '.clean.align.trimmed.fasta \
        --model GTR+G+FO \
        --prefix ' + prefix, shell = True, executable = executable)        
    
def best_tree(prefix, partition):
    
    ## If some trees failed to converge and the run was forced to quit,
    ## parse the log file to select highest log-likelihood among those
    ## trees that did converge.
    
    with open(prefix + '.raxml.log', 'r') as log_in:
        parse_tree = True
        ll = -1 * np.Inf
        for line in log_in:
            line = line.rstrip()
            line = line.split()
            if 'logLikelihood:' in line:
                if float(line[-1]) > ll:
                    ll = float(line[-1])
                    besttree = int(line[6].replace('#', '').replace(',', ''))
            if 'Final' in line:
                parse_tree = False
                
    if parse_tree == True:
        i = 0
        with open(prefix + '.raxml.mlTrees', 'r') as trees_in, open(prefix + '.raxml.bestTree', 'w') as tree_out:
            for tree in trees_in:
                i = i + 1
                if i == besttree:
                    print(tree, file = tree_out)
                    break
                
    ## rename the final trees
    
    os.rename(prefix + '.raxml.bestTree', prefix + '.final.bestTree')
    os.rename(prefix + '.raxml.log', prefix + '.final.log')
    
    ## if a multi-partition alignment, need to create a model
    ## that is useable with epa-ng.
    
    besttree_model = prefix + '.raxml.bestModel'
    
    if partition == False:
        os.rename(besttree_model, prefix + '.final.bestModel')
    else:
        with open(besttree_model, 'r') as best_model, open(prefix + '.final.bestModel', 'w') as model_out:
            
            ## Assumes desired partition is in line 1
            
            for line in best_model:
                line = re.sub('\+BU\{\d*\.\d*\}', '', line)
                break
            print(line, file = model_out)
            
            shutil.copyfile(besttree_model, prefix + '.original.bestModel')
            
#%% Define function to check that all references are present in the final tree.
## This is a diagnostic function and is hopefully no longer needed!
            
def check_tree(alignment_16S, final_tree):
    
    ## Parse the final raxml-ng tree and get all the terminal names.
    
    tree_ids = []
    tree = Phylo.read(final_tree, 'newick') 
    for clade in tree.get_terminals():
        tree_ids.append(clade.name)
    tree_ids = set(tree_ids)
    
    ## Parse the 16S alignment.
    
    good_16S = []
    
    for record in SeqIO.parse(alignment_16S, 'stockholm'):
        if record.id in tree_ids:
            good_16S.append(record)
            
    if len(good_16S) != len(tree_ids):
        print('The length of the alignment does not match the tree!')
        stop_here()
                    
#%% Define function to execute phylogenetic placement on a query fasta, or split query fasta
    
def query_align(query, ref, combined, cm):
    
    ## combined is the desired name of the combined sto file
    ## ref is the sto file for the reference sequences
    ## query is the query fasta
    
    basename = query.split('.')[0:-1]
    basename = '.'.join(basename)
        
    subprocess.call('seqmagick mogrify --ungap ' + query, shell = True, executable = executable)
    
    ## Conduct alignment with Infernal (cmalign) against the domain profile
    ## obtained from the Rfam website at http://rfam.xfam.org/family/RF00177/cm.
    
    subprocess.call('cmalign --dna -o ' + basename + '.align.sto --outformat Pfam ' + cm + ' ' + query, shell = True, executable = executable)   

    subprocess.call('esl-alimerge --outformat Pfam --dna \
    -o ' + combined + ' \
    ' + basename + '.align.sto \
    ' + ref, shell = True, executable = executable)
    
def place(query_alignment, ref_alignment, model, tree, temp_dir):
    
    ## model: *bestModel
    ## tree: *final.bestTree
    ## file_out: name for final jplace file, traditionally query + '.' + ref + '.clean.unique.align.jplace'
        
    subprocess.call('epa-ng --redo -q ' + temp_dir + query_alignment + ' \
    --model ' + model + ' \
    -s ' + temp_dir + ref_alignment + ' \
    -t ' + tree + ' \
    -w ' + temp_dir, shell = True, executable = executable)

def gappa(jplace, cwd):
    
    basename = jplace.split('.')[0:-1]
    basename = '.'.join(basename)
    
    subprocess.call('gappa examine edpl \
                    --allow-file-overwriting \
                    --out-dir ' + cwd + ' \
                    --file-prefix ' + basename + '.edpl \
                    --jplace-path ' + cwd + jplace, shell = True, executable = executable)
                    
    subprocess.call('gappa examine heat-tree \
                    --allow-file-overwriting \
                    --out-dir ' + cwd + ' \
                    --write-phyloxml-tree \
                    --tree-file-prefix ' + basename + ' \
                    --jplace-path ' + cwd + jplace, shell = True, executable = executable)

#%% Define function to convert json to csv.
    
def json_to_csv(jplace, subtree_ref):
    with open(jplace, 'r') as jfile:
        data = json.load(jfile)
        colnames = data['fields']
        
    ## Parse the tree and exchange curly brackets for brackets.
        
    tree = data['tree']
    tree = re.sub('{', '[', tree)
    tree = re.sub('}', ']', tree)
    tree = Phylo.read(StringIO(tree), 'newick') 
        
    ref_seqs = pd.DataFrame(columns = ['ref_name'])
    subtree = subtree_ref
    
    for clade in tree.get_terminals():
    	ref_name = clade.name
    	ref_seqs.loc[int(clade.comment)] = [ref_name]
        
    ## Parse placements.
  
    placements = pd.DataFrame(columns = colnames)

    for placement in data['placements']:
        if len(placement['p']) == 1:
            placements.loc[placement['n'][0]] = placement['p'][0]
        else:
            
            ## Pick the edge that has the highest ML value.
            
            edge_ml = []
            for pp in placement['p']:
                edge_ml.append(pp[1])
            
            edge_i = edge_ml.index(max(edge_ml))
            placements.loc[placement['n'][0]] = placement['p'][edge_i]
            
    placements = placements.astype({'edge_num': 'int32'})
    placements['subtree'] = subtree
    
    ## Parsing subtree to get a reasonable prefix for the global edge number
    ## isn't super robust, but should work so long as the subtree name doesn't
    ## have "." in it.  The global edge number must exactly reflect the indices
    ## of the database files, as this is how query reads are connected to 
    ## that information.
    
    subtree_suffix = subtree.split('.')[-1]
    placements['global_edge_num'] = subtree_suffix + '_' + placements.edge_num.astype('str') 
    
    ## Add ref seqs.
        
    ref_names = ref_seqs.reindex(placements.edge_num).ref_name
    ref_names.index = placements.index       
    placements['ref_name'] = ref_names   
        
    return(placements)

#%% Define function to calculate map_ratio.
    
def get_map_ratio(query_alignment, ref_alignment, placements):
    
    ## query and ref alignment should be *.newlength.fasta
    
    ref_names = set()
    ref_seqs = {}
    
    for ref_name in placements.ref_name.unique():
        ref_names.add(ref_name)
        
    ## Find all the reference seqs that had hits.
        
    for record in SeqIO.parse(ref_alignment, 'fasta'):
        if record.id in ref_names:
            ref_seqs[record.id] = str(record.seq)
            
    for record in SeqIO.parse(query_alignment, 'fasta'):
        ref_name = placements.loc[record.id, 'ref_name']
        if pd.notnull(ref_name):
            ref_seq = ref_seqs[ref_name]
            
            n_total = 0
            n_match = 0
            
            for i,j in enumerate(str(record.seq)):
                if j not in ['-', '.']:
                    n_total = n_total + 1
                    if j.capitalize() == ref_seq[i].capitalize():
                        n_match = n_match + 1
                        
            unaligned = str(record.seq).replace('-', '')
            unaligned = unaligned.replace('.', '').upper()
                        
            map_ratio = round(float(n_match) / n_total, 2)
            placements.loc[record.id, 'map_ratio'] = map_ratio
            placements.loc[record.id, 'map_id'] = n_match
            placements.loc[record.id, 'seq'] = unaligned
            placements.loc[record.id, 'origin'] = query_alignment
        
    return(placements)

## Define function to map nodes to subtrees in the master tree.
    
def map_master_tree(jplace):
    
    master_map = {}
    
    with open(jplace, 'r') as jfile:
        data = json.load(jfile)
        
    tree = data['tree']
    tree = re.sub('{', '[', tree)
    tree = re.sub('}', ']', tree)
    tree = Phylo.read(StringIO(tree), 'newick') 
        
    for clade in tree.get_terminals():
        clade_name = '_'.join(clade.name.split('_')[:-1])
        master_map[clade.comment] = clade_name
            
    for clade in tree.get_nonterminals():
        temp_subtrees = set()
        
        for terminal in clade.get_terminals():
            terminal_name = '_'.join(terminal.name.split('_')[:-1])
            temp_subtrees.add(terminal_name)
        
        if len(temp_subtrees) == 1:
            
            ## Then all terminal edges in clade belong to same subtree.
            
            clade_name = list(temp_subtrees)[0]
            master_map[clade.comment] = clade_name

    return(master_map)
            
#%% Define euk refs as needed.
    
euk_reps = {'Ochrophyta':'AY229897.1.1739_U',
            'Dinoflagellata':'FJ549370.1.1796_U',
            'Chlorophyta__Streptophyta':'AF506698.1.1787_U'}

#%% Execute main program.
    
if domain != 'eukarya':
    prefix_16S = ref_dir_domain + 'combined_16S.' + domain + '.tax'
    prefix_23S = ref_dir_domain + 'combined_23S.' + domain + '.tax'
    prefix_combined = ref_dir_domain + 'combined_16S.23S.' + domain + '.tax'
else:
    
    #!!! For eukarya prefix_16S and prefix_combined are the same.  To make this
    ## cleaner could change _16S to _ssu in future version, also for cm16S.
    
    prefix_16S = ref_dir_domain + 'combined_18S.eukarya.tax'
    prefix_combined = ref_dir_domain + 'combined_18S.eukarya.tax'
        
top_level_names = {'bacteria':'phylum_reps',
           'archaea':'archaea',
           'eukarya':'div_reps'}
    
## If the query flag is not given this is taken as instructions to build
## the reference package.
    
if 'query' not in list(command_args.keys()):
    
    if domain == 'archaea':
        
        ## For domain archaea it is not necessary (yet) to subdivide below the
        ## level of domain. However, to keep things consistent for downstream
        ## scripts, it shoud be named as though it's a subtree.
        
        ref_lineages, genome_data = classify_ref(ref_dir_domain)
       
        ## Make a copy of the reference file which will serve as the "subtree"
        ## file for domain archaea.
        
        shutil.copyfile(prefix_16S + '.fasta', prefix_16S + '.archaea.fasta')
        shutil.copyfile(prefix_23S + '.fasta', prefix_23S + '.archaea.fasta')
        
        fasta = 'archaea'
        
        clean_name(prefix_16S + '.' + fasta + '.fasta')
        clean_name(prefix_23S + '.' + fasta + '.fasta')
        
        reference_align(prefix_16S + '.' + fasta, cm16S)
        reference_align(prefix_23S + '.' + fasta, cm23S)
                
        nseqs16S = trim_alignment(prefix_16S + '.' + fasta)
        nseqs23S = trim_alignment(prefix_23S + '.' + fasta)
        
        if nseqs16S != nseqs23S:
            print('You have different numbers of 16S and 23S rRNA genes!  Check \
                  your input files!')
                  
            quit()
        
        subprocess.call('rm -f ' + ref_dir_domain + '*' + fasta + '*raxml*', shell = True, executable = executable)
        
        concatenate_alignments(prefix_16S + '.' + fasta + '.clean.align.trimmed.fasta',
                               prefix_23S + '.' + fasta + '.clean.align.trimmed.fasta',
                               prefix_combined + '.' + fasta + '.part',
                               prefix_combined + '.' + fasta + '.clean.align.trimmed.fasta')

        parse_alignment(prefix_combined + '.' + fasta, part = prefix_combined + '.' + fasta + '.part')
        run_raxml(prefix_combined + '.' + fasta, prefix_combined + '.' + fasta + '.raxml.rba')        
        best_tree(prefix_combined + '.' + fasta, partition = True)
        
        ## Clean up by moving the files that will be needed for paprica-run.sh
        ## to a dedicated directory.
            
        shutil.rmtree(ref_dir_domain + fasta, ignore_errors = True)
        os.makedirs(ref_dir_domain + fasta)
        
        for f in ['combined_16S.23S.archaea.tax.' + fasta + '.final.bestModel',
                  'combined_16S.23S.archaea.tax.' + fasta + '.part',
                  'combined_23S.archaea.tax.' + fasta + '.clean.align.sto',
                  'combined_16S.23S.archaea.tax.' + fasta + '.final.bestTree',
                  'combined_16S.archaea.tax.' + fasta + '.clean.align.sto',
                  'combined_16S.23S.archaea.tax.' + fasta + '.final.log']:
            
            shutil.copy(ref_dir_domain + f, ref_dir_domain + fasta + '/' + f)
        
        genome_data['subtree'] = domain        
        fastas = {fasta:nseqs16S}       
        n_aseqs = genome_data.shape[0]
        genome_data.to_csv(ref_dir_domain + 'genome_data.csv.gz')
        
    elif domain == 'bacteria':
    
        ## The reference taxonomy is used to subdivide the reference into phyla.
        ## Not all phyla have sufficient sequences to build a tree, however, 
        ## at this level there are several phyla that already have too many
        ## taxa per alignment. This level is a good compromise.
        
        ## ref_lineages is indexed by taxid, so need to map phylum to
        ## reference strains using taxid.
        
        ref_lineages, genome_data = classify_ref(ref_dir_domain)
        
        ## If phyla should be grouped together indicate that here.  Iteration
        ## is now over analysis_group rather than phylum.
        
        ref_lineages['analysis_group'] = ref_lineages.phylum
        ref_lineages.loc[ref_lineages.phylum == 'Firmicutes', 'analysis_group'] = 'Firmicutes_Tenericutes'
        ref_lineages.loc[ref_lineages.phylum == 'Tenericutes', 'analysis_group'] = 'Firmicutes_Tenericutes'
        
        fastas = {'phylum_reps':0}
                
        with open(prefix_16S + '.phylum_reps.fasta', 'w') as reps_out_16S, open(prefix_23S + '.phylum_reps.fasta', 'w') as reps_out_23S:  
            i  = 0
            
            for phylum in ref_lineages.analysis_group.unique():
                
                if pd.notnull(phylum):
                    print(phylum)
                    j = 0
                    phylum_taxids = set(ref_lineages.index[ref_lineages.analysis_group == phylum])
                    phylum_seqids = list(genome_data.loc[genome_data['taxid'].isin(phylum_taxids)].tax_name)
                    phylum_lineages = ref_lineages.loc[ref_lineages['analysis_group'] == phylum]
                    
                    rep_seqs = []
                    
                    ## If less than 30 reference sequences in the phylum use
                    ## all of them on the master tree.
                    
                    if len(phylum_seqids) < 30:
                        rep_seqs = phylum_seqids
                        
                    ## If more than 30 reference sequences select one from each
                    ## genus.  Right now reference sequences are not considered
                    ## if they don't have a genus.
                    
                    else:
                        for genus in phylum_lineages.genus.unique():
                            if pd.notnull(genus):
                                genus_taxid = phylum_lineages.loc[phylum_lineages['genus'] == genus].index[0]
                                genus_rep = genome_data.loc[genome_data['taxid'] == genus_taxid].tax_name[0]
                                rep_seqs.append(genus_rep)
                                
                    ## Convert phylum_seqids to set for faster searching.  Do
                    ## not do this to rep_seqs because you need to be able to
                    ## pull indices for specific values.
                    
                    phylum_seqids = set(phylum_seqids)
                                    
                    ## Now iterate across combined_16S.bacteria.tax.fasta and create phylum level
                    ## fastas and fasta of representatives. It's inefficient to iterate 
                    ## over combined_16S.bacteria.tax.fasta for each phylum, but this works fine
                    ## for now. The alternative would involve a lot of opening/closing
                    ## files, or having one file open for each phylum.
                    
                    phylum = re.sub(' ', '_', phylum)
                    genome_data.loc[genome_data['taxid'].isin(phylum_taxids), 'subtree'] = phylum
                    
                    with open(prefix_16S + '.' + phylum + '.fasta', 'w') as fasta_out_16S:
                        for record in SeqIO.parse(prefix_16S + '.fasta', 'fasta'):
                            if record.description in phylum_seqids:
                                j += 1
                                SeqIO.write(record, fasta_out_16S, 'fasta')
                            if record.description in rep_seqs:
                                record.id = phylum + '_' + str(rep_seqs.index(record.description))
                                record.description = ''
                                i += 1
                                SeqIO.write(record, reps_out_16S, 'fasta')
                                
                    with open(prefix_23S + '.' + phylum + '.fasta', 'w') as fasta_out_23S:
                        for record in SeqIO.parse(prefix_23S + '.fasta', 'fasta'):
                            if record.description in phylum_seqids:
                                SeqIO.write(record, fasta_out_23S, 'fasta')
                            if record.description in rep_seqs:
                                record.id = phylum + '_' + str(rep_seqs.index(record.description))
                                record.description = ''
                                SeqIO.write(record, reps_out_23S, 'fasta')
                                                            
                    fastas[phylum] = j  
                    
        fastas['phylum_reps'] = i
        
        ## Sorting fastas.keys just to get Proteobacteria out of the #2 slot
        ## which makes troubleshooting easier.
                                                        
        for fasta in sorted(fastas.keys()):

            nseqs = fastas[fasta]
            clean_name(prefix_16S + '.' + fasta + '.fasta')
            clean_name(prefix_23S + '.' + fasta + '.fasta')
            
            reference_align(prefix_16S + '.' + fasta, cm16S)
            reference_align(prefix_23S + '.' + fasta, cm23S)
            
            nseqs16S = trim_alignment(prefix_16S + '.' + fasta)
            nseqs23S = trim_alignment(prefix_23S + '.' + fasta)
            
            if nseqs16S != nseqs23S:
                print('You have different numbers of 16S and 23S rRNA genes!  Check \
                      your input files!')
                      
                quit()
        
            subprocess.call('rm -f ' + ref_dir_domain + '*' + fasta + '*raxml*', shell = True, executable = executable)
            
            concatenate_alignments(prefix_16S + '.' + fasta + '.clean.align.trimmed.fasta',
                                   prefix_23S + '.' + fasta + '.clean.align.trimmed.fasta',
                                   prefix_combined + '.' + fasta + '.part',
                                   prefix_combined + '.' + fasta + '.clean.align.trimmed.fasta')
        
            if nseqs > 3:
                
#            if fasta == 'phylum_reps':
     
            ## build trees
            
                parse_alignment(prefix_combined + '.' + fasta, part = prefix_combined + '.' + fasta + '.part')                   
                run_raxml(prefix_combined + '.' + fasta, prefix_combined + '.' + fasta + '.raxml.rba')
                best_tree(prefix_combined + '.' + fasta, partition = True)
                check_tree(prefix_16S + '.' + fasta + '.clean.align.sto', prefix_combined + '.' + fasta + '.final.bestTree')
                
                ## Clean up by moving the files that will be needed for paprica-run.sh
                ## to a dedicated directory.
                    
                shutil.rmtree(ref_dir_domain + fasta, ignore_errors = True)
                os.makedirs(ref_dir_domain + fasta)
                
                for f in ['combined_16S.23S.bacteria.tax.' + fasta + '.final.bestModel',
                          'combined_16S.23S.bacteria.tax.' + fasta + '.part',
                          'combined_16S.23S.bacteria.tax.' + fasta + '.final.bestTree',
                          'combined_16S.bacteria.tax.' + fasta + '.clean.align.sto',
                          'combined_16S.23S.bacteria.tax.' + fasta + '.final.log']:
                    
                    try:
                        shutil.copy(ref_dir_domain + f, ref_dir_domain + fasta + '/' + f)
                        
                    ## If not enough sequences were present to build a tree, files will not
                    ## exist.
                        
                    except FileNotFoundError:
                        continue
                
        n_aseqs = genome_data.shape[0]
        genome_data.to_csv(ref_dir_domain + 'genome_data.csv.gz')
                
    elif domain == 'eukarya':
        
        ## The PR2 database is used for the eukarya reference tree. The
        ## merged file found on the PR2 Github repository has carriage returns
        ## in some fields. These need to be purged before the merged file
        ## can be used here.
        
        pr2 = pd.read_csv(ref_dir_domain + 'pr2_version_4.12.0_merged_nocarriage.tsv', sep = '\t', index_col = 0)
        pr2 = pr2[pr2.supergroup != 'Opisthokonta']
        pr2 = pr2[pr2.sequence_length > 1700]
        pr2 = pr2[pr2.sequence_length < 1900] 
        pr2 = pr2[pr2.label != 'G']
        pr2 = pr2[pr2.gene == '18S_rRNA']
        pr2 = pr2[pr2.organelle == 'nucleus']
        pr2 = pr2[pr2.chimera != 1]
        pr2 = pr2[pr2.ambiguities == 0]               
        pr2.drop_duplicates(subset = ['sequence'], inplace = True)
        pr2 = pr2[pr2.reference_sequence == 1]
        pr2.dropna(subset = ['pr2_accession'], inplace = True)
        
        ## Selecting different taxonomic levels to include. Taxonomic levels
        ## can be grouped where necessary using "__".  Ideally this section
        ## would be re-written to follow the "analysis_group" ontology of
        ## bacteria.
        
        ## If Dinophyceae and Syndiniales are analyzed seperately Lingulodinium 
        ## ends up in Syndiniales, switching to Dinoflagellata
        
        pr2_units = {'Chlorophyta__Streptophyta':'division',
                     'Dinoflagellata':'division',
                     'Oxyrrhea':'class',
                     'Noctilucophyceae':'class',
                     'Phaeophyceae':'class',
                     'Chrysophyceae':'class',
                     'Pelagophyceae':'class',
                     'Bacillariophyta':'class',
                     'Bolidophyceae':'class',
                     'Synurophyceae':'class',
                     'Ochrophyta_X':'class',
                     'Raphidophyceae':'class',
                     'Chrysomerophyceae':'class',
                     'Eustigmatophyceae':'class',
                     'Xanthophyceae':'class',
                     'Dictyochophyceae':'class',
                     'Pinguiophyceae':'class',
                     'Aurearenophyceae':'class',
                     'MOCH-2':'class',
                     'Phaeothamniophyceae':'class',
                     'Picophagea':'class',
                     'Ciliophora':'division',
                     'Discoba':'division',
                     'Cercozoa':'division',
                     'Apusomonadidae':'division',
                     'Conosa':'division',
                     'Rhodophyta':'division',
                     'Cryptophyta':'division',
                     'Apicomplexa':'division',
                     'Radiolaria':'division',
                     'Perkinsea':'division',
                     'Lobosa':'division',
                     'Opalozoa':'division',
                     'Picozoa':'division',
                     'Pseudofungi':'division',
                     'Sagenista':'division',
                     'Telonemia':'division',
                     'Hilomonadea':'division',
                     'Malawimonadidae':'division',
                     'Metamonada':'division',
                     'Katablepharidophyta':'division',
                     'Glaucophyta':'division',
                     'Centroheliozoa':'division',
                     'Hacrobia_X':'division',
                     'Alveolata_X':'division',
                     'Breviatea':'division',
                     'Haptophyta':'division'}

        ## create a fasta file for each target clade

        fastas = {'div_reps':0}

        with open(prefix_16S + '.div_reps.fasta', 'w') as rep_out:
            i = 0
            
            for tc_name in pr2_units.keys():
                
                tax_level = pr2_units[tc_name]
                tc = tc_name.split('__')
                j = 0
                
                temp = pr2[pr2[tax_level].isin(tc)]
                
                pr2.loc[pr2[tax_level].isin(tc), 'subtree'] = tc_name
                
                ## Identify 20 representatives for this clade, or if < 20 members
                ## in clade select a representative from each genus.
                
                if temp.shape[0] < 20:
                    rep_seqs = temp
                    
                else:
                    
                    ## Select a representative from each genus in the clade.
                    
                    rep_seqs = pd.DataFrame()
                        
                    for genus in temp.genus.unique():
                        temp_genus = temp.loc[temp.genus == genus]
                        rep_seqs = pd.concat([rep_seqs, temp_genus.iloc[random.sample(range(0, temp_genus.shape[0]), 1),:]])
                    
                ## If euk rep is pre-determined for this clade, identify
                ## here and add to representatives if needed.
                    
                try:
                    rep_seq_id = euk_reps[tc_name]
                    
                    if rep_seq_id not in rep_seqs['pr2_accession']:
                        rep_seqs = pd.concat([rep_seqs, temp.loc[temp['pr2_accession'] == rep_seq_id]])

                except KeyError:
                    rep_seq_id = None
                
                for index, row in rep_seqs.iterrows():
                    print('>' + tc_name + '_' + str(i), file = rep_out)
                    print(row.sequence, file = rep_out)
                    i += 1
                
                with(open(prefix_16S + '.' + tc_name + '.fasta', 'w')) as fasta_out:
                    for index, row in temp.iterrows():
                        pr2.loc[index, 'tax_name'] = row.pr2_accession + '|' + row.species
                        print('>' + row.pr2_accession + '|' + row.species, file = fasta_out)
                        print(row.sequence, file = fasta_out)
                        j += 1
                        
                fastas[tc_name] = j
        fastas['div_reps'] = i
        
        for fasta in fastas.keys():
            clean_name(prefix_16S + '.' + fasta + '.fasta')
            reference_align(prefix_16S + '.' + fasta, cm16S)
            trim_alignment(prefix_16S + '.' + fasta)
            
            # subprocess.call('seqmagick mogrify --deduplicate-sequences ' + ref_dir_domain + fasta + '.fasta',
            #                 shell = True,
            #                 executable = executable)
            
            nseqs = len(re.findall('>', open(prefix_16S + '.' + fasta + '.fasta', 'r').read()))
            fastas[fasta] = nseqs
            
            subprocess.call('rm -f ' + ref_dir_domain + '*' + fasta + '*raxml*', shell = True, executable = executable)
            
            if nseqs > 3:
                
            ## This option is for optimizing specific trees.
            
            #if fasta in ['div_reps']:
     
            ## build trees
            
                parse_alignment(prefix_16S + '.' + fasta, None)   
                run_raxml(prefix_16S + '.' + fasta, prefix_16S + '.' + fasta + '.raxml.rba') 
                best_tree(prefix_16S + '.' + fasta, partition = False)
                
                ## Clean up by moving the files that will be needed for paprica-run.sh
                ## to a dedicated directory.
                    
                shutil.rmtree(ref_dir_domain + fasta, ignore_errors = True)
                os.makedirs(ref_dir_domain + fasta)
                
                for f in ['combined_18S.eukarya.tax.' + fasta + '.final.bestModel',
                          'combined_18S.eukarya.tax.' + fasta + '.clean.align.sto',
                          'combined_18S.eukarya.tax.' + fasta + '.final.bestTree',
                          'combined_18S.eukarya.tax.' + fasta + '.final.log']:
                    
                    try:
                        shutil.copy(ref_dir_domain + f, ref_dir_domain + fasta + '/' + f)
                        
                    ## If not enough sequences were present to build a tree, files will not
                    ## exist.
                        
                    except FileNotFoundError:
                        continue
                
        ## edge_lineages.csv and genome_data.csv need to be compatible with
        ## paprica-build_core_genomes.py.  The classify_ref function is not used
        ## for PR2 because not all entries have an NCBI taxid.
                
        ref_lineages = pr2[['kingdom', 'supergroup', 'division', 'class',
                            'order', 'family','genus', 'species', 'taxo_id']]
        ref_lineages.index = ref_lineages.taxo_id
        ref_lineages = ref_lineages.drop_duplicates()
        ref_lineages.to_csv(ref_dir_domain + 'edge_lineages.csv')
        
        pr2.rename(columns={'taxo_id':'taxid'}, inplace = True)
        pr2.index = pr2.pr2_accession
        pr2.to_csv(ref_dir_domain + 'genome_data.csv.gz')
        
        n_aseqs = pr2.shape[0]
        
    ## Because you're no longer using test.fasta to create a single reference 
    ## tree there is no way to get all the edge names from a single run with
    ## test.fasta. Best way to deal with this is to run the reps fasta file
    ## against all of the reference trees.  This is done for all domains.
        
    ## For archaea the ref sequences have the same name when being placed
    ## against themselves. This is a problem because easel adds a prefix to the
    ## sequence names when combining alignments. It's necessary to create a
    ## junk fasta to use for this purpose.
        
    query = ref + '.single'
    
    if domain != 'eukarya':
        short_prefix_ssu = 'combined_16S.' + domain + '.tax'
        short_prefix_ssu_lsu = 'combined_16S.23S.' + domain + '.tax'
    else:
        short_prefix_ssu = 'combined_18S.eukarya.tax'
        short_prefix_ssu_lsu = 'combined_18S.eukarya.tax'
    
    i = 0
    with open(ref_dir_domain + query + '.clean.fasta', 'w') as seq_out:
        for record in SeqIO.parse(prefix_16S + '.' + top_level_names[domain] + '.clean.fasta', 'fasta'):
            i += 1
            if i == 1:
                record.id = 'test_' + str(record.id)
                record.description = ''
                SeqIO.write(record, seq_out, 'fasta')
        
    for sub in fastas.keys():
        if fastas[sub] > 3:
            
            ## For 16S.23S you would need to align the query to a 16S alignment
            ## only.
            
            temp_dir = tempfile.mkdtemp(dir = ref_dir_domain) + '/'
            
            query_align(ref_dir_domain + query + '.clean.fasta',
                    ref_dir_domain + sub + '/' + short_prefix_ssu + '.' + sub + '.clean.align.sto',
                    temp_dir + query + '.' + ref + '.' + sub + '.clean.align.sto',
                    cm16S)
            
            ## Eukarya is not a partitioned alignment.
            
            if domain == 'eukarya':
                part_file = None
            else:
                part_file = ref_dir_domain + sub + '/' + short_prefix_ssu_lsu + '.' + sub + '.part'
            
            split_query_ref(ref_dir_domain + query + '.clean.align.sto',
                  ref_dir_domain + sub + '/' + short_prefix_ssu + '.' + sub + '.clean.align.sto',
                  temp_dir + query + '.' + ref + '.' + sub + '.clean.align.sto',
                  temp_dir,
                  part_file)
            
            place(query + '.clean.align.newlength.fasta',
                  ref + '.' + sub + '.clean.align.newlength.fasta',
                  ref_dir_domain + sub + '/' + short_prefix_ssu_lsu + '.' + sub + '.final.bestModel',
                  ref_dir_domain + sub + '/' + short_prefix_ssu_lsu + '.' + sub + '.final.bestTree',
                  temp_dir)
            
            shutil.copyfile(temp_dir + 'epa_result.jplace',
                             ref_dir_domain + sub + '/' + short_prefix_ssu_lsu + '.' + sub + '.jplace')
            
            shutil.rmtree(temp_dir)
        
    ## For all domains, create a file with the date/time of database creation.
    #!!! n_aseqs value not correct for all domains

    current_time = datetime.datetime.now().isoformat()
    
    with open(prefix_16S + '.database_info.txt', 'w') as database_info:
        print('ref tree built at:', current_time, file = database_info)
        print('nseqs in reference alignment:', n_aseqs, file = database_info) 
        
        for sub_ref in fastas.keys():
            print('*' + sub_ref + '\t' + str(fastas[sub_ref]), file = database_info)
            
else:
    
    ## If the query flag is given you are placing reads on the existing reference tree.
    ## Two features have been deprecated as they are no longer particularly relevant:
    ## the -n subsampling option and -splits.  It's recommended that subsampling
    ## be performed on the final output files in R or whatever software you're
    ## using for statistical analysis.  -splits does not offer a significant
    ## performance gain on unique input reads.
    
    ## Parse the database file to get the list of reference clades for which
    ## there should be subtrees.
    
    available_trees = set()
    
    with open(ref_dir_domain + ref + '.database_info.txt', 'r') as database_info:
        for line in database_info:
            line = line.rstrip()
            if line.startswith('*'):
                line = line.split('\t')
                if int(line[1]) > 3:
                    subtree = line[0]
                    subtree = subtree.strip('*')
                    available_trees.add(subtree)
                    
    ## Chose an appropriate name for top level reference tree based on
    ## domain.  phylum_ref is really only an appropriate variable name
    ## for domain Bacteria, but it does no harm.
    
    phylum_ref = top_level_names[domain]
    
    ## Since creating separate directories for each subtree the original
    ## prefix_16S and prefix_23S no longer point to correct locations.  So
    ## need to define a new short prefix variable here.  When used with
    ## ref_dir_domain and the name of the subtree it will point the correct
    ## path and prefix.
     
    if domain != 'eukarya':
        short_prefix_ssu = 'combined_16S.' + domain + '.tax'
        short_prefix_ssu_lsu = 'combined_16S.23S.' + domain + '.tax'
    else:
        short_prefix_ssu = 'combined_18S.eukarya.tax'
        short_prefix_ssu_lsu = 'combined_18S.eukarya.tax'
            
    clear_wspace = subprocess.call('rm -f ' + query + '.' + phylum_ref + '*',
                                   shell = True,
                                   executable = executable)
    
    ## reasonable name for results directory
    
    temp_dir = cwd + query
    temp_dir = temp_dir + '/'
    
    ## remove the previous results directory, if it exists
    
    shutil.rmtree(temp_dir, ignore_errors = True)
    
    ## create the new results directory
    
    os.makedirs(temp_dir, exist_ok = True)
    
    clean_name(cwd + query + '.fasta')
    make_unique(cwd + query + '.clean.fasta')

    query_align(cwd + query + '.clean.unique.fasta',
        ref_dir_domain + phylum_ref + '/' + short_prefix_ssu + '.' + phylum_ref + '.clean.align.sto',
        temp_dir + query + '.' + ref + '.' + phylum_ref + '.clean.unique.align.sto',
        cm16S)
    
    ## Eukarya is not a partitioned alignment.
            
    if domain == 'eukarya':
        part_file = None
    else:
        part_file = ref_dir_domain + phylum_ref + '/' + short_prefix_ssu_lsu + '.' + phylum_ref + '.part'

    split_query_ref(cwd + query + '.clean.unique.align.sto',
          ref_dir_domain + phylum_ref + '/' + short_prefix_ssu + '.' + phylum_ref + '.clean.align.sto',
          temp_dir + query + '.' + ref + '.' + phylum_ref + '.clean.unique.align.sto',
          temp_dir,
          part_file)

    place(query + '.clean.unique.align.newlength.fasta',
          ref + '.' + phylum_ref + '.clean.align.newlength.fasta',
          ref_dir_domain + phylum_ref + '/' + short_prefix_ssu_lsu + '.' + phylum_ref + '.final.bestModel',
          ref_dir_domain + phylum_ref + '/' + short_prefix_ssu_lsu + '.' + phylum_ref + '.final.bestTree',
          temp_dir)
    
    os.rename(temp_dir + 'epa_result.jplace', temp_dir + query + '.' + phylum_ref + '.jplace')
    
    placements = json_to_csv(temp_dir + query + '.' + phylum_ref + '.jplace', phylum_ref)
        
    placements = get_map_ratio(temp_dir + query + '.clean.unique.align.newlength.fasta',
        temp_dir + ref + '.' + phylum_ref + '.clean.align.newlength.fasta',
        placements)
                    
    gappa(query + '.' + phylum_ref + '.jplace', temp_dir)  
    edpl = pd.read_csv(temp_dir + query + '.' + phylum_ref + '.edpllist.csv', index_col = 1)
    placements = pd.concat([placements, edpl['EDPL']], axis = 1, sort = False)
    
    ## Add the count data.  This should also add all unique reads that didn't
    ## place to the original reference tree.
    
    if domain == 'archaea':
        combined_subtrees = placements
    
    else:
        
        ## If the domain is either bacteria or eukarya then placement to
        ## subtrees is necessary.
        
        ## First parse the master tree, mapping nodes to subtrees.
    
        master_map = map_master_tree(temp_dir + query + '.' + phylum_ref + '.jplace')
        
        ## Then create a fasta file of query reads for each subtree.
                        
        subtrees = set(master_map.get(edge) for edge in placements.edge_num.astype('str'))
        
        for subtree in subtrees:
            if pd.notnull(subtree):
                
                with open(cwd + subtree + '_' + query + '.clean.unique.fasta', 'w') as subtree_out:
                    for record in SeqIO.parse(query + '.clean.unique.fasta', 'fasta'):
                        try:
                            edge_num = placements.loc[record.id, 'global_edge_num']
                            edge_num = edge_num.split('_')[-1]
                            edge_subtree = master_map[edge_num]
                            
                            if subtree == edge_subtree:
                                SeqIO.write(record, subtree_out, 'fasta')
                                
                        except KeyError:
                            
                            ## This happens when placements.ref_name == nan,
                            ## indicating internal placement on master tree.
                            ## No action to be taken here, as those reads will
                            ## be captured with the subtree == nan iteration.
                            
                            continue
            else:
                with open(cwd + 'nosubtree_' + query + '.clean.unique.fasta', 'w') as subtree_out:
                    for record in SeqIO.parse(query + '.clean.unique.fasta', 'fasta'):
                        try:
                            edge_num = placements.loc[record.id, 'global_edge_num']
                            edge_num = edge_num.split('_')[-1]
                            edge_subtree = master_map[edge_num]
                            
                        except KeyError:
                            SeqIO.write(record, subtree_out, 'fasta') 
                            
        ## Then carry out the subtree analysis.
                            
        combined_subtrees = pd.DataFrame()
                        
        for subtree in subtrees:
            if subtree in available_trees:
                if pd.notnull(subtree):
                    
                    query_align(cwd + subtree + '_' + query + '.clean.unique.fasta',
                            ref_dir_domain + subtree + '/' + short_prefix_ssu + '.' + subtree + '.clean.align.sto',
                            temp_dir + query + '.' + ref + '.' + subtree + '.clean.unique.align.sto',
                            cm16S)
    
                    ## Eukarya is not a partitioned alignment.
                            
                    if domain == 'eukarya':
                        part_file = None
                    else:
                        part_file = ref_dir_domain + subtree + '/' + short_prefix_ssu_lsu + '.' + subtree + '.part'
                
                    split_query_ref(cwd + query + '.clean.unique.align.sto',
                          ref_dir_domain + subtree + '/' + short_prefix_ssu + '.' + subtree + '.clean.align.sto',
                          temp_dir + query + '.' + ref + '.' + subtree + '.clean.unique.align.sto',
                          temp_dir,
                          part_file)
                
                    place(query + '.clean.unique.align.newlength.fasta',
                          ref + '.' + subtree + '.clean.align.newlength.fasta',
                          ref_dir_domain + subtree + '/' + short_prefix_ssu_lsu + '.' + subtree + '.final.bestModel',
                          ref_dir_domain + subtree + '/' + short_prefix_ssu_lsu + '.' + subtree + '.final.bestTree',
                          temp_dir)
                    
                    os.rename(temp_dir + 'epa_result.jplace', temp_dir + query + '.' + subtree + '.jplace')
                    
                    subtree_csv = json_to_csv(temp_dir + query + '.' + subtree + '.jplace', subtree)
                        
                    subtree_csv = get_map_ratio(temp_dir + query + '.clean.unique.align.newlength.fasta',
                        temp_dir + ref + '.' + subtree + '.clean.align.newlength.fasta',
                        subtree_csv)
                    
                    gappa(query + '.' + subtree + '.jplace', temp_dir)  
                    edpl = pd.read_csv(temp_dir + query + '.' + subtree + '.edpllist.csv', index_col = 1)
                    subtree_csv = pd.concat([subtree_csv, edpl['EDPL']], axis = 1, sort = False)
                    subtree_csv['subtree'] = subtree
                    
                    combined_subtrees = pd.concat([combined_subtrees, subtree_csv])
                    subtree_csv.to_csv(temp_dir + query + '.' + ref + '.' + subtree + '.placements.csv') 
                    
        ## This is a pretty harsh purge of the intermediate files, but it works.
                    
        subprocess.call('rm -rf ' + cwd + '*_' + query + '*', shell = True, executable = executable)
        subprocess.call('rm -rf ' + cwd + 'nosubtree*' + query + '*', shell = True, executable = executable)
                    
        ## Identify all unique reads that were assigned to a clade that didn't
        ## have > 3 reference sequences, meaning that there is not reference tree.
                    
        no_subtrees = placements.index.difference(combined_subtrees.index)
        combined_subtrees = pd.concat([combined_subtrees, placements.loc[no_subtrees]], sort = False)               
        
    ## Add the count data.  This should also add all unique reads that didn't
    ## place to the original reference tree.  At this point you have three
    ## types of reads. 1) reads that did not place to the refernce tree, 2) reads
    ## that placed to the refernce tree but not a subtree, 3) reads that placed
    ## to the subtree.  All should be accounted for in the final placements.csv
    ## file.
            
    count = pd.read_csv(cwd + query + '.clean.unique.count', index_col = 0)
    combined_subtrees = pd.concat([combined_subtrees, count], axis = 1, sort = False)
    
    ## Then you need to add the sequences for those reads that did not place
    ## to the original reference tree.
    
    try:
        needed_seqs = list(combined_subtrees.index[pd.isnull(combined_subtrees.seq)])
    except AttributeError:
        
        ## This failure occurres on the test file, where no seq is ever defined
        ## for combined_subtrees.
        
        combined_subtrees.seq = np.nan
        needed_seqs = combined_subtrees.index[pd.isnull(combined_subtrees.seq)][0].tolist()
        
    needed_seqs = set(needed_seqs)
    
    for record in SeqIO.parse(query + '.clean.unique.fasta', 'fasta'):
        if record.id in needed_seqs:
            combined_subtrees.loc[record.id, 'seq'] = str(record.seq)
    
    combined_subtrees.index.rename('Pquery', inplace = True)
    combined_subtrees.to_csv(cwd + query + '.' + ref + '.placements.csv', index = 'Pquery')
