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
from Bio import Phylo, SeqIO
from io import StringIO
import tempfile
import shutil

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
    domain = 'eukarya'
try:
    ref = command_args['ref']  # The name of the reference package being used.
except KeyError:
    if domain == 'eukarya':
        ref = 'combined_18S.' + domain +'.tax'
    else:
        ref = 'combined_16S.' + domain +'.tax'

## If sys.argv == 1, you are probably running inside Python in testing mode.
## Provide some default values to make this possibe.  If > 1, parse command
## line for optional variables, providing reasonable defaults if not present.
## No default is currently provided for query.
    
if len(sys.argv) == 1:
    query = 'big_test.' + domain
#    command_args['query'] = query

## Figure out an appropriate number of cores for building trees.

system_cpus = multiprocessing.cpu_count()

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
    
def split_query_ref(query_alignment, ref_alignment, combined, temp_dir):
        
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
    
    ## Create a specific working directory in case multiple instances are
    ## being run at same time.
    
    with open(temp_dir + '/' + ref_out, 'w') as ref_out, open( temp_dir + '/' + query_out, 'w') as query_out:
        for record in SeqIO.parse(combined, 'stockholm'):
            if record.id in query_names:
                SeqIO.write(record, query_out, 'fasta')
            elif record.id in ref_names:
                SeqIO.write(record, ref_out, 'fasta')
            else:
                print('Error, record.id not in query or ref!')
                break    
    
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
                   
def classify_ref(ref_dir_domain):

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

def parse_alignment(prefix):
    subprocess.call('raxml-ng \
    --parse \
    --msa ' + prefix + '.clean.align.trimmed.fasta \
    --model GTR+G \
    --prefix ' + prefix, shell = True, executable = executable)
    
def best_tree(prefix, n_trees):
    
    ## Parse the log files and find the tree with the highest maximum likelihood.
        
    final_lls = pd.Series(dtype = 'float64')
        
    for i in range(1, n_trees + 1):
        with open(prefix + str(i) + '.raxml.log', 'r') as raxml_log:
            for line in raxml_log:
                if line.startswith('Final LogLikelihood'):
                    line = line.rstrip()
                    line = line.split()
                    final_ll = float(line[-1])
                    final_lls[prefix + str(i)] = final_ll
                        
    besttree = final_lls.idxmax()
    besttree_name = besttree + '.raxml.bestTree'
    besttree_log = besttree + '.raxml.log'
    besttree_model = besttree + '.raxml.bestModel'
    
    ## now delete the other trees and rename the winning
    
    os.rename(besttree_name, prefix + '.final.bestTree')
    os.rename(besttree_log, prefix + '.final.log')
    os.rename(besttree_model, prefix + '.final.bestModel')
    
    subprocess.run('rm ' + prefix + '*raxml*', shell = True, executable = '/bin/bash')

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
    
def place(query_alignment, ref_alignment, model, tree, file_out, temp_dir):
    
    ## model: *bestModel
    ## tree: *final.bestTree
    ## file_out: name for final jplace file, traditionally query + '.' + ref + '.clean.unique.align.jplace'
    
    subprocess.call('epa-ng --redo -q ' + temp_dir + '/' + query_alignment + ' \
    --model ' + model + ' \
    -s ' + temp_dir + '/' + ref_alignment + ' \
    -t ' + tree + ' \
    -w ' + temp_dir, shell = True, executable = executable)
    
    ## Move final files to primary working directory, and delete subdirectory.
    
    os.system('mv ' + temp_dir + '/epa_result.jplace ' + file_out)

def gappa(jplace, cwd):
    
    basename = jplace.split('.')[0:-1]
    basename = '.'.join(basename)
    
    subprocess.call('gappa examine edpl \
                    --allow-file-overwriting \
                    --file-prefix ' + basename + '.edpl. \
                    --jplace-path ' + cwd + jplace, shell = True, executable = executable)
                    
    subprocess.call('gappa examine heat-tree \
                    --allow-file-overwriting \
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
                        
            map_ratio = round(float(n_match) / n_total, 2)
            placements.loc[record.id, 'map_ratio'] = map_ratio
            placements.loc[record.id, 'map_id'] = n_match
        
    return(placements)
            
#%% Execute main program.
    
## If the query flag is not given this is taken as instruction to build
## the reference package.

if 'query' not in list(command_args.keys()):
    
    prefix = ref_dir_domain + ref
    
    ## For domain archaea it is not necessary (yet) to subdivide below the
    ## level of domain.
    
    if domain == 'archaea':
        
        ref_lineages, genome_data = classify_ref(ref_dir_domain)
        
        clean_name(prefix + '.fasta')
        reference_align(prefix, cm)
        nseqs = trim_alignment(prefix)
        
        subprocess.call('rm -f ' + prefix + '*raxml*', shell = True, executable = executable)
            
        parse_alignment(prefix)   
        Parallel(n_jobs = 12, verbose = 5)(delayed(run_raxml)
            (i, prefix, prefix + '.raxml.rba', 1) for i in range(1, 25))
        best_tree(prefix, 24)
        
        genome_data['subtree'] = domain        
        fastas = {ref:nseqs}       
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
        
        fastas = {ref + '.phylum_reps':0}
                
        with open(prefix + '.phylum_reps.fasta', 'w') as reps_out:  
            i  = 0
            
            for phylum in ref_lineages.phylum.unique():
                if pd.notnull(phylum):
                    j = 0
                    phylum_taxids = set(ref_lineages.index[ref_lineages.phylum == phylum])
                    phylum_seqids = genome_data.loc[genome_data['taxid'].isin(phylum_taxids)].tax_name
                    
                    rep_seq = phylum_seqids[random.randrange(len(phylum_seqids))]                
                    phylum_seqids = set(phylum_seqids)
                                    
                    ## Now iterate across combined_16S.fasta and create phylum level
                    ## fastas and fasta of representatives. It's inefficient to iterate 
                    ## over combined_16S.fasta for each phylum, but this works fine
                    ## for now. The alternative would involve a lot of opening/closing
                    ## files, or having one file open for each phylum.
                    
                    phylum = re.sub(' ', '_', phylum)
                    genome_data.loc[genome_data['taxid'].isin(phylum_taxids), 'subtree'] = phylum
                    
                    with open(prefix + '.' + phylum + '.fasta', 'w') as fasta_out:
                        for record in SeqIO.parse(prefix + '.fasta', 'fasta'):
                            if record.description in phylum_seqids:
                                j += 1
                                SeqIO.write(record, fasta_out, 'fasta')
                            if record.description == rep_seq:
                                record.id = phylum
                                record.description = ''
                                i += 1
                                SeqIO.write(record, reps_out, 'fasta')
                                                            
                    fastas[ref + '.' + phylum] = j                    
        fastas[ref + '.phylum_reps'] = i
                                                        
        for fasta in fastas.keys(): 
            nseqs = fastas[fasta]
            clean_name(ref_dir_domain + fasta + '.fasta')
            reference_align(ref_dir_domain + fasta, cm)
            trim_alignment(ref_dir_domain + fasta)
        
            subprocess.call('rm -f ' + fasta + '*raxml*', shell = True, executable = executable)
        
            if nseqs > 3:
     
            ## build trees
            
                parse_alignment(ref_dir_domain + fasta)   
                Parallel(n_jobs = 12, verbose = 5)(delayed(run_raxml)
                  (i, ref_dir_domain + fasta, ref_dir_domain + fasta + '.raxml.rba', 1) for i in range(1, 25))                
                best_tree(ref_dir_domain + fasta, 24)
                
        n_aseqs = genome_data.shape[0]
        genome_data.to_csv(ref_dir_domain + 'genome_data.csv.gz')
                
    elif domain == 'eukarya':
        
        ## The PR2 database is used for the eukarya reference tree. The
        ## merged file found on the PR2 Github repository has carriage returns
        ## in some fields. These need to be purged before the merged file
        ## can be used here.
        
        pr2 = pd.read_csv(ref_dir_domain + 'pr2_version_4.12.0_merged_nocarriage.tsv', sep = '\t')
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
        pr2['subtree'] = pr2.division

        ## create a fasta file for each target clade

        fastas = {ref + '.div_reps':0}

        with open(prefix + '.div_reps.fasta', 'w') as rep_out:
            i = 0
            
            for tc in pr2.division.unique():
                j = 0
                
                temp = pr2[pr2.division == tc]
                rep_seq_i = random.randrange(temp.shape[0])
                rep_seq = temp.iloc[rep_seq_i,:]
                
                print('>' + tc, file = rep_out)
                print(rep_seq.sequence, file = rep_out)
                i += 1
                
                with(open(prefix + '.' + tc + '.fasta', 'w')) as fasta_out:
                    for index, row in temp.iterrows():
                        pr2.loc[index, 'tax_name'] = row.pr2_accession + '|' + row.species
                        print('>' + row.pr2_accession + '|' + row.species, file = fasta_out)
                        print(row.sequence, file = fasta_out)
                        j += 1
                        
                fastas[ref + '.' + tc] = j
        fastas[ref + '.div_reps'] = i
        
        for fasta in fastas.keys():
            clean_name(ref_dir_domain + fasta + '.fasta')
            reference_align(ref_dir_domain + fasta, cm)
            trim_alignment(ref_dir_domain + fasta)
            
            # subprocess.call('seqmagick mogrify --deduplicate-sequences ' + ref_dir_domain + fasta + '.fasta',
            #                 shell = True,
            #                 executable = executable)
            
            nseqs = len(re.findall('>', open(ref_dir_domain + fasta + '.fasta', 'r').read()))
            fastas[fasta] = nseqs
            
            subprocess.call('rm -f ' + ref_dir_domain + fasta + '*raxml*', shell = True, executable = executable)
            
            if nseqs > 3:
     
            ## build trees
            
                parse_alignment(ref_dir_domain + fasta)   
                Parallel(n_jobs = 12, verbose = 5)(delayed(run_raxml)
                  (i, ref_dir_domain + fasta, ref_dir_domain + fasta + '.raxml.rba', 1) for i in range(1, 25))                
                best_tree(ref_dir_domain + fasta, 24)
                
        ## edge_lineages.csv and genome_data.csv need to be compatble with
        ## paprica-build_core_genomes.py.
                
        ref_lineages = pr2[['kingdom', 'supergroup', 'division', 'class', 'order', 'family','genus', 'species', 'taxo_id']]
        ref_lineages.index = ref_lineages.taxo_id
        ref_lineages.to_csv(ref_dir_domain + 'edge_lineages.csv')
        
        pr2.rename(columns={'taxo_id':'taxid'})
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
        
    top_level_names = {'bacteria':ref + '.phylum_reps',
                   'archaea':ref + '.single',
                   'eukarya':ref + '.div_reps'}
    
    query = top_level_names[domain]
    
    if domain == 'archaea':
        i = 0
        with open(ref_dir_domain + query + '.clean.fasta', 'w') as seq_out:
            for record in SeqIO.parse(ref_dir_domain + ref + '.clean.fasta', 'fasta'):
                i += 1
                if i == 1:
                    record.id = 'test_' + str(record.id)
                    record.description = ''
                    SeqIO.write(record, seq_out, 'fasta')
        
    for fasta in fastas.keys():
        if fastas[fasta] > 3:
        
            query_align(ref_dir_domain + query + '.clean.fasta',
                ref_dir_domain + fasta + '.clean.align.sto',
                ref_dir_domain + query + '.' + fasta + '.clean.align.sto',
                cm)
            
            temp_dir = tempfile.mkdtemp(dir = ref_dir_domain)
            
            split_query_ref(ref_dir_domain + query + '.clean.align.sto',
                ref_dir_domain + fasta + '.clean.align.sto',
                ref_dir_domain + query + '.' + fasta + '.clean.align.sto',
                temp_dir)
    
            place(query + '.clean.align.newlength.fasta',
                fasta + '.clean.align.newlength.fasta',
                ref_dir_domain + fasta + '.final.bestModel',
                ref_dir_domain + fasta + '.final.bestTree',
                ref_dir_domain + fasta + '.jplace',
                temp_dir)
            
            shutil.rmtree(temp_dir)
        
    ## For all domains, create a file with the date/time of database creation.
    #!!! n_aseqs value not correct for all domains

    current_time = datetime.datetime.now().isoformat()
    
    with open(ref_dir_domain + ref + '.database_info.txt', 'w') as database_info:
        print('ref tree built at:', current_time, file = database_info)
        print('nseqs in reference alignment:', n_aseqs, file = database_info) 
        
        for fasta in fastas.keys():
            print('*' + fasta + '\t' + str(fastas[fasta]), file = database_info)
            
else:
    
    ## If the query flag is give you are placing reads on the existing reference tree.
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
                line = line.split('.')[-1]
                line = line.split('\t')
                if int(line[1]) > 3:
                    available_trees.add(line[0])
                    
    ## Chose an appropriate name for top level reference tree based on
    ## domain.  phylum_ref is really only an appropriate variable name
    ## for domain Bacteria, but it does no harm.
                    
    top_level_names = {'bacteria':ref + '.phylum_reps',
                       'archaea':ref,
                       'eukarya':ref + '.div_reps'}
    
    phylum_ref = top_level_names[domain]
            
    clear_wspace = subprocess.call('rm -f ' + query + '.' + phylum_ref + '*',
                                   shell = True,
                                   executable = executable)
    
    clean_name(cwd + query + '.fasta')
    make_unique(cwd + query + '.clean.fasta')

    query_align(cwd + query + '.clean.unique.fasta',
          ref_dir_domain + phylum_ref + '.clean.align.sto',
          cwd + query + '.' + phylum_ref + '.clean.unique.align.sto',
          cm)

    temp_dir = tempfile.mkdtemp(dir = cwd)

    split_query_ref(cwd + query + '.clean.unique.align.sto',
          ref_dir_domain + phylum_ref + '.clean.align.sto',
          cwd + query + '.' + phylum_ref + '.clean.unique.align.sto',
          temp_dir)

    place(query + '.clean.unique.align.newlength.fasta',
          phylum_ref + '.clean.align.newlength.fasta',
          ref_dir_domain + phylum_ref + '.final.bestModel',
          ref_dir_domain + phylum_ref + '.final.bestTree',
          cwd + query + '.' + phylum_ref + '.jplace',
          temp_dir)
    
    placements = json_to_csv(query + '.' + phylum_ref + '.jplace', phylum_ref)
    
    placements = get_map_ratio(temp_dir + '/' + query + '.clean.unique.align.newlength.fasta',
        temp_dir + '/' + phylum_ref + '.clean.align.newlength.fasta',
        placements)
                    
    gappa(query + '.' + phylum_ref + '.jplace', cwd)  
    edpl = pd.read_csv(cwd + query + '.' + phylum_ref + '.edpl.list.csv', index_col = 1)
    placements = pd.concat([placements, edpl['EDPL']], axis = 1, sort = False)
    placements.to_csv(cwd + query + '.' + ref + '.placements.csv') 
    
    shutil.rmtree(temp_dir)
    
    if domain in ['bacteria', 'eukarya']:
        
        ## If the domain is either bacteria or eukarya then placement to
        ## subtrees is necessary.
        
        for subtree in placements.ref_name.unique():
            if pd.notnull(subtree):
                
                ## The phylum_ref variable is named as such because the
                ## domain Bacteria subtrees are constructed at the phylum
                ## level.  This is not the case for Eukarya, but it doesn't 
                ## matter here.
                
                phylum_ref = ref + '.' + subtree
                
                with open(cwd + phylum_ref + '_' + query + '.clean.unique.fasta', 'w') as subtree_out:
                    for record in SeqIO.parse(query + '.clean.unique.fasta', 'fasta'):
                        if placements.loc[record.id, 'ref_name'] == subtree:
                            SeqIO.write(record, subtree_out, 'fasta')  
            else:
                with open(cwd + ref + 'nosubtree_' + query + '.clean.unique.fasta', 'w') as subtree_out:
                    for record in SeqIO.parse(query + '.clean.unique.fasta', 'fasta'):
                        if pd.isnull(placements.loc[record.id, 'ref_name']):
                            SeqIO.write(record, subtree_out, 'fasta') 
                            
        combined_subtrees = pd.DataFrame()
                        
        for subtree in placements.ref_name.unique():
            if subtree in available_trees:
                if pd.notnull(subtree):
                    
                    phylum_ref = ref + '.' + subtree
                    
                    query_align(cwd + phylum_ref + '_' + query + '.clean.unique.fasta',
                          ref_dir_domain + phylum_ref + '.clean.align.sto',
                          cwd + query + '.' + phylum_ref + '.clean.unique.align.sto',
                          cm)
                            
                    temp_dir = tempfile.mkdtemp(dir = cwd)
                    
                    split_query_ref(cwd + phylum_ref + '_' + query + '.clean.unique.align.sto',
                          ref_dir_domain + phylum_ref + '.clean.align.sto',
                          cwd + query + '.' + phylum_ref + '.clean.unique.align.sto',
                          temp_dir)
                    
                    place(phylum_ref + '_' + query + '.clean.unique.align.newlength.fasta',
                          phylum_ref + '.clean.align.newlength.fasta',
                          ref_dir_domain + phylum_ref + '.final.bestModel',
                          ref_dir_domain + phylum_ref + '.final.bestTree',
                          cwd + query + '.' + phylum_ref + '.jplace',
                          temp_dir)
                    
                    subtree_csv = json_to_csv(cwd + query + '.' + phylum_ref + '.jplace', phylum_ref)
                    
                    subtree_csv = get_map_ratio(temp_dir + '/' + phylum_ref + '_' + query + '.clean.unique.align.newlength.fasta',
                                                temp_dir + '/' + phylum_ref + '.clean.align.newlength.fasta',
                                                subtree_csv)
                                
                    gappa(query + '.' + phylum_ref + '.jplace', cwd)
                    edpl = pd.read_csv(cwd + query + '.' + phylum_ref + '.edpl.list.csv', index_col = 1)
                    subtree_csv = pd.concat([subtree_csv, edpl['EDPL']], axis = 1, sort = False)
                    subtree_csv['subtree'] = subtree
                    
                    combined_subtrees = pd.concat([combined_subtrees, subtree_csv])
                    subtree_csv.to_csv(cwd + query + '.' + ref + '.placements.csv') 
              
                    shutil.rmtree(temp_dir)
                    
        ## Identify all unique reads that were assigned to a clade that didn't
        ## have > 3 reference sequences, meaning that there is not reference tree.
                    
        no_subtrees = placements.index.difference(combined_subtrees.index)
        combined_subtrees = pd.concat([combined_subtrees, placements.loc[no_subtrees]], sort = False)               
        
        ## Add the count data.  This should also add all unique reads that didn't
        ## place the original reference tree.
            
        count = pd.read_csv(cwd + query + '.clean.unique.count', index_col = 0)
        combined_subtrees = pd.concat([combined_subtrees, count], axis = 1, sort = False)
        combined_subtrees.to_csv(cwd + query + '.' + ref + '.placements.csv')
        
#!! edge_lineages still indexed by taxid, which seems inefficient, why not edge number?
