#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Tue Jan 06 09:50:07 2015

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
    Programs:
        Pathway-tools
        
    Python modules
        Pandas
        Bio
        Numpy
        Joblib
        
RUN AS:
    python paprica-build_core_genomes.py -tree [tree.phyloxml] -domain [bacteria|archaea|eukarya]
    
OPTIONS:
    -domain: The domain being analyzed (either bacteria, archaea, or eukarya)
    -pgdb_dir: The location where pathway-tools stores PGDBs
    -ref_dir: The directory containing the paprica database
    -tree: The phyloxml format tree that contains the clade numbers
    -cpus: The number of parallel calls to make to pathologic.  Defaults to 4,
        use -1 for all available.

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
from Bio import SeqIO
from Bio import Phylo

import subprocess
import os
import re
import sys
import shutil
from joblib import Parallel, delayed
import json
from io import StringIO

import pandas as pd
import numpy as np

executable = '/bin/bash'

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print('Manually stopped!')
    print(stop[1])
               
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
        
## Define some variables based on these arguments.  If nothing in sys.argv
## set some default values.  This is useful for testing.
        
if len(sys.argv) == 1:
    domain = 'bacteria'
    ref_dir = 'ref_genome_database'
    pgdb_dir = '/volumes/hd2/ptools-local/pgdbs/user/'
    cpus = 36
    database_info = 'combined_16S.bacteria.tax.database_info.txt'
    
else:        
    domain = command_args['domain']
    ref_dir = command_args['ref_dir']
    pgdb_dir = command_args['pgdb_dir']
    cpus = int(command_args['cpus'])
    database_info = command_args['database_info']
    
## Expand tilde manually.
    
pgdb_dir = os.path.expanduser(pgdb_dir)
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'

paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/' # The location of the actual paprica scripts.   
ref_dir_domain = paprica_path + ref_dir + domain + '/'

#%% Define some functions.

## Define a function to classify each node in a reference tree. Lineage is the
## df containing the edge_lineage data.  Tree is the tree object for the
## reference tree for which lineage is needed.  Ranks are the ranks desired
## for the domain.

def get_lineages_terminal(lineage, ranks,
                          genome_data, clade_number,
                          assembly):
    
    temp_lineages = pd.DataFrame(columns = ['consensus'] + ranks)
           
    print('getting lineage for', clade_number)

    temp_taxid = genome_data.loc[assembly, 'taxid']
    temp_lineage = lineage.loc[temp_taxid]
    
    ## Identify a meaningful lowest rank.  This will become "consensus".
    
    for i in range(1, len(ranks)):
        rank = ranks[-1 * i]
        consensus_taxa = temp_lineage[rank]
        if pd.notnull(consensus_taxa):
            break
        
    ##consensus_lineage['consensus'] = consensus_lineage[consensus_lineage.notnull()][-1]
        
    temp_lineages.loc[clade_number, 'consensus'] = consensus_taxa
    temp_lineages.loc[clade_number, ranks] = temp_lineage[ranks]
    
    return(temp_lineages)

def get_lineages_nonterminal(clade, subtree, genome_data, ranks, lineage):
            
    clade_number = subtree + '_' + clade.comment
    print('getting lineage for', clade_number)
    terminals = []
    
    ## Extract branch length here so that you don't have to loop
    ## over tree later in get_internals
    
    branch_length = clade.branch_length
    
    for terminal in clade.get_terminals():
        assembly = terminal.name.split('|')[0]
        terminals.append(assembly)
        
    temp_taxids = genome_data.loc[terminals, 'taxid']       
    temp_lineage = lineage.loc[temp_taxids]        
    
    ## Now iterate across columns, starting at superkingdom
    ## until you find the first mismatch.  The one
    ## before this is the consensus.
    
    found_nonconsensus = False
    consensus_lineage = pd.Series(index = ranks, name = clade_number)
    
    for i,rank in enumerate(temp_lineage.columns):
        
        if len(temp_lineage[rank].unique()) == 1:
            
            ## There is consensus at this level. Record the taxon name
            ## associated with this rank.
            
            consensus_lineage[rank] = temp_lineage[rank].unique()[0]
            
        else:
            
            ## There is not consensus at this level. Look up consensus
            ## taxonomy from previous level.
            
            consensus_rank = consensus_lineage.index[i - 1]
            consensus_lineage['consensus'] = consensus_lineage[consensus_rank]
            found_nonconsensus = True
            
            break
        
    ## In some cases the nodes in the clade will be identical. Pull out the lowest
    ## rank that is not NaN.
        
    if found_nonconsensus == False:
        consensus_lineage['consensus'] = consensus_lineage[consensus_lineage.notnull()][-1]
        
    return consensus_lineage, terminals, clade_number, branch_length
    
## Define a function to use pathway-tools to predict the metabolic pathways for
## each genome downloaded from Genbank, or transcriptome downloaded from MMTSP.
## The PGDBs are now named by assembly so that they can be re-used for each new
## version of the database.

def make_pgdb(d, ref_dir_domain):
    
    print(d, 'start prediction')
    
    subprocess.call('pathway-tools \
    -lisp \
    -no-cel-overview \
    -patho ' + ref_dir_domain + 'refseq/' + d + '/ \
    -disable-metadata-saving \
    &> ' + ref_dir_domain + 'pathos_' + d + '.log', shell = True, executable = executable)
     
    print(d, 'prediction complete')  
    
def parse_pwy_inference(pwy_inference_report, pathway_definitions):
    
    path = pwy_inference_report.split('/')[0:-1]
    path = '/'.join(path)
    keep = False
    pathways_found = []
    
    with open(pwy_inference_report, 'r') as report_in:
        for line in report_in:
            if keep == True:
                line = line.rstrip()
                
                if line.endswith(')'):
                    ## That means it's the ends of the pathways kept block.
                    keep = False
                    
                line = line.strip('(')
                line = line.rstrip(')')
                line = line.rstrip()
                line = line.split()
                pathways_found = pathways_found + line
            elif line.startswith('List of pathways kept'):
                keep = True
            
    with open(path + '/pathways-report.txt', 'w') as report_out:
        print('Pathway Name | Pathway Frame-id', file = report_out)
        for pathway in pathways_found:
            print(pathway_definitions.loc[pathway, 'NAME'] + '|' + pathway, file = report_out)
             
#%% Preparatory file generation and organization.
    
if domain == 'eukarya':
    ref = 'combined_18S.eukarya.tax'
else:
    ref = 'combined_16S.23S.' + domain + '.tax'
    
## Read in database information, this will indicate which subtrees need
## to be evaluated. 
    
available_trees = []

with open(ref_dir_domain + database_info, 'r') as database_file:
    for line in database_file:
        line = line.rstrip()
        if line.startswith('*'):
            line = line.strip('*')
            line = line.split('\t')
            if int(line[1]) > 3:
                available_trees.append(line[0])
                
## The top level tree should be excluded from this list
## for the bacteria and eukarya. Here we presume that the first entry is
## the top level tree which is suboptimal.
                
if domain in ['bacteria', 'eukarya']:
    available_trees = available_trees[1:-1]
                
## Collect information for determining lineage of each node in the reference
## tree.

if domain == 'eukarya':
    ranks = ['kingdom', 'supergroup', 'division', 'class', 'order', 'family',
             'genus', 'species']
else:
    ranks = ['superkingdom', 'phylum', 'clade', 'class',
         'order', 'family', 'genus', 'species', 'strain']
    
lineage = pd.read_csv(ref_dir_domain + 'edge_lineages.csv', index_col = 0)
lineage = lineage[ranks]
    
node_lineages = pd.DataFrame(columns = ['consensus'] + ranks)

## Read in the genome_data file.

genome_data = pd.read_csv(ref_dir_domain + 'genome_data.csv.gz', header = 0, index_col = 0)
genome_data['clade'] = np.nan
genome_data['tip_name'] = np.nan
genome_data['npaths_actual'] = np.nan
genome_data['branch_length'] = np.nan

## Iterate across jplace files here.  These come from available_trees.

assemblies = []
internal_terminals = {}
internal_branch_length = {}

for jplace in available_trees:
    
    ## Get the clade number of each assembly and add this information to 
    ## genome_data.
    
    with open(ref_dir_domain + ref + '.' + jplace + '.jplace', 'r') as jfile:
        data = json.load(jfile)
        colnames = data['fields']
    
    tree = data['tree']
    tree = re.sub('{', '[', tree)
    tree = re.sub('}', ']', tree)
    tree = Phylo.read(StringIO(tree), 'newick') 
    
    for clade in tree.get_terminals():
        ref_name = clade.name
        assembly = clade.name.split('|')[0]
        assemblies.append(assembly)
        
        subtree = genome_data.loc[assembly, 'subtree']       
        clade_number = subtree + '_' + clade.comment
        
        genome_data.loc[assembly, 'clade'] = clade_number
        genome_data.loc[assembly, 'tip_name'] = clade.name
        genome_data.loc[assembly, 'branch_length'] = clade.branch_length
        
        terminal_lineages = get_lineages_terminal(lineage = lineage, ranks = ranks,
                          genome_data = genome_data, clade_number = clade_number,
                          assembly = assembly)
        
        node_lineages = pd.concat([node_lineages, terminal_lineages])
        
    for clade in tree.get_nonterminals():
        if clade.comment != None:
            internal_lineages, clade_terminals, clade_number, branch_length = get_lineages_nonterminal(clade = clade,
                                                                                        subtree = subtree,
                                                                                        genome_data = genome_data,
                                                                                        ranks = ranks,
                                                                                        lineage = lineage)
            
            node_lineages = node_lineages.append(internal_lineages)
            internal_terminals[clade_number] = clade_terminals
            internal_branch_length[clade_number] = branch_length
        
## Write out lineage files, then exit if eukarya.
        
node_lineages.to_csv(ref_dir_domain + 'node_lineages.csv.gz') 

if domain == 'eukarya':
    genome_data.to_csv(ref_dir_domain + 'genome_data.final.csv.gz')
    quit()
    
#%% Prep for building PGDBs.
        
## For every existing PGDB directory (which might be none), determing if the
## pathways-report.txt file is present.  If it is not this means the previous
## pathway-tools pathologic run was incorrect.  In that case delete the
## directory and try again.
    
new_pgdbs = []
pathway_definitions = pd.read_csv(ref_dir + 'pathways.col', comment = '#', sep = '\t', index_col = 0)

#!!! This loop does take a long time and would be easy to parallelize

for i, d in enumerate(assemblies):
    
    ## If a previous build effort was unsuccessful rewrite the files needed by
    ## pathologic, in case the data files have been updatated in the public
    ## repository (Genbank or MMETSP) and this fixes the problem.  
    
    ## If an assembly was updated by NCBI the PGDB deleted by paprica-make_ref.py.
    ## This will force recreation of data files here.
    
    report_file = False
    
    try:
        
        for f in os.listdir(pgdb_dir + d.lower() + 'cyc/1.0/reports'):
            
            ## If an old version of pathways-report.txt is present remove it.
            ## This applies also to ad-hoc reports created by this script while
            ## we wait for a realistic solution to the pathologic errors. 
            
            #!!! currently the "old" style pathways-report file is not deleted.
            ## This is to prevent the current issue with pathologic.
            
            # if f == 'pathways-report.txt':
            #     os.remove(pgdb_dir + d.lower() + 'cyc/1.0/reports/' + f)
            #     report_file = True
                
            if f.startswith('pathways-report'):
                report_file = True
                                
        ## If there was no pathways-report file, or it was the old format and deleted,
        ## rewrite the files needed by pathologic.
                
        if report_file == False:
                            
            with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
                
                print('recreating pathologic files for', d, i, 'of', len(assemblies))
                
                print('ID' + '\t' + d, file=organism_params)
                print('Storage' + '\t' + 'File', file=organism_params)
                print('Name' + '\t' + d, file=organism_params)
                print('Rank' + '\t' + 'Strain', file=organism_params)
                print('Domain' + '\t' + 'TAX-2', file=organism_params)
                print('Create?' + '\t' + 't', file=organism_params)
                
                g = 0
                
                for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                    if gbk.endswith('genomic.gbk'):
                        
                        ## NCBI PGAP packages genetic elements in a single
                        ## genbank file which ptools doesn't like.  Split the
                        ## file.
                        
                        for record in SeqIO.parse(ref_dir_domain + 'refseq/' + d + '/' + gbk, 'genbank'):
                            
                            g = g + 1
                            basename = re.split('gbk', gbk)[0]
                            with open(ref_dir_domain + 'refseq/' + d + '/' + basename + str(g) + '.gbk', 'w') as seq_out:
                                SeqIO.write(record, seq_out, 'genbank')
                        
                            print('ID' + '\t' + d + '.' + str(g), file=genetic_elements)
                            print('NAME' + '\t' + d + '.' + str(g), file=genetic_elements)
                            print('TYPE' + '\t' + ':CHRSM', file=genetic_elements)
                            
                            if domain == 'eukarya':
                                print('CIRCULAR?' + '\t' + 'Y', file=genetic_elements)
                            else:
                                print('CIRCULAR?' + '\t' + 'N', file=genetic_elements)
                                
                            print('ANNOT-FILE' + '\t' + basename + str(g) + '.gbk', file=genetic_elements)
                            print('//', file=genetic_elements)
                        
            if g > 0:
                new_pgdbs.append(d)
                
    ## If there was no previous build attempt the directory will not exist and
    ## os will throw an error.
            
    except FileNotFoundError:
        
        print('creating pathologic files for', d, i, 'of', len(assemblies))
                    
        with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
                        
            print('ID' + '\t' + d, file=organism_params)
            print('Storage' + '\t' + 'File', file=organism_params)
            print('Name' + '\t' + d, file=organism_params)
            print('Rank' + '\t' + 'Strain', file=organism_params)
            print('Domain' + '\t' + 'TAX-2', file=organism_params)
            print('Create?' + '\t' + 't', file=organism_params)
            
            g = 0
            
            for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                if gbk.endswith('genomic.gbk'):
                    
                    ## NCBI PGAP packages genetic elements in a single
                    ## genbank file which ptools doesn't like.  Split the
                    ## file.
                    
                    for record in SeqIO.parse(ref_dir_domain + 'refseq/' + d + '/' + gbk, 'genbank'):
                        
                        g = g + 1
                        basename = re.split('gbk', gbk)[0]
                        with open(ref_dir_domain + 'refseq/' + d + '/' + basename + str(g) + '.gbk', 'w') as seq_out:
                            SeqIO.write(record, seq_out, 'genbank')
                    
                        print('ID' + '\t' + d + '.' + str(g), file=genetic_elements)
                        print('NAME' + '\t' + d + '.' + str(g), file=genetic_elements)
                        print('TYPE' + '\t' + ':CHRSM', file=genetic_elements)
                        
                        if domain == 'eukarya':
                            print('CIRCULAR?' + '\t' + 'Y', file=genetic_elements)
                        else:
                            print('CIRCULAR?' + '\t' + 'N', file=genetic_elements)
                            
                        print('ANNOT-FILE' + '\t' + basename + str(g) + '.gbk', file=genetic_elements)
                        print('//', file=genetic_elements)
                    
        if g > 0:
            new_pgdbs.append(d)
                
## Previous failed builds confuse pathway-tools.  Remove the directory.

print('removing previous failed PGDBs...')
                
for d in new_pgdbs:
    shutil.rmtree(pgdb_dir + d.lower() + 'cyc', ignore_errors = True)
    
#%% Generate the PGDBs for each assembly that does not have one.

print(len(new_pgdbs), 'new pgdbs will be created')

Parallel(n_jobs = cpus, verbose = 5)(delayed(make_pgdb)
(d, ref_dir_domain) for d in new_pgdbs)
    
#%% For each PGDB add the pathways to a new data_frame.

terminal_paths = pd.DataFrame(index = assemblies)

for i,d in enumerate(assemblies):
    n_paths = 0 # Number of pathways predicted.
    try:
        
        ## As of ptools v24 pathway-report.txt reformatted with date.  Need
        ## to find name of this file.
        
        report_file = None
        
        for f in os.listdir(pgdb_dir + d.lower() + 'cyc/1.0/reports'):
            if f.startswith('pathways-report_'):
                report_file = f
                print(f)
            elif f.startswith('pwy-inference-report'):
                inference_file = f
                
        #!!! For reasons that aren't clear, pathway-tools is failing at the
        ## end of the prediction, so the pathways-report file isn't being
        ## created.  However, the pwy-inference-report is created, and
        ## this can be parsed to create the pathways-report. This pathway
        ## report has the old-style name, so it is deleted each time, in the
        ## hope that eventually pathologic works correctly.
                
        if report_file == None:
            parse_pwy_inference(pgdb_dir + d.lower() + 'cyc/1.0/reports/' + inference_file, pathway_definitions)
            report_file = 'pathways-report.txt'
                
        ## Now the report file can be parsed.
        
        with open(pgdb_dir + d.lower() + 'cyc/1.0/reports/' + report_file, 'r') as report:            
            for line in report:
                if line.startswith('#') == False:
                    if line.startswith('Pathway Name') == False:
                        
                        ## PWY-WAS-NOT-DELETED no longer valid, but not doing any harm.
                        
                        if 'PWY-WAS-NOT-DELETED' not in line:
                            line = line.rstrip()
                            if line != '':
                                line = line.split('|')                                
                                if len(line) > 0:
                                    path = line[0]
    
                                    try:
                                        terminal_paths.loc[d, path] = 1
                                    except KeyError:
                                        terminal_paths[path] = np.nan
                                        terminal_paths.loc[d, path] = 1
                                        
                                    print('collecting paths for terminal node', d, i + 1, 'of', len(assemblies), path)
                                    n_paths = n_paths + 1
                                    
        genome_data.loc[d, 'npaths_actual'] = n_paths
        
    except (NameError, IOError):
        print(d, 'has no pathway report')
        
#%% Collect EC_numbers for each terminal node

## Read in user specified EC numbers.
        
user_ec = pd.read_csv(ref_dir + 'user/' + 'user_ec.csv', header = 0, index_col = 0, comment = '#')
        
terminal_ec = pd.DataFrame(index = assemblies)
#ec_names = pd.DataFrame()

## !!! This loop needs to be parallelized
        
for i,d in enumerate(assemblies):
    n_paths = 0
    
    for f in os.listdir(ref_dir_domain + '/refseq/' + d):
        if f.endswith('gbff'):
            try:
                for record in SeqIO.parse(ref_dir_domain + '/refseq/' + d + '/' + f, 'genbank'):
                    for feature in record.features:
                        if feature.type == 'CDS':
                            
                            try:
                                protein_id = feature.qualifiers['protein_id'][0]
                            except KeyError:
                                protein_id = 'no protein_id'
                            
                            if 'EC_number' in list(feature.qualifiers.keys()):
                                
                                n_paths = n_paths + 1
                                ec = feature.qualifiers['EC_number']
                                
                                ## Some draft genomes will not have a product qualifier.
                                
                                try:
                                    prod = feature.qualifiers['product'][0]
                                except KeyError:
                                    prod = 'product not specified'
                                
                                ## Because each EC number can appear multiple times
                                ## in a genome this information needs to be tallied.
                                
                                for each in ec:
                                    print('collecting EC numbers for terminal node', d, i + 1, 'of', str(len(assemblies)) + ',', protein_id + ':', each)
                                    
                                    try:
                                        temp = terminal_ec.loc[d, each]
                                        if pd.isnull(temp) == True:
                                            terminal_ec.loc[d, each] = 1
                                        else:
                                            terminal_ec.loc[d, each] = temp + 1
                                    except KeyError:
                                        terminal_ec.loc[d, each] = 1
                                        
#                                    ec_names.loc[each, 'name'] = prod

            ## For some assemblies an error is raised on a second(?) record identified
            ## in the Genbank file.  It isn't clear why this is happening, pass the error
            ## here.
            
            except AttributeError:
                pass
        
    ## For assembly d, add in any user specified EC numbers.
    
    if d in set(user_ec['GI_number']):
        temp_user_ec = user_ec[user_ec['GI_number'] == d]
        
        for entry in temp_user_ec.index:
            n_paths = n_paths + 1
            each = temp_user_ec.loc[entry, 'EC_number']
            print('collecting EC numbers for terminal node', d, i + 1, 'of', len(assemblies), each)
                                
            try:
                if pd.isnull(terminal_ec.loc[d, each]) == True:
                    terminal_ec.loc[d, each] = 1
                else:
                    print('There is already an entry for', d, each)
            except KeyError:
                terminal_ec.loc[d, each] = 1
                
            ## Check to make sure that the enzyme number has a name, add if it does not
                
#            try:
#                if pd.isnull(ec_names.loc[each, 'name']):
#                    ec_names.loc[each, 'name'] = temp_user_ec.loc[entry, 'product']
#            except KeyError:
#                ec_names.loc[each, 'name'] = temp_user_ec.loc[entry, 'product']

    ## Add the total number of enzymes for that assembly to genome_data.

    genome_data.loc[d, 'nec_actual'] = n_paths

#%% Collect pathway and other data for internal nodes.
## Make an initial pass over the tree to collect all the internal node numbers.

int_nodes = sorted(set(internal_terminals.keys()))
n_clades = len(int_nodes)

## Create a new dataframes to store pathway data, the fraction of daughters
## with each pathway, and genome data inferred from daughters for each internal
## node.

internal_probs_columns = list(terminal_paths.columns)
internal_data_columns = ['n16S', 'nge', 'ncds', 'genome_size', 'GC', 'phi', 'clade_size', 'npaths_terminal', 'nec_terminal', 'branch_length']
internal_ec_probs_columns = terminal_ec.columns
internal_ec_n_columns = terminal_ec.columns

internal_probs = np.memmap(open('internal_probs.mmap', 'w+b'), shape = (n_clades, len(internal_probs_columns)), dtype = 'f8')
internal_data = np.memmap(open('internal_data.mmap', 'w+b'), shape = (n_clades, len(internal_data_columns)), dtype = 'f8')
internal_ec_probs = np.memmap(open('internal_ec_probs.mmap', 'w+b'), shape = (n_clades, len(internal_ec_probs_columns)), dtype = 'f8')
internal_ec_n = np.memmap(open('internal_ec_n.mmap', 'w+b'), shape = (n_clades, len(internal_ec_n_columns)), dtype = 'f8')

## Define a function to iterate across all subtrees and collect information that will be saved in
## the "internal" memory-mapped arrays. Remember that these are not dataframes
## and cannot be indexed by column/row names!

def get_internals(clade_number,
                  internal_terminals,
                  internal_branch_length,
                  int_nodes,
                  genome_data,
                  terminal_paths,
                  terminal_ec,
                  internal_probs_columns,
                  internal_data_columns,
                  internal_ec_probs_columns,
                  internal_ec_n_columns,
                  internal_probs,
                  internal_data,
                  internal_ec_probs,
                  internal_ec_n):
        
    print('collecting data for internal node', clade_number)
    
    ## Data on the clade that you want later.
    
    clade_members = internal_terminals[clade_number]
    edge_i = int_nodes.index(clade_number)
    internal_data[edge_i, internal_data_columns.index('branch_length')] = internal_branch_length[clade_number]
    ntip = len(clade_members)
        
    ## Get data for all clade members.
        
    clade_data = genome_data.loc[clade_members]
    clade_paths = terminal_paths.loc[clade_members]
    clade_ec = terminal_ec.loc[clade_members]
        
    npaths = clade_paths.count(axis = 0, numeric_only = True)
    rpaths = npaths.div(ntip)
    internal_probs[edge_i, :] = rpaths
    
    nec = clade_ec.count(axis = 0, numeric_only = True)
    rec = nec.div(ntip)
    internal_ec_probs[edge_i, :] = rec
    
    ## The mean number of occurrences of EC_numbers are calculated for clade
    ## so that later on an estimate can be given of the number that will appear
    ## in an internal node.
    
    mec = clade_ec.mean(axis = 0, numeric_only = True)
    internal_ec_n[edge_i, :] = mec
    
    ## Calculate values for this edge.
        
    internal_data[edge_i, internal_data_columns.index('n16S')] = clade_data['n16S'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('nge')] = clade_data['nge'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('ncds')] = clade_data['ncds'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('genome_size')] = clade_data['genome_size'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('phi')] = clade_data['phi'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('GC')] = clade_data['GC'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('clade_size')] = ntip
    internal_data[edge_i, internal_data_columns.index('npaths_terminal')] = clade_data['npaths_actual'].dropna().mean()
    internal_data[edge_i, internal_data_columns.index('nec_terminal')] = clade_data['nec_actual'].dropna().mean()
        
#%% Execute the get_internals function in parallel, this massively speeds up
## the database build.
        
## For bacteria currently throws RuntimeError: maximum recursion depth exceeded
## if more than 1 processor used.  This is kind of a problem because the
## bacteria database is the one that takes a long time to build...
        
if domain == 'bacteria':
    njobs = -1
else:
    njobs = -1
    
## now iterating across internal_terminals.keys()
      
Parallel(n_jobs = njobs)(delayed(get_internals)
                         (clade_number,
                          internal_terminals,
                          internal_branch_length,
                          int_nodes,
                          genome_data,
                          terminal_paths,
                          terminal_ec,
                          internal_probs_columns,
                          internal_data_columns,
                          internal_ec_probs_columns,
                          internal_ec_n_columns,
                          internal_probs,
                          internal_data,
                          internal_ec_probs,
                          internal_ec_n)
for clade_number in internal_terminals.keys())
   
## Write out ya database files.
            
genome_data.to_csv(ref_dir_domain + 'genome_data.final.csv.gz')
terminal_paths.to_csv(ref_dir_domain + 'terminal_paths.csv.gz')
terminal_ec.to_csv(ref_dir_domain + 'terminal_ec.csv.gz')

internal_data = pd.DataFrame(internal_data, index = int_nodes, columns = internal_data_columns)
internal_data.to_csv(ref_dir_domain + 'internal_data.csv.gz')

internal_probs = pd.DataFrame(internal_probs, index = int_nodes, columns = internal_probs_columns)
internal_probs.to_csv(ref_dir_domain + 'internal_probs.csv.gz')

internal_ec_probs = pd.DataFrame(internal_ec_probs, index = int_nodes, columns = internal_ec_probs_columns)
internal_ec_probs.to_csv(ref_dir_domain + 'internal_ec_probs.csv.gz')

internal_ec_n = pd.DataFrame(internal_ec_n, index = int_nodes, columns = internal_ec_n_columns)
internal_ec_n.to_csv(ref_dir_domain + 'internal_ec_n.csv.gz')

## Clean up memory maps

os.system('rm *mmap')
