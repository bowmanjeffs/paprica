# -*- coding: utf-8 -*-
help_string = """
Created on Tue Jan 06 09:50:07 2015

@author: jeff

REQUIRES:
    Programs:
        Pathway-tools
        
    Python modules
        Pandas
        Bio
        Numpy
        
RUN AS:
    python paprica_build_core_genomes_v0.20.py -tree [tree.phyloxml] -domain [bacteria|archaea]
    
OPTIONS:
-pgdb_dir: The location where pathway-tools stores PGDBs
-ref_dir: The directory containing the paprica database
-tree: The phyloxml format tree that contains the clade numbers
-domain: The domain being analyzed (either bacteria or archaea)

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
from Bio import Phylo

import subprocess
import os
import re
import sys
import shutil
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

executable = '/bin/bash'

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]
               
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
        
## Define some variables based on these arguments.
        
if len(sys.argv) == 1:
    domain = 'archaea'
    tree = 'test.archaea.combined_16S.archaea.tax.clean.align.phyloxml'
    ref_dir = 'ref_genome_database'
    pgdb_dir = '/volumes/hd1/ptools-local/pgdbs/user'
    
else:        
    domain = command_args['domain']
    tree = command_args['tree']
    ref_dir = command_args['ref_dir']
    pgdb_dir = command_args['pgdb_dir']

ref_dir_domain = ref_dir + '/' + domain + '/'

## Read in the genome_data file.

genome_data = pd.DataFrame.from_csv(ref_dir_domain + 'genome_data.csv', header = 0, index_col = 0)
genome_data['clade'] = np.nan
genome_data['tip_name'] = np.nan
genome_data['npaths_actual'] = np.nan
genome_data['branch_length'] = np.nan
    
## Get the clade number of each assembly and add this information to 
## genome_data.
    
tree = Phylo.read(tree, 'phyloxml')
    
assemblies = []
    
for clade in tree.get_terminals():
    clade_number = int(clade.confidence)
    print clade_number
    
    assembly = clade.name
    assembly = assembly.strip('@')
    
    assembly = re.split('_', assembly)
    assembly = assembly[0] + '_' + assembly[1]
    
    genome_data.loc[assembly, 'clade'] = clade_number
    genome_data.loc[assembly, 'tip_name'] = clade.name
    genome_data.loc[assembly, 'branch_length'] = clade.branch_length
    
    assemblies.append(assembly)

#%% Generate the PGDBs for each assembly that does not have one.

## Use pathway-tools to predict the metabolic pathways for each genome downloaded
## from Genbank.  The PGDBs are now named by assembly so that they can be re-used
## for each new version of the database.

def make_pgdb(d, ref_dir_domain):
    
    print d, 'start prediction'
    predict_pathways = subprocess.Popen('pathway-tools -lisp -no-cel-overview -patho ' + ref_dir_domain + 'refseq/' + d + '/ -disable-metadata-saving &> ' + ref_dir_domain + 'pathos_' + d + '.log', shell = True, executable = executable)
    predict_pathways.communicate()   
    print d, 'prediction complete'
    
pgdbs = set(os.listdir(pgdb_dir))
new_pgdbs = []

for d in assemblies:
    
    ## If a previous build effort was unsuccessful, try again in case the files
    ## have been updatated in Genbank and this fixes the problem.    
    
    try:
        if 'pathways-report.txt' not in os.listdir(pgdb_dir + d.lower() + 'cyc/1.0/reports'):
            
            clade = genome_data.loc[d, 'clade']
                
            with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
        
                print >> organism_params, 'ID' + '\t' + d
                print >> organism_params, 'Storage' + '\t' + 'File'
                print >> organism_params, 'Name' + '\t' + d
                print >> organism_params, 'Rank' + '\t' + 'Strain'
                print >> organism_params, 'Domain' + '\t' + 'TAX-2'
                print >> organism_params, 'Create?' + '\t' + 't'
                
                g = 0
                
                for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                    if gbk.endswith('gbff'):
                        g = g + 1
                        
                        basename = re.split('gbff', gbk)[0]
                        subprocess.call('cd ' + ref_dir_domain + 'refseq/' + d + ';cp ' + gbk + ' ' + basename + 'gbk', shell = True, executable = executable)
                        
                        print >> genetic_elements, 'ID' + '\t' + d + '.' + str(g)
                        print >> genetic_elements, 'NAME' + '\t' + d + '.' + str(g)
                        print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                        print >> genetic_elements, 'CIRCULAR?' + '\t' + 'Y'
                        print >> genetic_elements, 'ANNOT-FILE' + '\t' + basename + 'gbk'
                        print >> genetic_elements, '//'
                        
            if g > 0:
                new_pgdbs.append(d)
                
    ## If there was no previous build attempt the directory will not exist and
    ## os wil throw an error.
            
    except OSError:
        
        clade = genome_data.loc[d, 'clade']
            
        with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
    
            print >> organism_params, 'ID' + '\t' + d
            print >> organism_params, 'Storage' + '\t' + 'File'
            print >> organism_params, 'Name' + '\t' + d
            print >> organism_params, 'Rank' + '\t' + 'Strain'
            print >> organism_params, 'Domain' + '\t' + 'TAX-2'
            print >> organism_params, 'Create?' + '\t' + 't'
            
            g = 0
            
            for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                if gbk.endswith('gbff'):
                    g = g + 1
                    
                    basename = re.split('gbff', gbk)[0]
                    subprocess.call('cd ' + ref_dir_domain + 'refseq/' + d + ';cp ' + gbk + ' ' + basename + 'gbk', shell = True, executable = executable)
                    
                    print >> genetic_elements, 'ID' + '\t' + d + '.' + str(g)
                    print >> genetic_elements, 'NAME' + '\t' + d + '.' + str(g)
                    print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                    print >> genetic_elements, 'CIRCULAR?' + '\t' + 'Y'
                    print >> genetic_elements, 'ANNOT-FILE' + '\t' + basename + 'gbk'
                    print >> genetic_elements, '//'
                    
            if g > 0:
                new_pgdbs.append(d)
                
## Previous failed builds confuse pathway-tools.  Remove the directory.

#### !!!!! remove when metadata saving disable is fixed
try:
    os.remove(pgdb_dir + 'PGDB-METADATA.ocelot')
    os.remove(pgdb_dir + 'PGDB-METADATA.ocelot~')
except OSError:
    pass
                
for d in new_pgdbs:
    shutil.rmtree(pgdb_dir + d.lower() + 'cyc', ignore_errors = True)

## Switched to using Parallel to reduce the number of dependencies.
##### !!!!!! change back to -1 when metadata issue is fixed

print len(new_pgdbs), 'new pgdbs will be created'

if __name__ == '__main__':  
    Parallel(n_jobs = 1, verbose = 5)(delayed(make_pgdb)
    (d, ref_dir_domain) for d in new_pgdbs)

#%% For each PGDB add the pathways to a new data_frame.

terminal_paths = pd.DataFrame(index = assemblies)

for d in assemblies:
    np = 0 # Number of pathways predicted.
    try:
        with open(pgdb_dir + d.lower() + 'cyc/1.0/reports/pathways-report.txt', 'r') as report:            
            for line in report:
                if line.startswith('#') == False:
                    if line.startswith('Pathway Name') == False:
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
                                    print 'collecting paths for terminal node', d, path
                                    np = np + 1
                                    
        genome_data.loc[d, 'npaths_actual'] = np
        
    except IOError:
        print d, 'has no pathway report'

#%% Collect pathway and other data for internal nodes.
## Make an initial pass over the tree to collect all the internal node numbers.

int_nodes = set()

for clade in tree.get_nonterminals():
    edge = int(clade.confidence)
    int_nodes.add(edge)
    
int_nodes = sorted(int_nodes)

## Create a new dataframes to store pathway data, the fraction of daughters
## with each pathway, and genome data inferred from daughters for each internal
## node.

internal_probs = pd.DataFrame(index = int_nodes, columns = terminal_paths.columns)
internal_data = pd.DataFrame(index = int_nodes, columns = ['n16S', 'nge', 'ncds', 'genome_size', 'GC', 'phi', 'clade_size', 'npaths_terminal', 'branch_length'])

for clade in tree.get_nonterminals():
    edge = int(clade.confidence)
    internal_data.loc[edge, 'branch_length'] = clade.branch_length
    temp_paths = pd.DataFrame(columns = terminal_paths.columns)
    
    n16S = []
    nge = []
    ncds = []
    genome_size = []
    phi = []
    npaths_terminal = []
    gc = []
    
    ntip = 0
    
    for tip in clade.get_terminals():
        
        ntip += 1
        name = tip.name
        assembly = genome_data[genome_data['tip_name'] == tip.name].index.tolist()[0]
        print 'collecting data for internal node', str(edge), 'member', assembly
        
        n16S.append(genome_data.loc[assembly, 'n16S'])
        nge.append(genome_data.loc[assembly, 'nge'])
        ncds.append(genome_data.loc[assembly, 'ncds'])
        genome_size.append(genome_data.loc[assembly, 'genome_size'])
        phi.append(genome_data.loc[assembly, 'phi'])
        gc.append(genome_data.loc[assembly, 'GC'])
        
        npaths_terminal.append(genome_data.loc[assembly, 'npaths_actual'])
        
        temp_paths = temp_paths.append(terminal_paths.loc[assembly,:])
        
    npaths = temp_paths.count(axis = 0, numeric_only = True)
    rpaths = npaths.div(ntip)
    internal_probs.loc[edge, :] = rpaths
    
    ## If any nas sneak in they will propgate through calculations, remove.
    
    n16S = pd.Series(n16S).dropna()
    nge = pd.Series(nge).dropna()
    ncds = pd.Series(nge).dropna()
    genome_size = pd.Series(genome_size).dropna()
    phi = pd.Series(phi).dropna()
    npaths_terminal = pd.Series(npaths_terminal).dropna()
    npaths_gc = pd.Series(gc).dropna()
    
    ## Calculate values for this edge.
            
    internal_data.loc[edge, 'n16S'] = n16S.mean()
    internal_data.loc[edge, 'nge'] = nge.mean()
    internal_data.loc[edge, 'ncds'] = ncds.mean()
    internal_data.loc[edge, 'genome_size'] = genome_size.mean()
    internal_data.loc[edge, 'phi'] = phi.mean()
    internal_data.loc[edge, 'clade_size'] = ntip
    internal_data.loc[edge, 'npaths_terminal'] = npaths_terminal.mean()
    internal_data.loc[edge, 'GC'] = npaths_gc.mean()
    
## Write out ya database files.
        
genome_data.to_csv(ref_dir_domain + 'genome_data.final.csv')
terminal_paths.to_csv(ref_dir_domain + 'terminal_paths.csv')
internal_probs.to_csv(ref_dir_domain + 'internal_probs.csv')
internal_data.to_csv(ref_dir_domain + 'internal_data.csv')