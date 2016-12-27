#!/usr/bin/env python
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

This script must be located in the 'paprica' directory as it makes use of relative
paths.

"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqFeature
from Bio import Phylo

import subprocess
import os
import re
import sys
import shutil
from joblib import Parallel, delayed
import urllib2

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
        
## Define some variables based on these arguments.  If nothing in sys.argv
## set some default values.  This is useful for testing.
        
if len(sys.argv) == 1:
    domain = 'eukarya'
    tree = 'test.eukarya.combined_18S.eukarya.tax.clean.align.phyloxml'
    ref_dir = 'ref_genome_database'
    pgdb_dir = '~/ptools-local/pgdbs/user/'
    
else:        
    domain = command_args['domain']
    tree = command_args['tree']
    ref_dir = command_args['ref_dir']
    pgdb_dir = command_args['pgdb_dir']
    
## Expand tilde manually.
    
pgdb_dir = os.path.expanduser(pgdb_dir)
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'

paprica_path = os.path.dirname(os.path.abspath("__file__")) + '/' # The location of the actual paprica scripts.   
ref_dir_domain = paprica_path + ref_dir + domain + '/'

#%% Define some functions.

## Define a function to use pathway-tools to predict the metabolic pathways for
## each genome downloaded from Genbank, or transcriptome downloaded from MMTSP.
## The PGDBs are now named by assembly so that they can be re-used for each new
## version of the database.

def make_pgdb(d, ref_dir_domain):
    
    print d, 'start prediction'
    predict_pathways = subprocess.Popen('pathway-tools -lisp -no-cel-overview -patho ' + ref_dir_domain + 'refseq/' + d + '/ -disable-metadata-saving &> ' + ref_dir_domain + 'pathos_' + d + '.log', shell = True, executable = executable)
    predict_pathways.communicate()   
    print d, 'prediction complete'
    
## The directory structure for the domain eukaryote does not come automatically
## with Genbank format files.  Define a function that uses the pep.fa and
## swissprot.gff3 files to create a Genbank format file that can be read by
## pathologic during pgdb creation.
    
def create_euk_files(d, ref_dir_domain):
    
    ## First create a df mapping protein id to SwissProt accession number.
    
    columns = ['prot_id', 'swissprot', 'description']
    spt = pd.read_csv(ref_dir_domain + 'refseq/' + d + '/swissprot.gff3', index_col = 0, comment = '#', names = columns, usecols = [0,8,10], sep = ';|\t', engine = 'python')
    spt['swissprot'] = spt['swissprot'].str.replace('Name=Swiss-Prot:', '')
    spt['description'] = spt['description'].str.replace('Description=Swiss-Prot:', '')
    
    ## Create empty list to hold gene features pulled from cds.fa.
    
    features = []
    
    ## Artificial start, stops are needed.
    
    combined_length = 1
    
    print 'generating genbank format files for', d + '...'
    
    for record in SeqIO.parse(ref_dir_domain + 'refseq/' + d + '/' + d + '.pep.fa', 'fasta'):
            
        new_name = str(record.description).split('NCGR_PEP_ID=')[1]
        new_name = new_name.split(' /')[0]
        record.id = new_name
        	
        try:
            temp_spt = spt.loc[record.id, 'swissprot']
        except KeyError:
            continue
        		
        temp_sprot = sprot_df[sprot_df.index.isin(list(temp_spt))]
        	
        ecs = list(set(temp_sprot.ec))
        descriptions = list(set(temp_sprot.name))
        	
        ## Embed all information necessary to create the Genbank file as qualifiers, then
        ## append to this list of records for that genome.
        	
        qualifiers = {'locus_tag':str(record.id), 'EC_number':ecs, 'product':descriptions, 'translation':str(record.seq)}
        new_feature = SeqFeature.SeqFeature(type = 'CDS', qualifiers = qualifiers)
        new_feature.location = SeqFeature.FeatureLocation(combined_length, combined_length + len(str(record.seq)))
        features.append(new_feature)
        
        combined_length = combined_length + len(str(record.seq))
        
    ## Write the records in Genbank format.  Even though you will ultimately
    ## want to use the gbk extension, to match the (silly) Genbank convention
    ## use gbff.
        
    new_record = SeqRecord(Seq('nnnn', alphabet = IUPAC.ambiguous_dna), id = d, name = d, features = features)   
    SeqIO.write(new_record, open(ref_dir_domain + 'refseq/' + d + '/' + d + '.gbff', 'w'), 'genbank')
    
#%% Preparatory file generation and organization.

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
    
    if domain == 'eukarya':
        assembly = re.split('_', assembly)[0]

    else:        
        assembly = re.split('_', assembly)
        assembly = assembly[0] + '_' + assembly[1]
    
    genome_data.loc[assembly, 'clade'] = clade_number
    genome_data.loc[assembly, 'tip_name'] = clade.name
    genome_data.loc[assembly, 'branch_length'] = clade.branch_length
    
    assemblies.append(assembly)
    
## For eukaryotes, the EC number associated with reads from the MMETSP database
## comes from mapping the swissprot hits for each read to the annotation file
## for the swissprot database.  To conduct this mapping, if domain == eukaryotes
## download and parse enzyme.dat from ftp://ftp.expasy.org/databases/enzyme/enzyme.dat.
    
if domain == 'eukarya':
    print 'Downloading enzyme.dat from ftp.expasy.org...'
    enzyme = urllib2.urlopen('ftp://ftp.expasy.org/databases/enzyme/enzyme.dat')
    
    print 'Parsing enzyme.dat to enzyme_table.dat...'
    with open('enzyme_table.dat', 'w') as enzyme_out:
        for line in enzyme:
            print >> enzyme_out, 'accession', 'ec', 'name'
            
            for line in enzyme:
                if line.startswith('ID'):
                    ec = line.split()[1]
                    ec = ec.rstrip()
                    
                if line.startswith('DE'):
                    name = line.split()[1]
                    name = name.rstrip()
                    
                if line.startswith('DR'):
                    line = line.strip('DR')
                    line = line.strip()
                    line = line.rstrip()
                    line = line.rstrip(';')
                    line = line.split(';')
                    
                    for sprot in line:
                        sprot = sprot.split(',')[0]
                        sprot = sprot.strip()
                        sprot = sprot.rstrip()
                        
                        print >> enzyme_out, sprot, ec, name
                    
    enzyme.close()

    print 'Reading enzyme_table.dat...'
    sprot_df = pd.read_csv('enzyme_table.dat', header = 0, index_col = 0, sep = ' ')

## For eukarya, generate gbff files from pep.fa.  This takes a long time, so
## only do this for directories that don't already have this file.  This should
## only happen if something went wrong on a previous build, or you deleted some
## or all of the Genbank files.

if domain == 'eukarya':
    
    ## Determine which directories need the gbff file.
    
    gbff_needed = []
    
    for d in os.listdir(ref_dir_domain + 'refseq'):
        need = True
        for f in os.listdir(ref_dir_domain + 'refseq/' + d):
            if f.endswith('gbff'):
                need = False
                
        if need == True:
            gbff_needed.append(d)
        
    ## Execute create_euk_files on those directories that need gbff.
        
    if __name__ == '__main__':  
        Parallel(n_jobs = -1, verbose = 5)(delayed(create_euk_files)
        (d, ref_dir_domain) for d in gbff_needed)
    
## For every existing PGDB directory (which might be none), determing if the
## pathways-report.txt file is present.  If it is not this means the previous
## pathway-tools pathologic run was incorrect.  In that case delete the
## directory and try again.
    
new_pgdbs = []

for d in assemblies:
    
    ## If a previous build effort was unsuccessful rewrite the files needed by
    ## pathologic, in case the data files have been updatated in the public
    ## repository (Genbank or MMETSP) and this fixes the problem.    
    
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
                        
                        if domain == 'eukarya':
                            print >> genetic_elements, 'CIRCULAR?' + '\t' + 'Y'
                        else:
                            print >> genetic_elements, 'CIRCULAR?' + '\t' + 'N'
                            
                        print >> genetic_elements, 'ANNOT-FILE' + '\t' + basename + 'gbk'
                        print >> genetic_elements, '//'
                        
            if g > 0:
                new_pgdbs.append(d)
                
    ## If there was no previous build attempt the directory will not exist and
    ## os will throw an error.
            
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
                
for d in new_pgdbs:
    shutil.rmtree(pgdb_dir + d.lower() + 'cyc', ignore_errors = True)
    
#%% Generate the PGDBs for each assembly that does not have one.

## Previously this was done externally with Gnu Parallel. Switched to using joblib
## to reduce the number of dependencies.

print len(new_pgdbs), 'new pgdbs will be created'

if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(make_pgdb)
    (d, ref_dir_domain) for d in new_pgdbs)
    
#%% For each PGDB add the pathways to a new data_frame.

terminal_paths = pd.DataFrame(index = assemblies)

for i,d in enumerate(assemblies):
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
                                        
                                    print 'collecting paths for terminal node', d, i + 1, 'of', len(assemblies), path
                                    np = np + 1
                                    
        genome_data.loc[d, 'npaths_actual'] = np
        
    except IOError:
        print d, 'has no pathway report'
        
#%% Collect EC_numbers for each terminal node

## Read in user specified EC numbers.
        
user_ec = pd.read_csv(ref_dir + 'user/' + 'user_ec.csv', header = 0, index_col = 0, comment = '#')
        
terminal_ec = pd.DataFrame(index = assemblies)
ec_names = pd.DataFrame()
        
for i,d in enumerate(assemblies):
    np = 0
    
    for f in os.listdir(ref_dir_domain + '/refseq/' + d):
        if f.endswith('gbff'):
            for record in SeqIO.parse(ref_dir_domain + '/refseq/' + d + '/' + f, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        if 'EC_number' in feature.qualifiers.keys():
                            
                            np = np + 1
                            ec = feature.qualifiers['EC_number']
                            
                            ## Some draft genomes will not have a product qualifier.
                            
                            try:
                                prod = feature.qualifiers['product'][0]
                            except KeyError:
                                prod = 'product not specified'
                            
                            ## Because each EC number can appear multiple times
                            ## in a genome this information needs to be tallied.
                            
                            for each in ec:
                                print 'collecting EC numbers for terminal node', d, i + 1, 'of', len(assemblies), each
                                
                                try:
                                    temp = terminal_ec.loc[d, each]
                                    if pd.isnull(temp) == True:
                                        terminal_ec.loc[d, each] = 1
                                    else:
                                        terminal_ec.loc[d, each] = temp + 1
                                except KeyError:
                                    terminal_ec.loc[d, each] = 1
                                    
                                ec_names.loc[each, 'name'] = prod
        
    ## For assembly d, add in any user specified EC numbers.
    
    if d in set(user_ec['GI_number']):
        temp_user_ec = user_ec[user_ec['GI_number'] == d]
        
        for entry in temp_user_ec.index:
            np = np + 1
            each = temp_user_ec.loc[entry, 'EC_number']
            print 'collecting EC numbers for terminal node', d, i + 1, 'of', len(assemblies), each
                                
            try:
                if pd.isnull(terminal_ec.loc[d, each]) == True:
                    terminal_ec.loc[d, each] = 1
                else:
                    print 'There is already an entry for', d, each
            except KeyError:
                terminal_ec.loc[d, each] = 1
                
            ## Check to make sure that the enzyme number has a name, add if it does not
                
            try:
                if pd.isnull(ec_names.loc[each, 'name']):
                    ec_names.loc[each, 'name'] = temp_user_ec.loc[entry, 'product']
            except KeyError:
                ec_names.loc[each, 'name'] = temp_user_ec.loc[entry, 'product']

    ## Add the total number of enzymes for that assembly to genome_data.

    genome_data.loc[d, 'nec_actual'] = np

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
internal_ec_probs = pd.DataFrame(index = int_nodes, columns = terminal_ec.columns)
internal_ec_n = pd.DataFrame(index = int_nodes, columns = terminal_ec.columns)

## Determine how many clade members you need to iterate across.

n_clades = len(tree.get_nonterminals())
i_clade = 0

for clade in tree.get_nonterminals():
    
    i_clade += 1
    print 'collecting data for internal node', str(edge) + ',', i_clade, 'of', n_clades
    
    ## Data on the clade that you want later.
    
    edge = int(clade.confidence)
    internal_data.loc[edge, 'branch_length'] = clade.branch_length
    
    ## Iterate across all terminal nodes in clades to get the corresponding
    ## assemblies.
    
    ntip = len(clade.get_terminals())
    clade_members = []
    
    for tip in clade.get_terminals():
        
        name = tip.name
        assembly = genome_data[genome_data['tip_name'] == tip.name].index.tolist()[0]        
        clade_members.append(assembly)
        
    ## Get data for all clade members.
        
    clade_data = genome_data.loc[clade_members]
    clade_paths = terminal_paths.loc[clade_members]
    clade_ec = terminal_ec.loc[clade_members]
        
    npaths = clade_paths.count(axis = 0, numeric_only = True)
    rpaths = npaths.div(ntip)
    internal_probs.loc[edge, :] = rpaths
    
    nec = clade_ec.count(axis = 0, numeric_only = True)
    rec = nec.div(ntip)
    internal_ec_probs.loc[edge, :] = rec
    
    ## The mean number of occurrences of EC_numbers are calculated for clade
    ## so that later on an estimate can be given of the number that will appear
    ## in an internal node.
    
    mec = clade_ec.mean(axis = 0, numeric_only = True)
    internal_ec_n.loc[edge, :] = mec
    
    ## Calculate values for this edge.
    
    if domain != 'eukarya':
        
        ## These parameters cannot be calculated from transcriptomes.
            
        internal_data.loc[edge, 'n16S'] = clade_data['n16S'].dropna().mean()
        internal_data.loc[edge, 'nge'] = clade_data['nge'].dropna().mean()
        internal_data.loc[edge, 'ncds'] = clade_data['ncds'].dropna().mean()
        internal_data.loc[edge, 'genome_size'] = clade_data['genome_size'].dropna().mean()
        internal_data.loc[edge, 'phi'] = clade_data['phi'].dropna().mean()
        internal_data.loc[edge, 'GC'] = clade_data['GC'].dropna().mean()
        
    ## These parameters are relevant to all domains.
        
    internal_data.loc[edge, 'clade_size'] = ntip
    internal_data.loc[edge, 'npaths_terminal'] = clade_data['npaths_actual'].dropna().mean()
    internal_data.loc[edge, 'nec_terminal'] = clade_data['nec_actual'].dropna().mean()
    
## Write out ya database files.
        
genome_data.to_csv(ref_dir_domain + 'genome_data.final.csv')
internal_data.to_csv(ref_dir_domain + 'internal_data.csv')

terminal_paths.to_csv(ref_dir_domain + 'terminal_paths.csv')
internal_probs.to_csv(ref_dir_domain + 'internal_probs.csv')

internal_ec_probs.to_csv(ref_dir_domain + 'internal_ec_probs.csv')
internal_ec_n.to_csv(ref_dir_domain + 'internal_ec_n.csv')
terminal_ec.to_csv(ref_dir_domain + 'terminal_ec.csv')