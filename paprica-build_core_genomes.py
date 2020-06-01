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
import urllib.request, urllib.error, urllib.parse

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
    domain = 'archaea'
    tree_file = 'test.' + domain + '.combined_16S.' + domain + '.tax.clean.unique.align.phyloxml'
    ref_dir = 'ref_genome_database'
    pgdb_dir = '/volumes/hd2/ptools-local/pgdbs/user/'
    cpus = 36
    
else:        
    domain = command_args['domain']
    tree_file = command_args['tree']
    ref_dir = command_args['ref_dir']
    pgdb_dir = command_args['pgdb_dir']
    cpus = int(command_args['cpus'])
    
## Expand tilde manually.
    
pgdb_dir = os.path.expanduser(pgdb_dir)
    
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'

paprica_path = os.path.dirname(os.path.realpath("__file__")) + '/' # The location of the actual paprica scripts.   
ref_dir_domain = paprica_path + ref_dir + domain + '/'

#%% Define some functions.

## Define a function to use pathway-tools to predict the metabolic pathways for
## each genome downloaded from Genbank, or transcriptome downloaded from MMTSP.
## The PGDBs are now named by assembly so that they can be re-used for each new
## version of the database.

def make_pgdb(d, ref_dir_domain):
    
    print(d, 'start prediction')
    
    predict_pathways = subprocess.Popen('pathway-tools \
    -lisp \
    -no-cel-overview \
    -patho ' + ref_dir_domain + 'refseq/' + d + '/ \
    -disable-metadata-saving \
    &> ' + ref_dir_domain + 'pathos_' + d + '.log', shell = True, executable = executable)
    
    predict_pathways.communicate()   
    print(d, 'prediction complete')
    
## The directory structure for the domain eukaryote does not come automatically
## with Genbank format files.  Define a function that uses the pep.fa and
## swissprot.gff3 files to create a Genbank format file that can be read by
## pathologic during pgdb creation.
    
def create_euk_files(d):
    
    ## First create a df mapping protein id to SwissProt accession number.
    
    columns = ['prot_id', 'swissprot', 'description']
    spt = pd.read_csv(ref_dir_domain + 'refseq/' + d + '/swissprot.gff3', index_col = 0, comment = '#', names = columns, usecols = [0,8,10], sep = ';|\t', engine = 'python')
    spt['swissprot'] = spt['swissprot'].str.replace('Name=Swiss-Prot:', '')
    spt['description'] = spt['description'].str.replace('Description=Swiss-Prot:', '')
    
    ## Create empty list to hold gene features pulled from cds.fa.
    
    features = []
    
    ## Artificial start, stops are needed.
    
    combined_length = 1
    
    print('generating genbank format files for', d + '...')
    
    ## Some directory names differ from the accession number.  Rename these
    ## directories to match the accession number.
    
    for f in os.listdir(ref_dir_domain + 'refseq/' + d):
        if f.endswith('.pep.fa'):
            a = f.split('.pep.fa')[0]
            
    if a != d:
        os.rename(ref_dir_domain + 'refseq/' + d, ref_dir_domain + 'refseq/' + a)
        print('directory', d, 'is now', a)
    
    for record in SeqIO.parse(ref_dir_domain + 'refseq/' + a + '/' + a + '.pep.fa', 'fasta'):
            
        ## The swissprot annotations are indexed by MMETSP record locator, not
        ## by the actual record.id.
        
        sprot_name = str(record.description).split('NCGR_PEP_ID=')[1]
        sprot_name = sprot_name.split(' /')[0]
        	
        try:
            temp_spt = spt.loc[sprot_name, 'swissprot']
        except KeyError:
            continue
        		
        temp_sprot = sprot_df[sprot_df.index.isin(list(temp_spt))]
        	
        ecs = list(set(temp_sprot.ec))
        descriptions = list(set(temp_sprot.name))
        	
        ## Embed all information necessary to create the Genbank file as qualifiers, then
        ## append to this list of records for that genome.
        	
        qualifiers = {'protein_id':sprot_name, 'locus_tag':str(record.id), 'EC_number':ecs, 'product':descriptions, 'translation':str(record.seq)}
        new_feature = SeqFeature.SeqFeature(type = 'CDS', qualifiers = qualifiers)
        new_feature.location = SeqFeature.FeatureLocation(combined_length, combined_length + len(str(record.seq)))
        features.append(new_feature)
        
        combined_length = combined_length + len(str(record.seq))
        
    ## Write the records in Genbank format.  Even though you will ultimately
    ## want to use the gbk extension, to match the (silly) Genbank convention
    ## use gbff.
        
    new_record = SeqRecord(Seq('nnnn', alphabet = IUPAC.ambiguous_dna), id = a, name = a, features = features)   
    SeqIO.write(new_record, open(ref_dir_domain + 'refseq/' + a + '/' + a + '.gbff', 'w'), 'genbank')
    
#%% Preparatory file generation and organization.

## Read in the genome_data file.

genome_data = pd.read_csv(ref_dir_domain + 'genome_data.csv.gz', header = 0, index_col = 0)
genome_data['clade'] = np.nan
genome_data['tip_name'] = np.nan
genome_data['npaths_actual'] = np.nan
genome_data['branch_length'] = np.nan
    
## Get the clade number of each assembly and add this information to 
## genome_data.
    
tree = Phylo.read(tree_file, 'phyloxml')
    
assemblies = []
    
for clade in tree.get_terminals():
    clade_number = int(clade.confidence)
    print(clade_number)
    
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
    print('Downloading enzyme.dat from ftp.expasy.org...')
    enzyme = urllib.request.urlopen('ftp://ftp.expasy.org/databases/enzyme/enzyme.dat').read()
    
    wget0 = subprocess.Popen('cd ' + ref_dir_domain + ';wget ftp://ftp.expasy.org/databases/enzyme/enzyme.dat', shell = True, executable = executable)
    wget0.communicate()
    
    print('Parsing enzyme.dat to enzyme_table.dat...')
    with open(ref_dir_domain + 'enzyme.dat', 'r') as enzyme, open('enzyme_table.dat', 'w') as enzyme_out:
        for line in enzyme:
            print('accession', 'ec', 'name', file=enzyme_out)
            
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
                        
                        print(sprot, ec, name, file=enzyme_out)

    print('Reading enzyme_table.dat...')
    sprot_df = pd.read_csv('enzyme_table.dat', header = 0, index_col = 0, sep = ' ')

    ## For eukarya, generate gbff files from pep.fa.  This takes a long time, so
    ## only do this for directories that don't already have this file.  This should
    ## only happen if something went wrong on a previous build, or you deleted some
    ## or all of the Genbank files.
    
    ## !!! For eukarya paprica-make_ref.py always deletes the database and
    ## !!! downloads everything from scratch, the gbff files will always need
    ## !!! to be rebuilt.  This needs to be fixed.
    
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
        (d) for d in gbff_needed)
    
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
        
        ## As of ptools v24 pathway-report.txt reformatted with date.  Need
        ## to find name of this file.
        
        report_file = 'none'
        
        for f in os.listdir(pgdb_dir + d.lower() + 'cyc/1.0/reports'):
            
            ## If an old version of pathways-report.txt is present remove it.            
            
            if f == 'pathways-report.txt':
                os.remove(pgdb_dir + d.lower() + 'cyc/1.0/reports/' + f)
                report_file = 'none'
                
            elif f.startswith('pathways-report_'):
                report_file = f
                
        ## If there was no pathways-report file, or it was the old format and deleted,
        ## rewrite the files needed by pathologic.
                
        if report_file == 'none':
            
            clade = genome_data.loc[d, 'clade']
                
            with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
                
                print('recreating pathologic files for', d)
                
                print('ID' + '\t' + d, file=organism_params)
                print('Storage' + '\t' + 'File', file=organism_params)
                print('Name' + '\t' + d, file=organism_params)
                print('Rank' + '\t' + 'Strain', file=organism_params)
                print('Domain' + '\t' + 'TAX-2', file=organism_params)
                print('Create?' + '\t' + 't', file=organism_params)
                
                g = 0
                
                for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                    if gbk.endswith('gbff'):
                        g = g + 1
                        
                        basename = re.split('gbff', gbk)[0]
                        subprocess.call('cd ' + ref_dir_domain + 'refseq/' + d + ';cp ' + gbk + ' ' + basename + 'gbk', shell = True, executable = executable)
                        
                        print('ID' + '\t' + d + '.' + str(g), file=genetic_elements)
                        print('NAME' + '\t' + d + '.' + str(g), file=genetic_elements)
                        print('TYPE' + '\t' + ':CHRSM', file=genetic_elements)
                        
                        if domain == 'eukarya':
                            print('CIRCULAR?' + '\t' + 'Y', file=genetic_elements)
                        else:
                            print('CIRCULAR?' + '\t' + 'N', file=genetic_elements)
                            
                        print('ANNOT-FILE' + '\t' + basename + 'gbk', file=genetic_elements)
                        print('//', file=genetic_elements)
                        
            if g > 0:
                new_pgdbs.append(d)
                
    ## If there was no previous build attempt the directory will not exist and
    ## os will throw an error.
            
    except FileNotFoundError:
        
        clade = genome_data.loc[d, 'clade']
            
        with open(ref_dir_domain + 'refseq/' + d  + '/organism-params.dat', 'w') as organism_params, open(ref_dir_domain + 'refseq/' + d  + '/genetic-elements.dat', 'w') as genetic_elements:                
    
            print('ID' + '\t' + d, file=organism_params)
            print('Storage' + '\t' + 'File', file=organism_params)
            print('Name' + '\t' + d, file=organism_params)
            print('Rank' + '\t' + 'Strain', file=organism_params)
            print('Domain' + '\t' + 'TAX-2', file=organism_params)
            print('Create?' + '\t' + 't', file=organism_params)
            
            g = 0
            
            for gbk in os.listdir(ref_dir_domain + 'refseq/' + d):
                if gbk.endswith('gbff'):
                    g = g + 1
                    
                    basename = re.split('gbff', gbk)[0]
                    subprocess.call('cd ' + ref_dir_domain + 'refseq/' + d + ';cp ' + gbk + ' ' + basename + 'gbk', shell = True, executable = executable)
                    
                    print('ID' + '\t' + d + '.' + str(g), file=genetic_elements)
                    print('NAME' + '\t' + d + '.' + str(g), file=genetic_elements)
                    print('TYPE' + '\t' + ':CHRSM', file=genetic_elements)
                    print('CIRCULAR?' + '\t' + 'Y', file=genetic_elements)
                    print('ANNOT-FILE' + '\t' + basename + 'gbk', file=genetic_elements)
                    print('//', file=genetic_elements)
                    
            if g > 0:
                new_pgdbs.append(d)
                
## Previous failed builds confuse pathway-tools.  Remove the directory.

print('removing previous failed PGDBs...')
                
for d in new_pgdbs:
    shutil.rmtree(pgdb_dir + d.lower() + 'cyc', ignore_errors = True)
    
#%% Generate the PGDBs for each assembly that does not have one.

## Previously this was done externally with Gnu Parallel. Switched to using joblib
## to reduce the number of dependencies.

print(len(new_pgdbs), 'new pgdbs will be created')

if __name__ == '__main__':  
    Parallel(n_jobs = cpus, verbose = 5)(delayed(make_pgdb)
    (d, ref_dir_domain) for d in new_pgdbs)
    
#%% For each PGDB add the pathways to a new data_frame.

terminal_paths = pd.DataFrame(index = assemblies)

for i,d in enumerate(assemblies):
    n_paths = 0 # Number of pathways predicted.
    try:
        
        ## As of ptools v24 pathway-report.txt reformatted with date.  Need
        ## to find name of this file.
        
        for f in os.listdir(pgdb_dir + d.lower() + 'cyc/1.0/reports'):
            if f.startswith('pathways-report'):
                report_file = f
                
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
        
    except IOError:
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

int_nodes = set()

for clade in tree.get_nonterminals():
    try:
        edge = int(clade.confidence)
        int_nodes.add(edge)
    except TypeError:
        continue
    
int_nodes = sorted(int_nodes)
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

## Define a funciton to iterate across all subtrees and collect information that will be saved in
## the "internal" memory-mapped arrays. Remember that these are not dataframes
## and cannot be indexed by column/row names!

def get_internals(clade,
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

    try:
        int(clade.confidence)
    except TypeError:
        return
           
    if clade.confidence > 0:
        
        edge = int(clade.confidence)
        print('collecting data for internal node', str(edge))
        
        ## Data on the clade that you want later.
        
        edge_i = int_nodes.index(edge)
        internal_data[edge_i, internal_data_columns.index('branch_length')] = clade.branch_length
        
        ## Iterate across all terminal nodes in clades to get the corresponding
        ## assemblies.
        
        ntip = len(clade.get_terminals())
        clade_members = []
        
        for tip in clade.get_terminals():

            assembly = genome_data[genome_data['tip_name'] == tip.name].index.tolist()[0]        
            clade_members.append(assembly)
            
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
        
        if domain != 'eukarya':
            
            ## These parameters cannot be calculated from transcriptomes.
            
            internal_data[edge_i, internal_data_columns.index('n16S')] = clade_data['n16S'].dropna().mean()
            internal_data[edge_i, internal_data_columns.index('nge')] = clade_data['nge'].dropna().mean()
            internal_data[edge_i, internal_data_columns.index('ncds')] = clade_data['ncds'].dropna().mean()
            internal_data[edge_i, internal_data_columns.index('genome_size')] = clade_data['genome_size'].dropna().mean()
            internal_data[edge_i, internal_data_columns.index('phi')] = clade_data['phi'].dropna().mean()
            internal_data[edge_i, internal_data_columns.index('GC')] = clade_data['GC'].dropna().mean()
            
        ## These parameters are relevant to all domains.
        
        internal_data[edge_i, internal_data_columns.index('clade_size')] = ntip
        internal_data[edge_i, internal_data_columns.index('npaths_terminal')] = clade_data['npaths_actual'].dropna().mean()
        internal_data[edge_i, internal_data_columns.index('nec_terminal')] = clade_data['nec_actual'].dropna().mean()
        
#%% Execute the get_internals function in parallel, this massively speeds up
## the database build.
        
## For bacteria currently throws RuntimeError: maximum recursion depth exceeded
## if more than 1 processor used.  This is kind of a problem because the
## bacteria database is the one that takes a long time to build...
        
if domain == 'bacteria':
    njobs = 1
else:
    njobs = -1

if __name__ == '__main__':         
    Parallel(n_jobs = njobs)(delayed(get_internals)(clade, int_nodes, genome_data, terminal_paths, terminal_ec, internal_probs_columns, internal_data_columns, internal_ec_probs_columns, internal_ec_n_columns, internal_probs, internal_data, internal_ec_probs, internal_ec_n)
    for clade in tree.get_nonterminals())
        
#%% Collect taxonomy information for each of the nodes in the reference tree.

lineage = pd.read_csv(ref_dir_domain + 'taxa.csv', index_col = 0)
ref_taxa = pd.read_csv(ref_dir_domain + 'seq_info.updated.csv', index_col = 0)

ranks = lineage.columns

node_lineages = pd.DataFrame(columns = lineage.columns)
node_lineages_index = []
    
for clade in tree.get_nonterminals():
    try:
        clade_number = int(clade.confidence)
        print('getting lineage for', clade_number)
        terminals = []
        
        for terminal in clade.get_terminals():
            terminals.append(terminal.name.strip('@'))
            
        temp_taxids = ref_taxa.loc[terminals, 'tax_id']
        
        temp_lineage = lineage.loc[temp_taxids]
        temp_lineage = temp_lineage.dropna(1, thresh = 1)
        temp_lineage.drop(['parent_id', 'rank', 'tax_name'], inplace = True, axis = 1)
        
        ## Now iterate across columns, starting at root
        ## until you find the first mismatch.  The one
        ## before this is the consensus.
        
        for i,rank in enumerate(temp_lineage.columns):
            if i != 0:
                if len(temp_lineage[rank].unique()) > 1:
                    consensus_rank = temp_lineage.columns[i - 1]
                    consensus_taxid = temp_lineage.loc[temp_lineage.index[0], consensus_rank]
                    break
                    
        ## Now look up the consensus taxonomy.
            
        consensus_taxa = lineage.loc[consensus_taxid, 'tax_name']
        consensus_lineage = lineage.loc[consensus_taxid]
        
        if len(consensus_lineage.shape) > 1:
            consensus_lineage = consensus_lineage.drop_duplicates()
        
        ## Save consensus lineage.
        
        node_lineages = node_lineages.append(consensus_lineage)
        node_lineages_index.append(clade_number)

    ## Error exception here is problematic, for some reason KeyError
    ## exception does not catch error raised by temp_taxids not being present
    ## in temp_lineage.  Works without being explicit on error type, but this
    ## is bad.
            
    except:
        print('none')
        
for clade in tree.get_terminals():
    
    try:        
        clade_number = int(clade.confidence)
        print('getting lineage for', clade_number)
        terminal = clade.name.strip('@')
    
        temp_taxid = ref_taxa.loc[terminal, 'tax_id']
        temp_lineage = lineage.loc[temp_taxid]
        node_lineages = node_lineages.append(temp_lineage)
        node_lineages_index.append(clade_number)
    
    except:
        print('none')
           
node_lineages.index = node_lineages_index
    
for rank in ranks:
    for index in node_lineages.index:
        temp = node_lineages.loc[index, rank]
        
        if temp in lineage.index:
            node_lineages.loc[index, rank] = lineage.loc[temp, 'tax_name']
            
## Check that all edges have a minimal entry in lineages, this is important
## as the euks get a bit whonky and sometimes there's no tax info.
            
for edge in int_nodes:
    if edge not in node_lineages.index:
        node_lineages.loc[edge, 'parent_id'] = 1
        
for edge in terminal_paths.index:
    if edge not in genome_data.clade:
        node_lineages.loc[edge, 'parent_id'] = 1
   
## Write out ya database files.
            
node_lineages.to_csv(ref_dir_domain + 'node_lineages.csv.gz')      
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