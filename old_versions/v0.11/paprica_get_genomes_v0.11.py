# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 11:05:04 2015

@author: jeff

"""

### user setable variables ###

ref_dir = '/home/user/genome_finder/ref_genome_database_a/' # location of genome database created with paprika_build_core_genomes.py
pgdb_dir = '/home/user/ptools-local/pgdbs/user/' # location of pathway-tools pgdbs
strain_dir = '/home/user/ref_genome_database_a/ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' # location of strain genomes 
version = 'a' # change version number if you are rebuilding the database

### end user setable variables ###

import sys
import os
from Bio import Phylo
import numpy as np

ref_package = 'combined_16S' # name of reference package

name = sys.argv[1]
name_out = sys.argv[2]

## diagnostic only !
#name = 'PAL_219_20131207_F02.combined_16S.pplacer.filter.jplace.csv'
#name_out = 'PAL_219_20131207_F02'

wd = os.getcwd() + '/'
internal_nodes = set(os.listdir(ref_dir + ref_package + '.core_genomes'))

## find all terminal node uids

ref_edges = {}

tree_name = name.rstrip('csv')
tree_name = tree_name + 'fat'
tree = Phylo.read(tree_name, 'phyloxml')

for tip in tree.get_terminals():
    edge =  str(int(tip.confidence.value))
    tip_name = tip.name
    tip_name = tip_name.strip('@ref_')
    uid = tip_name[tip_name.find('uid'):]
    ref_edges[edge] = uid

## generate a dictionary mapping uid to strain for all reference strains

ref_uid_strain = {}

for ref in os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria'):
    ref_uid = ref[ref.find('uid'):]
    ref_uid_strain[ref_uid] = ref

#### collect statistics on all reference genomes ####

## map phi to uid

ref_uid_strain = {}
phi_dict = {}
min_phi = 1
n16S_dict = {}

with open(ref_dir + 'genome_data.txt', 'r') as phi_file:
    for line in phi_file:
        line = line.rstrip()
        line = line.split()
        
        line_uid = line[0]
        line_strain = line[1]
        ref_uid_strain[line_uid] = line_strain
                
        try:
            if float(line[2]) < min_phi:
                min_phi = float(line[2])
                phi_dict[line_uid] = line[2]
            else:
                phi_dict[line_uid] = line[2]
        except ValueError:
            continue
            
        n16S_dict[line_uid] = line[3]

## for core genomes, collect statistics for each node number

stats_dict = {}

for d in internal_nodes:
    with open(ref_dir + ref_package + '.core_genomes' + '/' + d + '/' + d + '.stats', 'r') as stats:
        temp_size = []
        temp_phi = []
        temp_16S = []
        
        for line in stats:
            line = line.rstrip()
            line = line.split()
            
            if line[0] != 'core':
                line_uid = line[0].rstrip('_combined.fasta')
                temp_size.append(line[1])
                
                ## there are some strains for which phi cannot be calculated, use lowest value for these                
                
                try:
                    temp_phi.append(phi_dict[line_uid])
                except KeyError:
                    temp_phi.append(min_phi)
                    
                try:
                    temp_16S.append(n16S_dict[line_uid])
                except KeyError:
                    continue
                
            else:
                temp_phi = map(float, temp_phi)
                temp_phi = np.array(temp_phi)
                temp_mean_phi = np.mean(temp_phi)
                
                temp_size = map(float, temp_size)
                temp_size = np.array(temp_size)
                temp_mean = np.mean(temp_size)
                temp_sd = np.std(temp_size)
                core_size = line[1]
                
                temp_16S = map(int, temp_16S)
                temp_16S = np.array(temp_16S)
                temp_mean_16S = np.mean(temp_16S)
                
                stats_dict[d] = core_size, temp_mean, temp_sd, temp_mean_phi, temp_mean_16S

## tally edge placements

query_edges = {}
edge_tally = {}

read = False # skip the first line

with open(name, 'r') as csv_in:
    for line in csv_in:
        if read != False:
            line = line.rstrip()
            line = line.split(',')
            seq = line[1]
            if seq.startswith('ref_'):
                pass
            else:
                query_edges[line[3]] = seq
                
                ## count the number of times each edge appears in sample                
                
                try:
                    temp = edge_tally[line[3]]
                    temp = temp + 1
                    edge_tally[line[3]] = temp
                
                except KeyError:
                    edge_tally[line[3]] = 1

        read = True
        
with open(wd + name_out + '.edge_tally.txt', 'w') as tally_out:
    sample_score = []
    opt_sample_score = []
    print >> tally_out, 'edge' + '\t' + 'nplacements' + '\t' + 'nplacements.corrected' + '\t' + 'core.size' + '\t' + 'mean.size' + '\t' + 'sd.size ' + '\t' + '1-phi' + '\t' + '16S.copy'
    
    for key in edge_tally:
        nplace = edge_tally[key] # number of placements to edge
        
        try:
            edge_stats = stats_dict[str(key)]
            n16S = edge_stats[4]
            ncplace = round(nplace / n16S)
            
            ## it makes no sense to have less than one observed 16S, so round up to 1
            if ncplace < 1:
                ncplace = 1
                
            opt_sample_score.append(nplace) # optimum score is number of placements * 1
            
            ## sample score = nplacements * (core genome size / mean cluster genome size) * 1 - meanphi
            sample_score.append(ncplace * (float(edge_stats[0]) / float(edge_stats[1])) * (1 - float(edge_stats[3])))

        ## terminal nodes will not have an entry in stat_dict
        except KeyError:
            try:            
                edge_uid = ref_edges[key]
                edge_phi = phi_dict[edge_uid]
                n16S = float(n16S_dict[edge_uid])
                ncplace = round(nplace / n16S)
                
                ## it makes no sense to have less than one observed 16S, so round up to 1
                if ncplace < 1:
                    ncplace = 1
                    
                edge_stats = 'NA', 'NA', 'NA', edge_phi, n16S
                opt_sample_score.append(ncplace) # optimum score is number of placements * 1

            except KeyError:
                edge_stats = 'NA', 'NA', 'NA', min_phi, 'NA'
                ncplace = nplace
        
            ## in this case edge score is just nplacements * (phi)
            sample_score.append(ncplace * (1 - float(edge_stats[3])))
            edge_name = query_edges[key]
           
        print >> tally_out, str(key) + '\t' + str(nplace) + '\t' + str(ncplace) + '\t' + str(edge_stats[0]) + '\t' + str(edge_stats[1]) + '\t' + str(edge_stats[2]) + '\t' + str(1 - float(edge_stats[3])) + '\t' + str(edge_stats[4])

    try:
        sample_score = sum(sample_score) / sum(opt_sample_score)
    except ZeroDivisionError:
        sample_score = 'NA'
        
    print >> tally_out, 'sample.score' + '\t' + str(sample_score)
    
pgdbs = set(os.listdir(pgdb_dir))

## setup core and strain genome directories for pgdb creation, if needed 

with open('generate_pgdbs.sh', 'w') as run_pgdb:
    print >> run_pgdb, '#!/bin/bash'
    
    for key in edge_tally.keys():
        if version + key + 'cyc' not in pgdbs:
            if key in internal_nodes:
                
                with open(ref_dir + ref_package + '.core_genomes/' + key + '/' + 'organism-params.dat', 'w') as organism_params:
                    print >> organism_params, 'ID' + '\t' + version + key
                    print >> organism_params, 'Storage' + '\t' + 'File'
                    print >> organism_params, 'Name' + '\t' + version + key
                    print >> organism_params, 'Rank' + '\t' + 'Strain'
                    print >> organism_params, 'Domain' + '\t' + 'TAX-2'
                    print >> organism_params, 'Create?' + '\t' + 't'
                    
                with open(ref_dir + ref_package + '.core_genomes/' + key + '/' + 'genetic-elements.dat', 'w') as genetic_elements:            
                    print >> genetic_elements, 'ID' + '\t' + version + key + '.1'
                    print >> genetic_elements, 'NAME' + '\t' + version + key + '.1'
                    print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                    print >> genetic_elements, 'CIRCULAR?' + '\t' + 'N'
                    print >> genetic_elements, 'ANNOT-FILE' + '\t' + key + '_core_genome.gbk'
                    print >> genetic_elements, '//'
                    
                print >> run_pgdb, 'pathway-tools -lisp -no-cel-overview -patho ' + ref_dir + ref_package + '.core_genomes/' + key + '/ -disable-metadata-saving &> pathos_' + version + key + '.log'   
                    
            else:
                uid = ref_edges[key]
                strain = ref_uid_strain[uid]
                
                with open(strain_dir + strain + '/' + 'organism-params.dat', 'w') as organism_params, open(strain_dir + strain + '/' + 'genetic-elements.dat', 'w') as genetic_elements:                

                    print >> organism_params, 'ID' + '\t' + version + key
                    print >> organism_params, 'Storage' + '\t' + 'File'
                    print >> organism_params, 'Name' + '\t' + version + key
                    print >> organism_params, 'Rank' + '\t' + 'Strain'
                    print >> organism_params, 'Domain' + '\t' + 'TAX-2'
                    print >> organism_params, 'Create?' + '\t' + 't'
                    
                    g = 0
                    for gbk in os.listdir(strain_dir + strain):
                        if gbk.endswith('gbk'):
                            g = g + 1
                            
                            print >> genetic_elements, 'ID' + '\t' + version + key + '.' + str(g)
                            print >> genetic_elements, 'NAME' + '\t' + version + key + '.' + str(g)
                            print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                            print >> genetic_elements, 'CIRCULAR?' + '\t' + 'N'
                            print >> genetic_elements, 'ANNOT-FILE' + '\t' + gbk
                            print >> genetic_elements, '//'
                    
                print >> run_pgdb, 'pathway-tools -lisp -no-cel-overview -patho ' + strain_dir + strain + '/ -disable-metadata-saving &> pathos_' + version + key + '.log'   
                
  