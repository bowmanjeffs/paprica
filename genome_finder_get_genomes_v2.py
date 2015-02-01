# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 11:05:04 2015

@author: jeff

NOW REQUIRES mean_divg_values.txt
"""

### user setable variables ###

ref_dir = '/volumes/hd1/ref_genome_database/' # location of genome database created with genome_finder_build_core_genomes.py
ref_package = 'combined_16S' # name of reference package
pgdb_dir = '/home/jeff/ptools-local/pgdbs/user/' # location of pathway-tools pgdbs
strain_dir = '/volumes/hd1/ref_genome_database/ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' # location of strain genomes 

### end user setable variables ###

import sys
import os
import subprocess
from Bio import SeqIO, Phylo
import re
import numpy as np

name = sys.argv[1]
name_out = sys.argv[2]

## diagnostic only !
#name = 'palmer.good.filter.subsample.SRR584328_combined_16S.pplacer.filter.jplace.csv'
#name_out = 'SRR584328'

wd = os.getcwd() + '/'
query_dir = wd + '/' + name_out + '.query_genomes/'

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
    uid = re.split('uid', tip_name)
    uid = uid[1]
    uid = 'uid' + uid
    ref_edges[edge] = uid

## generate a dictionary mapping uid to strain for all reference strains

ref_uid_strain = {}

for ref in os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria'):
    ref_uid = re.split('uid', ref)
    ref_uid = ref_uid[1]
    ref_uid = 'uid' + ref_uid 
    ref_uid_strain[ref_uid] = ref

#### collect statistics on all reference genomes ####

## map divg to uid

divg_dict = {}
min_divg = 1

with open('mean_divg_values.txt', 'r') as divg_file:
    for line in divg_file:
        line = line.rstrip()
        line = line.split()
        
        line_uid = re.split('uid', line[0])
        line_uid = line_uid[1].rstrip('.combined')
        line_uid = 'uid' + line_uid
        
        divg_dict[line_uid] = line[1]

        if float(line[1]) < min_divg:
            min_divg = float(line[1])

## for core genomes, map node number to stats

stats_dict = {}

for d in internal_nodes:
    with open(ref_dir + ref_package + '.core_genomes' + '/' + d + '/' + d + '.stats', 'r') as stats:
        temp_size = []
        temp_divg = []
        
        for line in stats:
            line = line.rstrip()
            line = line.split()
            
            if line[0] != 'core':
                line_uid = line[0].rstrip('_combined.fasta')
                temp_size.append(line[1])
                
                ## there are some strains for which divg cannot be calculated, use lowest value for these                
                
                try:
                    temp_divg.append(divg_dict[line_uid])
                except KeyError:
                    temp_divg.append(min_divg)
                
            else:
                temp_divg = map(float, temp_divg)
                temp_divg = np.array(temp_divg)
                temp_mean_divg = np.mean(temp_divg)
                
                temp_size = map(float, temp_size)
                temp_size = np.array(temp_size)
                temp_mean = np.mean(temp_size)
                temp_sd = np.std(temp_size)
                core_size = line[1]
                
                stats_dict[d] = core_size, temp_mean, temp_sd, temp_mean_divg

## find edges with a reference sequence

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
    print >> tally_out, 'edge' + '\t' + 'nplacements' + '\t' + 'core.size' + '\t' + 'mean.size' + '\t' + 'sd.size ' + '\t' + '1-divg'
    
    for key in edge_tally:
        
        ## optimum score is number of placements * 1
        opt_sample_score.append(float(edge_tally[key]))
        
        try:
            edge_stats = stats_dict[str(key)]
            
            ## sample score = nplacements * (core genome size / mean cluster genome size) * 1 - meandivg
            sample_score.append(float(edge_tally[key]) * (float(edge_stats[0]) / float(edge_stats[1])) * (1 - float(edge_stats[3])))
            
        ## terminal nodes will not have an entry in stat_dict
        except KeyError:            
            edge_uid = ref_edges[key]
            edge_divg = divg_dict[edge_uid]
            edge_stats = 'NA', 'NA', 'NA', edge_divg
        
            ## in this case sample score is just nplacements * (1 - divg)
            sample_score.append(float(edge_tally[key]) * (1 - float(edge_stats[3])))
            edge_name = query_edges[key]
    
        print >> tally_out, str(key) + '\t' + str(edge_tally[key]) + '\t' + str(edge_stats[0]) + '\t' + str(edge_stats[1]) + '\t' + str(edge_stats[2]) + '\t' + str(1 - float(edge_stats[3]))

    sample_score = sum(sample_score) / sum(opt_sample_score)
    print >> tally_out, 'sample.score' + '\t' + str(sample_score)
    
pgdbs = set(os.listdir(pgdb_dir))

## setup core and strain genome directories for pgdb creation, if needed 

with open('generate_pgdbs.sh', 'w') as run_pgdb:
    print >> run_pgdb, '#!/bin/bash'
    
    for key in edge_tally.keys():
        if key + 'cyc' not in pgdbs:
            if key in internal_nodes:
                
                with open(ref_dir + ref_package + '.core_genomes/' + key + '/' + 'organism-params.dat', 'w') as organism_params:
                    print >> organism_params, 'ID' + '\t' + key
                    print >> organism_params, 'Storage' + '\t' + 'File'
                    print >> organism_params, 'Name' + '\t' + key
                    print >> organism_params, 'Rank' + '\t' + 'Strain'
                    print >> organism_params, 'Domain' + '\t' + 'TAX-2'
                    print >> organism_params, 'Create?' + '\t' + 't'
                    
                with open(ref_dir + ref_package + '.core_genomes/' + key + '/' + 'genetic-elements.dat', 'w') as genetic_elements:            
                    print >> genetic_elements, 'ID' + '\t' + key + '.1'
                    print >> genetic_elements, 'NAME' + '\t' + key + '.1'
                    print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                    print >> genetic_elements, 'CIRCULAR?' + '\t' + 'N'
                    print >> genetic_elements, 'ANNOT-FILE' + '\t' + key + '_core_genome.gbk'
                    print >> genetic_elements, '//'
                    
                print >> run_pgdb, 'pathway-tools -lisp -patho ' + ref_dir + ref_package + '.core_genomes/' + key + '/ -disable-metadata-saving &> pathos_' + key + '.log'   
                    
            else:
                uid = ref_edges[key]
                strain = ref_uid_strain[uid]
                
                with open(strain_dir + strain + '/' + 'organism-params.dat', 'w') as organism_params, open(strain_dir + strain + '/' + 'genetic-elements.dat', 'w') as genetic_elements:                

                    print >> organism_params, 'ID' + '\t' + key
                    print >> organism_params, 'Storage' + '\t' + 'File'
                    print >> organism_params, 'Name' + '\t' + key
                    print >> organism_params, 'Rank' + '\t' + 'Strain'
                    print >> organism_params, 'Domain' + '\t' + 'TAX-2'
                    print >> organism_params, 'Create?' + '\t' + 't'
                    
                    g = 0
                    for gbk in os.listdir(strain_dir + strain):
                        if gbk.endswith('gbk'):
                            g = g + 1
                            
                            print >> genetic_elements, 'ID' + '\t' + key + '.' + str(g)
                            print >> genetic_elements, 'NAME' + '\t' + key + '.' + str(g)
                            print >> genetic_elements, 'TYPE' + '\t' + ':CHRSM'
                            print >> genetic_elements, 'CIRCULAR?' + '\t' + 'N'
                            print >> genetic_elements, 'ANNOT-FILE' + '\t' + gbk
                            print >> genetic_elements, '//'
                    
                print >> run_pgdb, 'pathway-tools -lisp -patho ' + strain_dir + strain + '/ -disable-metadata-saving &> pathos_' + key + '.log'   
                
  