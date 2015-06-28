# -*- coding: utf-8 -*-
"""
Created on Tue Jan 06 09:50:07 2015

@author: jeff
"""

### user setable variables ###

ref_dir = '/volumes/hd1/ref_genome_database_v1/' # location of the database directory

### end user setable variables ###

from Bio import Phylo, SeqIO
import subprocess
import os
import re
from joblib import Parallel, delayed
import sys

executable = '/bin/bash'
query = sys.argv[1]
ref = 'combined_16S'

print 'getting edges from reference tree'

tree = Phylo.read(query, 'phyloxml')

nodes = {}

for clade in tree.get_nonterminals():
    
    members = []    
    try:
        edge = int(clade.confidence)
        for tip in clade.get_terminals():
            name = tip.name
            name = name.strip('@ref_')
            uid = re.split('uid', name)
            uid = uid[1]
            uid = 'uid' + uid
            print edge, name, uid
            members.append(uid)
        nodes[edge] = members
        
    except TypeError:
        pass
    
## make dictionary mapping uid to directory
## necessary because genome names in ref no longer match original directory names
    
print 'getting uids from reference genome names'
    
dirs = {}
    
for d in os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria'):
    uid = d[d.find('uid'):]
    dirs[uid] = d
    
## make a directory to hold the clade directories, then a directory for each clade
    
print 'building reference directory structure'
    
rm = subprocess.Popen('rm -r ' + ref_dir + ref + '.core_genomes', shell = True, executable = executable)
rm.communicate()
    
mkdir = subprocess.Popen('mkdir ' + ref_dir + ref + '.core_genomes', shell = True, executable = executable)
mkdir.communicate()
        
## for each clade, blast one genome against all others.  genes appearing in all are
## defined as core genes
        
def blast(ref_dir, ref, clade, nodes, dirs):
                
    d = ref_dir + ref + '.core_genomes/' + str(clade) + '/'
    name_dic = {}
    blast_dic = {}
    
    mkdir = subprocess.Popen('mkdir ' + d, shell = True, executable = executable)
    mkdir.communicate()
    
    genomes = nodes[clade]
    
    for genome in genomes:
        d_genome = dirs[genome]
        d_genome = ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' + d_genome
        cat = subprocess.Popen('cat ' + d_genome + '/*.gbk > ' + d + genome + '_combined.gbk', shell = True, executable = executable)
        cat.communicate()
        
    good = []
    clade_stats = open(d + str(clade) + '.stats', 'w')
    
    for f in os.listdir(d):
        if f.endswith('.gbk'):
            
            gbk_name = f
            fasta_name = gbk_name.rstrip('.gbk') + '.fasta'
            found = set()
            
            with open(d + fasta_name, 'w') as fasta_out:
                for record in SeqIO.parse(d + gbk_name, 'genbank'):
                    parent = record.seq
                    for feature in record.features:
                        
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        feature_id = (start, end, strand)
                        if feature_id not in found:
                            
                            if strand == -1:
                                sequence = parent[start:end].reverse_complement()
                            else:
                                sequence = parent[start:end]
                            
                            print start, end, sequence[0:10]
                            
                            seqname = record.id + '_' + str(start) + '_' + str(end) + '_' + str(strand)            
                            name_dic[record.id, start, end, strand] = seqname
                            
                            print >> fasta_out, '>' + seqname + '\n' + sequence
                            
                            found.add(feature_id)
            
            print >> clade_stats, fasta_name, len(found)            
            make_db = subprocess.Popen('cd ' + d + ';makeblastdb -in ' + fasta_name + ' -dbtype nucl', stderr = subprocess.PIPE, shell = True, executable = executable)
            error = make_db.communicate()[1]
            
            ## sometime database build fails inexplicably.  if this happens and total number of searches is not corrected
            ## there will be no core genes found            
            
            if error == '':
                good.append(fasta_name)
    
    use_file = True # only need to do one blast! 
    n = len(good)
           
    for f in good:
        if use_file == True:
            blast_query_file = f
            for f2 in good:
                if f2 != f:
                    blastn = subprocess.Popen('cd ' + d + ';blastn -task dc-megablast -num_alignments 1 -evalue 1 -query ' + f + ' -db ' + f2 + ' -out ' + f2 + '.txt -outfmt 6', shell = True, stderr = subprocess.PIPE, executable = executable)
                    error = blastn.communicate()[1]
                    
                    if error != '':
                        n = n - 1
                    
                    with open(d + f2 + '.txt', 'r') as blast_results:
                        for line in blast_results:
                            line = line.rstrip()
                            line = line.split('\t')
                            query = line[0]
                            ref = line[1]
                            try:
                                temp = blast_dic[query]
                                temp.append(ref)
                                blast_dic[query] = temp
                            except KeyError:
                                blast_dic[query] = [ref]
                        
            use_file = False
    
    ## determine which query genes have n - 1 hits, indicating they are present in all genomes
    
    core_gene = set()
    core_gene_name = set()
    
    with open(d + str(clade) + '_blast_results.txt', 'w') as blast_out:                        
        for query in blast_dic.keys():
            if len(blast_dic[query]) == n - 1:
                core_gene.add(query)
                print >> blast_out, query, '\t', blast_dic[query]
            
    ## identify the features that are represented by core genes
            
    for name in name_dic.keys():
        if name_dic[name] in core_gene:
            core_gene_name.add(name)
            
    ## find those features in the representative genbank file for this group
            
    core_records = []
    core_cds_n = 0
    
    with open(d + str(clade) + '_core_search.txt', 'w') as core_out:
        for record in SeqIO.parse(d + blast_query_file.rstrip('fasta') + 'gbk', 'genbank'): 
            i = 0
            temp_features = []
            
            for feature in record.features:
                
                i = i + 1
        
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                        
                if (record.id, start, end, strand) in core_gene_name:
                    temp_features.append(feature)
                    print >> core_out, i, feature.type, len(record.features), record.id, feature.location, 'core'
                    if feature.type == 'CDS':
                        core_cds_n = core_cds_n + 1
                
                elif feature.type == 'source':
                    temp_features.append(feature)
                    print >> core_out, i, feature.type, len(record.features), record.id, feature.location, 'core'
    
                else:
                    print >> core_out, i, feature.type, len(record.features), record.id, feature.location, 'not core'
            
            record.features = temp_features            
            core_records.append(record)
            
    print blast_query_file.rstrip('fasta') + 'gbk', core_cds_n, len(core_gene_name)
        
    SeqIO.write(core_records, open(d + str(clade) + '_core_genome.gbk', 'w'), 'genbank')
    
    ## delete nonessential files
    
    rm = subprocess.Popen('cd ' + d + ';rm *fasta;rm *nhr;rm *nin;rm *nsq;rm *combined.gbk', shell = True, executable = executable)
    rm.communicate()
    
    print >> clade_stats, 'core', core_cds_n    
    clade_stats.close()
    
if __name__ == '__main__':  
    Parallel(n_jobs = -1, verbose = 5)(delayed(blast)
    (ref_dir, ref, clade, nodes, dirs) for clade in nodes.keys())

## for testing    
#blast(ref_dir, ref, 511, nodes, dirs)
    



