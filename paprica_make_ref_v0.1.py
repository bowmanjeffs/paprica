# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 17:06:39 2015

@author: jeff
"""

### user setable variables ###

cpus = str(24) # number of cpus available
ref_dir = '/volumes/hd1/ref_genome_database_v1/' # location of the database directory
tax_dir = '/volumes/deming/databases/' # location of the ncbi 16SMicrobial database
download = False # set to true to initiate fresh download of genomes

## if there are genomes that you would like to exclude from analysis
## put them in a list here.  Salmonella uid61197, for example, is not
## actually in finished state.  It throws of some stats if included.
bad = ['Salmonella_enterica_serovar_Weltevreden_2007_60_3289_1_uid61197']

### end user setable variables ###

import os
import subprocess

executable = '/bin/bash'

## download genomes from Genbank

if download == True:

    wget1 = subprocess.Popen('cd ' + ref_dir + ';wget -r -A *.gbk -nc ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/', executable = executable, shell = True)
    wget1.communicate()
    
    wget2 = subprocess.Popen('cd ' + ref_dir + ';wget -r -A *.fna -nc ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/', executable = executable, shell = True)
    wget2.communicate()
    
    ## eliminate bad genomes
    if len(bad) > 0:    
        for b in bad:        
            rm = subprocess.Popen('rm -r ' + ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' + b, executable = executable, shell = True)
            subprocess.communicate()
    
## read phi values
    
phi_dict = {}
    
with open('mean_phi_values.txt', 'r') as phi_file:
    for line in phi_file:
        line = line.rstrip()
        line = line.split()
        phi_name = line[0].rstrip('.combined')
        phi_uid = phi_name[phi_name.find('uid'):]
        phi_dict[phi_uid] = line[1]
        
## generate genome data and combined 16S file

with open(ref_dir + 'combined_16S.fasta', 'w') as fasta_out, open(ref_dir + 'genome_data.txt', 'w') as tally_out:
    for d in os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria'):
        cat = subprocess.Popen('cat ' + ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' + d + '/*.fna > temp.fna', shell = True, executable = executable)
        cat.communicate()
        
        name = d
        uid = name[name.find('uid'):]
        phi = 'not_found'
        
        try:
            phi = phi_dict[uid]
        except KeyError:
            print name, 'phi value not found'
            
        ## count genetic elements
            
        n_elements = 0
        elements = os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' + d)
        for element in elements:
            if element.endswith('gbk'):
                n_elements = n_elements + 1
                
        ## get genome size
                
        gs = ''                
        with open('temp.fna', 'r') as fna_in:
            for line in fna_in:
                if line.startswith('>') == False:
                    line = line.rstrip()
                    gs = gs + line
 
        ## find 16S rRNA gene
 
        print 'finding 16S genes'
            
        print name
        blast = subprocess.Popen('blastn ' \
        '-task megablast ' \
        '-num_threads ' + cpus + ' '\
        '-max_target_seqs 1 ' \
        '-evalue 1e-50 ' \
        '-db ' + tax_dir + '16SMicrobial ' \
        '-outfmt 5 ' \
        '-out temp.xml ' \
        '-query temp.fna;' \
        'python paprica_parse_blast_xml_v0.1.py temp.xml', shell = True, executable = executable)
        blast.communicate()
        
        ## write 16S rRNA gene to combined 16S rRNA gene file
        
        with open('temp.fasta', 'r') as fasta_in:
            s = ''
            for line in fasta_in:
                line = line.rstrip('\n')
                if line.startswith('>') == False:
                    s = s + line
            print >> fasta_out, '>ref_'+name
            print >> fasta_out, s
            
        ## write out info on genome
        
        t = 0    
        with open('temp.txt', 'r') as tally_in:
            for line in tally_in:
                t = t + 1
        print 'total 16S =', t - 1
        
        print >> tally_out, uid + '\t' + name + '\t' + str(phi) + '\t' + str(t - 1) + '\t' + str(len(gs)) + '\t' + str(n_elements)
            
        rm = subprocess.Popen('rm temp.*', shell = True, executable = executable)
        rm.communicate()                   
