# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 17:06:39 2015

@author: jeff
"""

### user setable variables ###

cpus = str(24) # number of cpus available
ref_dir = '/volumes/hd1/ref_genome_database/' # location of the database directory
tax_dir = '/volumes/deming/databases/' # location of the ncbi 16SMicrobial database

### end user setable variables ###

import os
import subprocess

executable = '/bin/bash'

## download genomes from Genbank

wget1 = subprocess.Popen('wget -r -A *.gbk -nc ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/ ' + ref_dir, executable = executable, shell = True)
wget2 = subprocess.Popen('wget -r -A *.fna -nc ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/ ' + ref_dir, executable = executable, shell = True)

with open('combined_16S.fasta', 'w') as fasta_out:
    for d in os.listdir(ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria'):
        cat = subprocess.Popen('cat ' + ref_dir + 'ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/' + d + '/*.fna > temp.fna', shell = True, executable = executable)
        cat.communicate()
        
        name = d            
 
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
        'python genome_finder_parse_blast_xml_v5.py temp.xml', shell = True, executable = executable)
        blast.communicate()
        
        with open('temp.fasta', 'r') as fasta_in:
            s = ''
            for line in fasta_in:
                line = line.rstrip('\n')
                if line.startswith('>') == False:
                    s = s + line
            print >> fasta_out, '>ref_'+name
            print >> fasta_out, s
            
        rm = subprocess.Popen('rm temp.xml;rm temp.fna', shell = True, executable = executable)
        rm.communicate()                   
