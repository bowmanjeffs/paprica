# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:34:24 2013

@author: Jeff

#### WARNING #### Script needs to be updated in accordance with changes
reflected in parse_blast_xml_gz_v3_for_network... namely use of re.sub

UPDATE: The blast xml format can return multiple alignments for a single
hit, this awk step in this script removes redundant alignments.  The
fasta file returned by this script should be unique.

This script takes as input a gzip xml file produced from a blast search
and returns a table of hits that is easier to read and tally.  If fasta_out
is set to True it also rights out a fasta file of the query seqs, as 
amino acid if blastx is used.

parse_blast_xml_gz.py [file_in].xml.gz

outputs: [file_in].fasta and [file_in].txt
"""
#####################################################################

fasta_out = True

#####################################################################

import re
import sys
import gzip
import subprocess

base = sys.argv[1]
base = base.rstrip('xml.gz')

if fasta_out == True:
    seq_out = open(base+'.fasta', 'w')

output = open(base+'.txt','w')
n = 0
print >> output, 'read'+'\t'+'hit_def'+'\t'+'hit_acc'+'\t'+'e'
read_def = set()

if sys.argv[1].endswith('gz'):
    with gzip.open(sys.argv[1],'rb') as xml:
            for line in xml:
                if re.search('<Iteration_query-def>', line) != None:
                    n = n + 1
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Iteration_query-def>', '', line)
                    line = re.sub('</Iteration_query-def>', '', line)
                    query_def = line
                if re.search('No hits found', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Iteration_message>', '', line)
                    line = re.sub('</Iteration_message>', '', line)
                    print n
                    print >> output, query_def+'\t'+line
                if re.search('<Hit_def>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hit_def>', '', line)
                    line = re.sub('</Hit_def>', '', line)
                    hit_def = line
                if re.search('<Hit_accession>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hit_accession>', '', line)
                    line = re.sub('</Hit_accession>', '', line)
                    hit_acc = line       
                if re.search('<Hsp_evalue>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hsp_evalue>', '', line)
                    line = re.sub('</Hsp_evalue>', '', line)
                    e_val = line
                if re.search('<Hsp_qseq>', line) != None:
                    if query_def not in read_def:
                        read_def.add(query_def)
                        print n
                        print >> output, query_def+'\t'+hit_def+'\t'+hit_acc+'\t'+e_val
                        line = line.strip()
                        line = line.rstrip()
                        line = re.sub('<Hsp_qseq>', '', line)
                        line = re.sub('</Hsp_qseq>', '', line)               
                        seq = line
                        if fasta_out == True:
                            print >> seq_out, '>'+query_def+'\n'+seq
                            
else:
    with open(sys.argv[1],'rb') as xml:
            for line in xml:
                if re.search('<Iteration_query-def>', line) != None:
                    n = n + 1
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Iteration_query-def>', '', line)
                    line = re.sub('</Iteration_query-def>', '', line)
                    query_def = line
                if re.search('No hits found', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Iteration_message>', '', line)
                    line = re.sub('</Iteration_message>', '', line)
                    print n
                    print >> output, query_def+'\t'+line
                if re.search('<Hit_def>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hit_def>', '', line)
                    line = re.sub('</Hit_def>', '', line)
                    hit_def = line
                if re.search('<Hit_accession>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hit_accession>', '', line)
                    line = re.sub('</Hit_accession>', '', line)
                    hit_acc = line       
                if re.search('<Hsp_evalue>', line) != None:
                    line = line.strip()
                    line = line.rstrip()
                    line = re.sub('<Hsp_evalue>', '', line)
                    line = re.sub('</Hsp_evalue>', '', line)
                    e_val = line
                if re.search('<Hsp_qseq>', line) != None:
                    if query_def not in read_def:
                        read_def.add(query_def)
                        print n
                        print >> output, query_def+'\t'+hit_def+'\t'+hit_acc+'\t'+e_val
                        line = line.strip()
                        line = line.rstrip()
                        line = re.sub('<Hsp_qseq>', '', line)
                        line = re.sub('</Hsp_qseq>', '', line)               
                        seq = line
                        if fasta_out == True:
                            print >> seq_out, '>'+query_def+'\n'+seq    

output.close()
if fasta_out == True:
    seq_out.close()
    
subprocess.call('awk \'!($0 in a) {a[$0];print}\' '+base+'.txt > '+base+'_unique.txt', shell = True)
subprocess.call('fgrep -v \'No hits found\' '+base+'_unique.txt > '+base+'_unique_hits.txt', shell = True)
                
                
