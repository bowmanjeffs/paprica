# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:34:24 2013

@author: Jeff

"""
#####################################################################

keep = 1 ## number of queries with hit to write out

#####################################################################

import re
import sys
import subprocess

base = sys.argv[1]
base = base.rstrip('.xml')

n = 0
                            
with open(sys.argv[1],'rb') as xml, open(base+'.fasta', 'w') as seq_out, open(base+'.txt','w') as output:
    print >> output, 'read'+'\t'+'hit_def'+'\t'+'hit_acc'+'\t'+'e'
    for line in xml:
        if re.search('<Iteration_query-def>', line) != None:
            line = line.strip()
            line = line.rstrip()
            line = re.sub('<Iteration_query-def>', '', line)
            line = re.sub('</Iteration_query-def>', '', line)
            query_def = line
#        if re.search('No hits found', line) != None:
#            line = line.strip()
#            line = line.rstrip()
#            line = re.sub('<Iteration_message>', '', line)
#            line = re.sub('</Iteration_message>', '', line)
#            print >> output, query_def+'\t'+line
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
            n = n + 1
            print n
            print >> output, query_def+'\t'+hit_def+'\t'+hit_acc+'\t'+e_val
            line = line.strip()
            line = line.rstrip()
            line = re.sub('<Hsp_qseq>', '', line)
            line = re.sub('</Hsp_qseq>', '', line)               
            seq = line
            if n <= keep:
                print >> seq_out, '>'+query_def+'\n'+seq
				
#subprocess.call('fgrep -v \'No hits found\' '+base+'.txt > '+base+'_hits.txt', shell = True)
                
                
