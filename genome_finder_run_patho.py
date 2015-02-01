# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 11:51:28 2015

@author: jeff
"""





with open('SRP016030_accession.txt', 'r') as accession_file, open('SRP016030_pathos.sh', 'w') as output:
    print >> output, '#!/bin/bash'
    for line in accession_file:
        line = line.rstrip()
        print >> output, 'screen -S ' + line + ' ~/pathway-tools/pathway-tools -lisp -patho /volumes/hd1/palmer_hypothetical_metagenomes/' + line + '.query_genomes/;screen -d ' + line

