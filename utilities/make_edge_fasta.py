# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 16:58:12 2015

@author: jeff

CALL AS:
    python extract_pplacer_seqs.py [something.csv] [something.fasta] [start-inclusive] [stop-not inclusive]
    
This script makes a fasta file of the reads associated with a given range of edges.  If you
would like the reads associated within only one edge (such as 1234), you would use edge edge+1
as the arguments (1234 1235).

The input csv file should be the csv file produced by guppy (called in that paprica_place_it script).
It probably has a name that follows the form yoursample.combined_16S.tax.clean.align.csv.

The input fasta should be any fasta file that contains the cleaned reads names.
yoursample.combined_16S.tax.clean.fasta is a good bet.
    
"""

from Bio import SeqIO
import sys

csv = sys.argv[1]
fasta = sys.argv[2]
start = int(sys.argv[3])
stop = int(sys.argv[4])

sample = fasta.rstrip('.fasta')

get = set(map(str, range(start, stop))) # not inclusive of last number
to_get = set()

with open(csv, 'r') as csv_in:
    for line in csv_in:
        line = line.split(',')
        if line[3] in get:
            name = line[1]
            to_get.add(name)
            
with open(sample + '_' + str(start) + '_' + str(stop) + '.fasta', 'w') as fasta_out:
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id in to_get:
            print 'found', record.id
            SeqIO.write(record, fasta_out, 'fasta')
        
        