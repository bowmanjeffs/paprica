# -*- coding: utf-8 -*-
help_string = """
Created on Thu Jun 25 16:58:12 2015

@author: jeff

CALL AS:
    python make_edge_fasta.py -csv [something.csv] -fasta [something.fasta] -start [start-inclusive] -stop [stop-not inclusive]
    
This script makes a fasta file of the reads associated with a given range of edges.  If you
would like the reads associated within only one edge (such as 1234), you would use edge edge+1
as the arguments (1234 1235).

The input csv file should be the csv file produced by guppy (called in that paprica_place_it script).
It probably has a name that follows the form yoursample.combined_16S.tax.clean.align.csv.

The input fasta should be any fasta file that contains the cleaned reads names.
yoursample.combined_16S.tax.clean.fasta is a good bet.

Identified reads will be named edgenumber_original read name
    
"""

from Bio import SeqIO
import sys

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in command_args.keys():
    print help_string
    quit()

csv = command_args['csv']
fasta = command_args['fasta']
start = int(command_args['start'])
stop = int(command_args['stop'])

fasta_name = fasta.rstrip('.fasta')

get = set(map(str, range(start, stop))) # not inclusive of last number
to_get = set()
names = {}

with open(csv, 'r') as csv_in:
    for line in csv_in:
        line = line.split(',')
        if line[3] in get:
            name = line[1]
            to_get.add(name)
            names[name] = str(line[3]) + '_' + name
            
with open(fasta_name + '_' + str(start) + '_' + str(stop) + '.fasta', 'w') as fasta_out:
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id in to_get:
            print 'found', record.id
            record.id = names[str(record.id)]
            SeqIO.write(record, fasta_out, 'fasta')
        
        