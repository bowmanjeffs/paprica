#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_string = """
Created on Thu Jun 25 16:58:12 2015

@author: jeff

CALL AS:
    make_edge_fasta.py -csv [unique.csv] -start [start-inclusive] -stop [stop-not inclusive]
    
This script makes a fasta file of the reads associated with a given range of edges.  If you
would like the reads associated within only one edge (such as 1234), you would use edge edge+1
as the arguments (1234 1235).

The input csv file should be the unique csv file produced by paprica (called in that paprica_place_it script).
It probably has a name that follows the form [yoursample].[domain].unique_seqs.csv.
    
"""

import sys
import pandas as pd

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in list(command_args.keys()):
    print(help_string)
    quit()
    
if len(command_args) == 0:
    command_args['csv'] = 'B234-P_B31_18S.18S.exp.eukarya.combined_18S.eukarya.tax.clean.unique.align.csv'
    command_args['start'] = 275
    command_args['stop'] = 277

csv = command_args['csv']
start = int(command_args['start'])
stop = int(command_args['stop'])

csv_in = pd.read_csv(csv, header = 0, index_col = 0)
select = set(csv_in.index[csv_in['edge_num'].between(start, stop)])
select = csv_in.loc[csv_in['edge_num'].between(start, stop)]
fasta = csv.split('.')
fasta = fasta[:-3]
fasta = '.'.join(fasta)
            
with open(fasta + '_' + str(start) + '_' + str(stop) + '.clean.unique.fasta', 'w') as fasta_out:
    for read in select.index:
        print('>' + read, file = fasta_out)
        print(select.loc[read, 'seq'], file = fasta_out)
        
        