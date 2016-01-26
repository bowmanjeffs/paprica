# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 21:47:16 2015

@author: jeff

This script parses a phyloxml format tree created by guppy fat for the names of
edges.  The names returned will be the names of the sequences used to build the
reference tree.  You can either specify a specific edge to lookup or create a
csv file of all terminal edges and the corresponding taxa.

Creating the csv file is deprecated for paprica, as this information is
included in the *edge_data.csv output.  However this feature has been retained
because it is useful in pplacer analysis outside of paprica.

Run as:

python paprica_ref_tree_quick_parse.py [csv|integer] [phyloxml format tree]

Using csv as argument creates a csv format file of all the terminal edges and
the strain names.  Using an integer as an argument prints out either the strain
name if a terminal edge or, if an internal edge, the strain name of all terminal
edges originating from that edge.

"""

import sys
from Bio import Phylo
    
create_csv = False 

if sys.argv[1] == 'csv':
    create_csv = True # generate a csv file of all terminal edges and taxa   
else:
    x = int(sys.argv[1]) # the edge you would like to identify

tree = sys.argv[2] # the name of a phyloxml tree generated using Guppy fat
tree = Phylo.read(tree, 'phyloxml')

if create_csv == True:
    with open('edge_taxa.csv', 'w') as csv_out:
        for clade in tree.get_terminals():
            try:
                print int(clade.confidence), clade.name
                print >> csv_out, str(int(clade.confidence)) + ',' + clade.name
            except TypeError:
                continue

else:    
    for clade in tree.get_terminals():
        try:
            if int(clade.confidence) == x:
                found = clade
        except TypeError:
            pass
    
    for clade in tree.get_nonterminals():
        try:
            #print int(clade.confidence)
            if int(clade.confidence) == x:
                found = clade.get_terminals()
        except TypeError:
            pass
    
    print found




