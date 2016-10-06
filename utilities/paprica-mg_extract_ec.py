# -*- coding: utf-8 -*-
"""
Created on Thu Oct 06 10:57:13 2016

@author: jeff

This script extracts all of the gene translations in paprica-mg.ec.csv that
match a given EC number and makes a nice fasta file.  This script should be
executed from within the paprica/utilities folder.

Executing as:

python paprica-mg_extract_ec.py 4.1.1.21

Will produce the file paprica-mg_4.1.1.21.fasta, which contains all the
nonredundant gene products with that EC number in the paprica-mg database.
"""

import pandas as pd
import sys

target_ec = sys.argv[1]

mg_database = pd.read_csv('../ref_genome_database/paprica-mg.ec.csv', index_col = 0, header = 0)
mg_database_select = mg_database[mg_database['EC_number'] == target_ec]

with open('paprica-mg_' + target_ec + '.fasta', 'w') as fasta_out:
    for protein_id in mg_database_select.index:
        genome = mg_database.loc[protein_id, 'genome']
        domain = mg_database.loc[protein_id, 'domain']
        print 'found EC', target_ec, 'as', protein_id, 'in', genome
        print >> fasta_out, '>' + protein_id + '_' + genome + '_' + domain
        print >> fasta_out, mg_database_select.loc[protein_id, 'translation']
