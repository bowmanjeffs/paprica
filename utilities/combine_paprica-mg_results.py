# -*- coding: utf-8 -*-
"""
Created on Thu Jun 09 15:56:03 2016

@author: jeff

This script aggregates information from multiple '.mg.sum_ec.csv' files
produces by running paprica on multiple samples.  It produces a matrix of edges
by sample, and a matrix of mean edge parameters, by sample.

Run as: python combine_paprica-mg_results.py [prefix for output]

It will automatically loop through all files in the directory with the .mg.sum_ec.csv suffix.

"""

import pandas as pd
import re
import os
import sys

try:
    prefix = sys.argv[1]
except IndexError:
    prefix = 'test'

ec_tally = pd.DataFrame()
columns = []

for f in os.listdir('.'):
    if f.endswith('.mg.sum_ec.csv'):
        
        name = re.sub('.mg.sum_ec.csv', '', f)
        temp_ec = pd.read_csv(f, index_col = 0)               
        ec_tally = ec_tally.merge(temp_ec, how = 'outer', left_index = True, right_index = True)     
        columns.append(name)
        print name

ec_tally.columns = columns            
pd.DataFrame.to_csv(ec_tally.transpose(), prefix + '.mg.ec_tally.csv')
