# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 15:45:43 2015

@author: jeff

This script aggregates information from multiple '.edge_data.csv' files
produces by running paprica on multiple samples.  It produces a matrix of edges
by sample, and a matrix of mean edge parameters, by sample.

Run as: python combine_edge_results.py [prefix for output]

It will automatically loop through all .edge_data.csv files in the directory.

"""

import os
import pandas as pd
import re
import math
import sys

prefix = sys.argv[1]

edge_tally = pd.DataFrame()
edge_data = pd.DataFrame()

def fill_edge_data(param, name, df_in):
    temp = []
    for index in df_in.index:
        
        ## Exception is necessary for old version of build_core_genomes which
        ## did not provide a n16S value for draft genomes.
        
        try:
            n = range(int(math.ceil(df_in.loc[index, 'nedge_corrected'])))
        except ValueError:
            n = range(int(df_in.loc[index, 'nedge']))
            
        for i in n:
            temp.append(df_in.loc[index, param])
            
    temp = pd.Series(temp)
            
    mean = temp.mean()
    sd = temp.std()
    
    print name, param, mean, sd
    return mean, sd
        
for f in os.listdir('.'):
    if f.endswith('edge_data.csv'):
        
        temp_edge = pd.DataFrame.from_csv(f, index_col = 0)
        name = re.sub('.edge_data.csv', '', f)
        
        for param in ['n16S', 'nge', 'ncds', 'genome_size', 'GC', 'phi', 'confidence']:
            edge_data.loc[name, param + '.mean'], edge_data.loc[name, param + '.sd'] = fill_edge_data(param, name, temp_edge)
        
        for index in temp_edge.index:
            edge_tally.loc[name, index] = temp_edge.loc[index, 'nedge_corrected']
            
pd.DataFrame.to_csv(edge_tally, prefix + '.edge_tally.csv')
pd.DataFrame.to_csv(edge_data, prefix + '.edge_data.csv')        