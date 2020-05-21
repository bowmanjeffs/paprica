#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:17:10 2020

@author: jeff

This script iterates across all pathways.csv or ec.csv files for a given domain
in a directory.  For a target enzyme or pathway, it returns the number
contributed by each edge.

Run as:
    ./map_edge_pathways.py -domain [bacteria | archaea | eukarya] -[pathway | enzyme] [pathway or enzyme commission number]
Example:
    ./map_edge_pathways.py -domain bacteria -pathway "2,4,6-trichlorophenol degradation"
    
Note that many pathway names contain funny characters that may not behave well
on the command line.  We have not made any effort to deal with these.

"""

import pandas as pd
import os
import sys
import re

## Read in command line arguments.

command_args = {}

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        try:
            command_args[arg] = sys.argv[i + 1]
        except IndexError:
            command_args[arg] = ''
        
if 'h' in list(command_args.keys()):
    pass
    #print help_string ## no help sting
    quit()
    
if 'domain' in command_args.keys():
    domain = command_args['domain']

else:
    domain = 'bacteria'
    
if 'pathway' in command_args.keys():
    target = command_args['pathway']
    search_suffix = domain + '.pathways.csv'
    
elif 'enzyme' in command_args.keys():
    target = command_args['enzyme']
    search_suffix = domain + '.ec.csv'
    
else:
    target = '2,4,6-trichlorophenol degradation'
    search_suffix = domain + '.pathways.csv'
    
data = pd.DataFrame()

for f in os.listdir('.'):
    if f.endswith(search_suffix):
        print(f)
        temp = pd.read_csv(f, index_col = 0)
        temp_select = temp.loc[target]
        name = re.sub(search_suffix, '', f)
        temp_select.rename(name, inplace = True)
        data = pd.concat([data, temp_select], axis = 1)

data.fillna(0, inplace = True)        
pd.DataFrame.to_csv(data, target + '.' + domain + '.edge_mapping.csv')