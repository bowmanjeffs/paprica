# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 21:39:12 2022

@author: jeff
"""

import pandas as pd
import os

ref_dir = 'ref_genome_database'
pgdb_dir = '/volumes/hd2/ptools-local/pgdbs/user/'
    
pathway_definitions = pd.read_csv(ref_dir + '/pathways.col', comment = '#', sep = '\t', index_col = 0, low_memory = False)

for i,d in enumerate(os.listdir(pgdb_dir)):
#for d in ['gca_013364155.1cyc']:
    n_paths = 0 # Number of pathways predicted.
    try:
        
        ## As of ptools v24 pathway-report.txt reformatted with date.  Need
        ## to find name of this file.
        
        report_file = None
        
        for f in os.listdir(pgdb_dir + d + '/1.0/reports'):
            if f.startswith('pathways-report_'):
                report_file = f
            elif f.startswith('pwy-inference-report'):
                inference_file = f
                
        #!!! For reasons that aren't clear, pathway-tools is failing at the
        ## end of the prediction, so the pathways-report file isn't being
        ## created.  However, the pwy-inference-report is created, and
        ## this can be parsed to create the pathways-report. This pathway
        ## report has the old-style name, so it is deleted each time, in the
        ## hope that eventually pathologic works correctly.
                
        if report_file == None:
            report_file = 'pathways-report.txt'
                
        ## Now the report file can be parsed.
            
        temp = []
        
        with open(pgdb_dir + d + '/1.0/reports/' + report_file, 'r') as report:            
            for line in report:
                if line.startswith('#') == False & line.startswith('Pathway Name') == False:
                    line = line.rstrip()
                    if line != '':
                        line = line.split('|')                                
                        if len(line) > 0:
                            path = line[0]
                            path_id = line[1]
                            temp.append(path_id)
                            pathway_definitions.loc[path_id, 'NAME'] = path
                            
        print(d, i)
                            
    except (NameError, IOError):
        print(d, 'has no pathway report')
        
pathway_definitions.to_csv(ref_dir + '/pathways.col', sep = '\t')