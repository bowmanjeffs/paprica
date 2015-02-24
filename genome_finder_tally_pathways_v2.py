# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 18:26:20 2015

@author: jeff
"""
pgdb_dir = '/home/jeff/ptools-local/pgdbs/user/' # location of pathway-tools pgdbs

import sys

name = sys.argv[1]
#name = 'SRR584327'
all_pathways = {}

with open(name + '.edge_tally.txt', 'r') as edge_tally_file, open(name + '.pathway_summary.txt', 'w') as path_out, open(name + '.pathway_detail.txt', 'w') as detail_out:
    for line in edge_tally_file:
        line = line.rstrip()
        line = line.split()
        if line[0] == 'edge':
            pass
        elif line[0] == 'sample.score':
            sscore = line[1]
        else:
            edge = line[0]
            nedge = line[1]
            
            try:
                with open(pgdb_dir + edge + 'cyc/1.0/reports/pathways-report.txt', 'r') as report:
                    print >> detail_out, edge + '\t' + nedge + '\t',
                    
                    for line in report:
                        if line.startswith('#') == False:
                            if line.startswith('Pathway Name') == False:
                                line = line.rstrip()
                                if line != '':
                                    line = line.split('|')                                
                                    if len(line) > 0:
                                        path = line[0]
        
                                        try:
                                            temp = all_pathways[path]
                                            temp = temp + int(nedge)
                                            all_pathways[path] = temp
                                        except KeyError:
                                            all_pathways[path] = int(nedge)
                                        
                                        print >> detail_out, path + '\t',
                    print >> detail_out, '\n',

            except IOError:
                print edge, 'has no pathway report'
                print >> detail_out, edge + '\t' + 'has no pathway report'

    for path in sorted(all_pathways.keys()):
        print >> path_out, path + '\t' + str(all_pathways[path]) + '\t' + sscore
        
                            
        
        