#!/usr/bin/env python
# -*- coding: utf-8 -*-

help_string = """
Created on Sun Jan 04 17:06:39 2015

@author: Jeff Bowman, jsbowman@ucsd.edu

paprica is licensed under a Creative Commons Attribution-NonCommercial
4.0 International License.  IF you use any portion of paprica in your
work please cite:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described
by Metabolic Structure: A General Framework and Application to a Seasonally
Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic
Peninsula." PloS one 10.8 (2015): e0135868.

If your analysis makes specific use of pplacer, Infernal, pathway-tools, or
any other software please make sure that you also cite the relevant publications.

REQUIRES:
    Files:
        kmer_top_1e5.txt
        bacterial_ssu.cm
        archaea_ssu.cm
        eukarya_ssu.cm

    Programs:
        raxmlHPC-PTHREADS-AVX2
        infernal
        seqmagick
        taxtastic

    Python modules:
        pandas
        Bio
        joblib
        Scipy
        Numpy
        
RUN AS:
    python paprica-make_ref.py [options]
    
OPTIONS:
    -cpus: The number of cpus for RAxML to use.
    -domain: Which domain are you analyzing?  Either bacteria, archaea, or eukarya.
    -download: Initiate a fresh download from Genbank?  Either T, F, or test.  Test
    allows you to use the small test set of genomes provided here: http://www.polarmicrobes.org/extras/ref_genome_database.tgz.
    -ref_dir: The name for the database you are building.  The default is "ref_genome_database".

This script must be located in the 'paprica' directory as it makes use of relative
paths.
    
"""

### User setable variables. ###

## If there are assemblies that you would like to exclude from analysis
## put them in the list 'bad' below.  Bedellvibrio, for example, causes
## errors for the placement of some reads.  For the eukarya, several of the
## 18S genes fall well outside the scope of the covariance model, presumably
## due to limitations in the training set.  This all need to be excluded or
## the tree is useless.

bad_bacteria=["spnn","suc","ppt","ppq","lsi","fpa","caj","car", \
"ayw","mtb","mtk","ppx","lce","lcs","hap","coo","lpm","bsl","lbl", \
"lbj","lbf","rai","tna","neu","tki","cso","hhm","abj","abh","srm", \
"cct","pay","bmw","lhl","lhk","lpz","psa","rim","bfg","aap","ste", \
"fto","ftm","ftn","gba","pgb","eca","cti","apw","kbt","apt","apu", \
"nmd","ctt","kct","bca","bue","pzu","balm","rbt","bbg"]


bad_archaea = []

bad_eukarya = ['MMETSP0017', \
'MMETSP0027', \
'MMETSP0103', \
'MMETSP0151', \
'MMETSP0200', \
'MMETSP0267', \
'MMETSP0403', \
'MMETSP0409', \
'MMETSP0414', \
'MMETSP0434', \
'MMETSP0448', \
'MMETSP0468', \
'MMETSP0562', \
'MMETSP0595', \
'MMETSP0689', \
'MMETSP0898', \
'MMETSP0918', \
'MMETSP0945', \
'MMETSP1015', \
'MMETSP1074', \
'MMETSP1317', \
'MMETSP1392', \
'MMETSP1446']

badGuys=bad_bacteria+bad_archaea+bad_eukarya

### End user setable variables. Mucking with anything below this point might ###
### break the script.  Of course it might also improve it :) ###

import os
import subprocess
from joblib import Parallel, delayed
import gzip
import re
import sys
import shutil
import sqlite3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio.Seq import Seq

import pandas as pd
import numpy as np
#from scipy import spatial

executable = '/bin/bash'
                
## Read in command line arguments.

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
        
## Define some variables based on these arguments.  If any arguments are
## missing replace with a default value.
        
try:        
    domain = command_args['domain']
except KeyError:
    domain = 'eukarya'
try:
    cpus = str(command_args['cpus'])
except KeyError:
    cpus = str(8)
try:
    download = command_args['download']
except KeyError:
    download = 'T'
try:
    ref_dir = command_args['ref_dir']
except KeyError:
    ref_dir = 'ref_genome_database'

        
## Make sure that ref_dir ends with /.
    
if ref_dir.endswith('/') == False:
    ref_dir = ref_dir + '/'

os.mkdir(ref_dir)
os.mkdir(ref_dir+"bacteria/")
os.mkdir(ref_dir+"archaea/")
os.mkdir(ref_dir+"eukarya/")


paprica_path = os.path.dirname(os.path.abspath("__file__")) + '/' # The location of the actual paprica scripts.    

ref_dir_domain = paprica_path + ref_dir + domain + '/'

## Select model based on domain specified.
    
if domain == 'bacteria':
    cm = paprica_path + 'models/bacteria_ssu.cm'
elif domain == 'archaea':
    cm = paprica_path + 'models/archaea_ssu.cm'
elif domain == 'eukarya':
    cm = paprica_path + 'models/eukarya_ssu.cm'
else:
    print 'Error, you must specify either -domain bacteria, -domain archaea, or -domain eukarya!'
    #quit()

#%% Define functions

## Define a stop function for diagnostic use only.

def stop_here():
    stop = []
    print 'Manually stopped!'
    print stop[1]


def getOrganisms():
    ## Define a dictionary all organisms. def get
    genomeDomain={}
    scientificName={}
    genomeList=os.popen("curl rest.kegg.jp/list/organism").read()

    #Get dictionary of all organisms and domain. Should skip bad organisms here
    for line in genomeList.split('\n'):
        if len(line)>0:
            info= line.split('\t')
        if info[1] not in badGuys:
            genomeDomain[info[1]]= info[3].split(';')[1]
            scientificName[info[1]]=info[2]
    return(genomeDomain,scientificName)


 
def findGenes():
    #Need to grab a 16S gene identifier. Just grab the first one. 
    sixteenSDict={}
    for org in genomeDomain:
        sixteenS=os.popen("curl rest.kegg.jp/link/ko/"+org).read()
        for line in sixteenS.split('\n'):
            if len(line)>0:
                values=line.split('\t')
                ko= values[1].strip('ko:')
                if ko=="K01977":
                    gene=values[0].split(':')	
                    sixteenSDict[org]=gene[1]

#Grab 16S sequences and split into respective domains fasta file

## added these files after fixing the sys.stdout redirect and assume it will work, but has not been tested. 
def buildFasta():
    bacteriaFasta = open(ref_dir+"bacteria/16S_bacteria_kegg.fasta", 'w+')
    archaeaFasta = open(ref_dir+"archaea/16S_archaea_kegg.fasta", 'w+')
    for org,gene in sixteenSDict.items():
        info=os.popen("curl rest.kegg.jp/get/"+org+":"+gene).read()
        seq=False
        seqCharacters=""
        for line in info.split('\n'):
            if "NTSEQ" in line:
                seq=True
                continue
            if seq==True:
                if "///" in line and genomeDomain.get(org)=="Bacteria":
                    bacteriaFasta.write(">"+org+"\n")
                    bacteriaFasta.write(seqCharacters+"\n")
                if "///" in line and genomeDomain.get(org)=="Archaea":
                    archaeaFasta.write(">"+org+"\n")
                    archaeaFasta.write(seqCharacters+"\n")
                seqCharacters=seqCharacters+line.strip(" ")

def buildEnzymeDB(enzyme,genomeDomain):
    DB=open(ref_dir+"paprica"+enzyme+"_DB.txt",'w+')
    for org in genomeDomain.keys():
        genes={}
        result=os.popen("curl http://rest.kegg.jp/link/"+enzyme.lower()+"/"+org).read()
        for row in result.split("\n"):
            if len(row) <1:
                break
            gene_dat=row.strip("\n").split("\t")[1].strip(enzyme.lower()+":")
            if gene_dat not in genes:
                genes[gene_dat]=1
            else:
                genes[gene_dat]=genes.get(gene_dat)+1
        for keys,values in genes.items():
            DB.write(org+"\t"+keys+"\t"+str(values)+"\n")
    DB.close()

def collectNamesNCBI():
    ## This is my way of resolving ncbi and KEGG. It is a terribly long step. It works, but its long. 
    ## Feel free to fix if you like. 
    os.popen('taxit new_database '+ref_dir+'ncbi_taxonomy.db -p '+ref_dir)
    conn = sqlite3.connect(ref_dir+'ncbi_taxonomy.db')
    c=conn.cursor()
    c.execute('SELECT tax_id,tax_name from names')
    names=c.fetchall()
    badChar=re.compile('[^a-zA-Z]')
    strain={}
    species={}
    genus={}
    for row in names:
        length= len(row[1].split(' '))
        g=row[1].split(' ')[0]
        sp=[]
        str=[]
        if length==2:
            sp=row[1].split(' ')[1]
            species[g+" "+sp]=row[0]
        if length>2:
            str=badChar.sub('',row[1]).lower()
            strain[str]=row[0]
        genus[g]=row[0]


def matchTaxaNames():
    seq_info={}	
    for org,name in scientificName.items():
        check=False
        taxa= name.split(' ')
        print org, name
        for i in taxa:
            if i in genus.keys():
                g=i
                seq_info[org]=genus.get(i)
                check=True
            if g+" "+i in species.keys():
                seq_info[org]=species.get(g+" "+i)
        name= badChar.sub('',name).lower()
        if name in strain.keys():
            seq_info[org]=strain.get(name)
        if not check and genomeDomain.get(org)=="Archaea":
            seq_info[org]='2157'
        if not check and genomeDomain.get(org)=="Bacteria":
            seq_info[org]='2'

def writeSeq_InfoAndTaxa(domain,genomeDomain,seq_info):
    seqInfo=open(ref_dir+domain+"/seq_info.csv",'w+')
    tax_ids=open(ref_dir+domain+"/tax_ids.txt",'w+')
    seqInfo.write("seqname,accesion,tax_id,species_name,is_type\n")
    for keys,values in seq_info.items():
        if genomeDomain.get(keys)==domain:
            seqInfo.write(keys+',,'+values+',,\n')
            tax_ids.write(values+'\n')
        else:
            continue


def main():
    (genomeomain,scientificName)=getOrganisms()
    findGenes()
    buildFasta()
    buildEnzymeDB("KO",genomeDomain)
    matchTaxaNames()
    for domain in ["Archaea","Bacteria"]:
    	writeSeq_InfoAndTaxa(domain,genomeDomain,seq_info)

if __name__ == "__main__":
    main()
