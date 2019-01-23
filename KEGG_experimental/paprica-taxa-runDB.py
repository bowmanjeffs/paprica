#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re
import argparse
import os
import numpy as np
import datetime
import subprocess
parser=argparse.ArgumentParser()
parser.add_argument("-i","--fasta")
parser.add_argument("-c","--counts")
parser.add_argument("-o","--outputDir",default="outFolder/")
parser.add_argument("-enz","--enzyme",default="ko")
parser.add_argument("-cut","--cutoff",default=0.5)
parser.add_argument("-m","--meta",default="all")

parser.add_argument('-f', action='store_true')

args=parser.parse_args()


def stop_here():
    sys.exit(1)

def writeSbatch(domain):
    sbatchFile=open(args.fasta.strip('fasta')+domain+"SBATCH.sh", 'w+')
    sbatchFile.write("#!/bin/bash\n"\
                    "## job options/settings\n"\
                    "#SBATCH --job-name="+args.fasta.strip('fasta')+domain+"\n"\
                    "#SBATCH --partition=batch\n"
                    "#SBATCH --nodes=1\n"
                    "#SBATCH --mem=30G\n"
                    "paprica-run.sh "+args.fasta.strip('fasta')+domain+".fasta "+domain)

nextline_1='\n'
nextline_2='\r'
fastaString=args.fasta.rsplit(".").pop()


def organizeCountTable(countTable):
    print "Storing OTU table"
    otuTaxa={}
    otuList=[]
    counts=[]
    with open(countTable, 'rb') as seqtable:
        sep='\t'
        line=seqtable.readline()
        if "#OTU" not in line:
            print 'Removing header line'
            line=seqtable.readline()
            if "#OTU" not in line:
                print '"#OTU" is detected in table. Double check header for compatability with biom files'
                stop_here()
        samples=line.strip(nextline_1).strip(nextline_2).split(sep)[1:]
        if len(samples)<2:
            sep=','
            samples=line.strip(nextline_1).strip(nextline_2).split(sep)[1:]

        for row in seqtable:
            info=row.strip(nextline_1).strip(nextline_2).split(sep)
            taxa=info[len(info)-1].split(';')[0]
            try:
                float(taxa)
                print "Encountered numeric taxonomy. This usually occurs if there is no taxonomy column."
                print "In the case there is taxonomy. This is what we are trying to classify: ",taxa
                continue
            except ValueError:
                otuTaxa[info[0]]=taxa.strip('k__')
                otuList.append(info[0])
                intCounts=[float(i) for i in info[1:(len(info)-1)]]
                counts.append(intCounts)

    domainList=[]
    for keys,values in otuTaxa.items():
        domainList.append(values)
    domainList=list(set(domainList))
    domainList=[x.lower() for x in domainList]

    if "bacteria" in domainList and "archaea" in domainList:
        print "Two domains detected. Paprica will run twice."
        print "Once for Bacteria and once for Archaea."
    if len([x for x in domainList if x !='bacteria' and x!='archaea'])!=0:
        oddballs=[x for x in domainList if x !='bacteria' and x!='archaea']
        print ','.join(oddballs)+' detected in taxonomy. Will try '+','.join(oddballs)+' running under bacteria covariance model and prediction'
        domainList=list(set(["bacteria" if x !="archaea" else x for x in domainList]))

    return counts,otuList,otuTaxa,domainList,samples

def splitFasta(otuTaxa,fasta,out):
    print 'Writing domain specific fasta files'
    name=os.path.basename(fasta).rsplit(fastaString,1)[0]
    bacteriaFasta=open(out+'/'+name+"bacteria.fasta", 'w+')
    archaeaFasta=open(out+'/'+name+"archaea.fasta", 'w+')

    with open(fasta,'rb') as fasta:
        for line in fasta:
            if ">" in line:
                otu=line.strip('>').strip(nextline_1).strip(nextline_2).split('_')[0]
                if otuTaxa.get(otu)=="Bacteria" or otuTaxa.get(otu)!="Archaea":
                    bacteriaFasta.write(">"+otu+"\n")
                if otuTaxa.get(otu)=="Archaea":
                    archaeaFasta.write(">"+otu+"\n")
            if ">" not in line:
                if otuTaxa.get(otu)=="Bacteria" or otuTaxa.get(otu)!="Archaea":
                    bacteriaFasta.write(line)
                if otuTaxa.get(otu)=="Archaea":
                    archaeaFasta.write(line)
        bacteriaFasta.close()
        archaeaFasta.close()

#writeSbatch("bacteria")#sbatching doesn't work because Descartes
#writeSbatch("archaea")#sbatching doesn't work because Descartes

def runPlaceIt(fasta,domain,out,log):
    place_it='paprica-place-it.py -query '+fasta+' -domain '+domain+' -o '+out
    print 'calling place-it as:'
    print place_it
    #os.system(place_it)
    placeIt=subprocess.Popen(place_it,stdout=subprocess.PIPE,shell=True)
    out,err=placeIt.communicate()
    print out,'' if err==None else err
    log.write(out)
	

def runTallyPathways(fasta,domain,out,enzyme,cutoff,log):
    print 'calling tally pathways as:'
    print 'paprica-tally-pathways-dev.py -ref_dir ref_genome_database -i '+fasta+'.taxit.clean.align.csv -domain '+domain+' -cutoff '+str(cutoff)+' -enzyme '+enzyme+' -o '+out
    tallyPathways=subprocess.Popen('paprica-tally-pathways-dev.py -ref_dir ref_genome_database -i '+fasta+'.taxit.clean.align.csv -domain '+domain+' -cutoff '+str(cutoff)+' -enzyme '+enzyme+' -o '+out,stdout=subprocess.PIPE,shell=True)
    out,err=tallyPathways.communicate()
    print out,'' if err==None else err
    log.write(out)

#{fasta filename}.{domain}.unique.seqs
#this is the map of sequence identifiers (or OTU #) to edge


def map_edges_seqs(domain,fasta):
    seq_edge=".unique.seqs.csv"
    seq_path=os.path.dirname(fasta)+'/placeIt/'
    seq_edge_map=seq_path+os.path.basename(fasta).rsplit(fastaString,1)[0]+domain+seq_edge
    seq_edge_dict={}
    edgeSeqDict={}
    with open(seq_edge_map,'rb') as csvfile:
        csvfile.readline()
        for row in csvfile:
            temp_array=row.strip('\n').split(',')
            seq_edge_dict[temp_array[3].split('_')[0]]=temp_array[2]
            if temp_array[2] not in edgeSeqDict:    
                edgeSeqDict[temp_array[2]]=[temp_array[3].split('_')[0]]
            else:
                edgeSeqDict[temp_array[2]].append(temp_array[3].split('_')[0])
    return edgeSeqDict


def map_edges_enzyme(domain,fasta,enz):
    #{fasta filename}.{domain}.ec e.g. seq.bacteria.ec uses seq.fasta for bacteria
    #this is the map of mapped edges to enzyme commision numbers
    edge_enz="."+enz+".csv"
    edge_path=os.path.dirname(fasta)+'/tallyPathways/'
    edge_enz_map=edge_path+os.path.basename(fasta).rsplit(fastaString,1)[0]+domain+edge_enz

#this is the map of total counts of EC's identified within the loaded sequence data as mapped per edge. As multiple sequences map to the same edge, you will eventually have to divide by the number of sequences that mapped to the edge to get EC per edge.
    with open(edge_enz_map,'rb') as csvfile:
        new_list=[row.strip('\n').split(',') for row in csvfile]

    enz_edge=zip(*new_list)
    #print ec_edge[1]
    enz_edge_dict={}
    for j in range(0,len(enz_edge)):
        enz_edge_dict[enz_edge[j][0]]=j
    return enz_edge,enz_edge_dict
#this creates a map between the sequence and the edge within paprica


def consolidateOtuGenes(fasta,domain,counts,otuList,enz):
    print 'Converting OTU counts into gene counts for '+domain
    edge_seq_dict=map_edges_seqs(domain,fasta)
    (enz_edge,edge_enz_dict)=map_edges_enzyme(domain,fasta,enz)

    koTable=[]
    for edge,otus in edge_seq_dict.items():
        n_edges=len(otus)
        edgeOTU=[]
        [edgeOTU.append(counts[k]) for k, org in enumerate(otuList) if org in otus]
        print edge, otus
        edge_counts=[float(k)/n_edges for k in enz_edge[edge_enz_dict.get(edge)][1:]]
        sumEdgeOTU=[(sum(x)) for x in zip(*edgeOTU)]
        ko_table=[]
        [ko_table.append([ko*float(count) for count in sumEdgeOTU]) for ko in edge_counts]
        koTable.append(ko_table)
    return(enz_edge,koTable)

def combineMetagenomes(enzEdgeAll, koTableAll):
    enzPosition=set(enzEdgeAll[0]).union(set(enzEdgeAll[1]))
    koTable=[]
    enz_edge=[['']]
    for enz in (enzPosition):
        if enz in enzEdgeAll[0] and enz not in enzEdgeAll[1]:
            koTable.append(koTableAll[0][enzEdgeAll[0].index(enz)-1])
            enz_edge[0].append(enz)

        if enz in enzEdgeAll[1] and enz not in enzEdgeAll[0]:
            koTable.append(koTableAll[1][enzEdgeAll[1].index(enz)-1])
            enz_edge[0].append(enz)

        if enz in enzEdgeAll[0] and enz in enzEdgeAll[1]:
            koTable.append([x+y for x,y in zip(koTableAll[0][enzEdgeAll[0].index(enz)-1],koTableAll[1][enzEdgeAll[1].index(enz)-1])])
            enz_edge[0].append(enz)

    return enz_edge,koTable


def writeItOut(out,koTable,enz_edge,samples):	
    metagenomeFile=open(out,'w+')
    metagenomeFile.write('#Metagenome Table\t'+'\t'.join(samples[:(len(samples)-1)])+'\n')
    [metagenomeFile.write(enz_edge[0][ko+1]+'\t'+'\t'.join([str(count) for count in koTable[ko]])+'\n') for ko in range(0,len(koTable))]

	

def main():
    outDir=args.outputDir
    if outDir.endswith('/') == False:
        outDir = outDir + '/'

    print 'Data being stored in '+outDir
    
    try :
        os.makedirs(outDir)
        print "making output directory"
    except OSError:
        if args.f:
            pass
        else:
            print 'Directory already exists. Use flag "-f" to force overwrite files.'
            stop_here()
    logFile=open(outDir+'paprica'+datetime.datetime.now().strftime("%d%m%y%H%M")+'.txt','w+')

    (counts,otuList,otuTaxa,domainList,samples) = organizeCountTable(countTable=args.counts)					
    koTableAll=[]
    enzEdgeAll=[]
 
    splitFasta(otuTaxa,fasta=args.fasta,out=outDir)
    for domain in domainList:
        if domain=="unassigned":
            continue
        #print os.path.basename(args.fasta)
        #print os.path.basename(args.fasta).rsplit(fastaString,1)
        runPlaceIt(fasta=outDir+os.path.basename(args.fasta).rsplit(fastaString,1)[0]+domain, domain=domain, out=outDir+'placeIt/',log=logFile)
        runTallyPathways(fasta=outDir+os.path.basename(args.fasta).rsplit(fastaString,1)[0]+domain, domain = domain,out=outDir+'tallyPathways/',enzyme=args.enzyme,cutoff=args.cutoff,log=logFile)

        (enz_edge,koTable)=consolidateOtuGenes(fasta=outDir+'/'+os.path.basename(args.fasta),domain=domain,counts=counts,otuList=otuList,enz=args.enzyme)
        if args.meta=="all":
            print "Combining all the counts for "+domain+". This may take a while depending on the number of samples and OTUs."
            koTableAll.append([[sum(x) for x in zip(*row)] for row in zip(*koTable)])
            enzEdgeAll.append([enz for enz in (enz_edge[0])])
            writeItOut(out=outDir+'/'+domain+'metagenome.txt',koTable=koTableAll[len(koTableAll)-1],enz_edge=enz_edge,samples=samples)

        if args.meta=="both":
            print "Combining all the counts for "+domain+". This may take a while depending on the number of samples and OTUs."
            koTableAll.append([[sum(x) for x in zip(*row)] for row in zip(*koTable)])
            enzEdgeAll.append([enz for enz in (enz_edge[0])])
    if len(domainList)>1:
        (enz_edge,koTable)=combineMetagenomes(enzEdgeAll=enzEdgeAll,koTableAll=koTableAll)
        writeItOut(out=outDir+'/'+'ALLmetagenome.txt',koTable=koTable,enz_edge=enz_edge,samples=samples)


if __name__=="__main__":
    main()
