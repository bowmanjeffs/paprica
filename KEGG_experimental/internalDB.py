from Bio import Phylo

import subprocess
import os
import re
import sys
import json

ref_dir = 'ref_genome_database/'

def getTree(domain):
    cwd=os.getcwd()+'/'
    ref_dir_domain=ref_dir+domain+'/'
    json_data=open(ref_dir_domain+'taxit/CONTENTS.json').read()
    data=json.loads(json_data)
    treeFile=data["files"]["tree"]
    print "getting tree: "+treeFile
    tree=Phylo.read(ref_dir_domain+'taxit/'+treeFile,'newick')
    for idx,clade in enumerate(tree.find_clades(order="postorder")):
        clade.confidence= idx
    return(tree)






def getTerminalData(domain,gene,tree):
    print "collecting terminal data"
    cwd=os.getcwd()+'/'
    ref_dir_domain=ref_dir+domain+'/'
    assemblies = []
    
    for clade in tree.get_terminals():
            clade_number = int(clade.confidence)
            assembly = clade.name
            assemblies.append(assembly)

    gene_dict={}
    with open(ref_dir+'paprica'+gene+'_DB.txt') as tabfile:
        for row in tabfile:
            info=row.strip('\n').split('\t')
            if info[0] in gene_dict and info[0] in assemblies:
                gene_dict[info[0]].append(info[1:])
            elif info[0] in assemblies:
                gene_dict[info[0]]=[info[1:]]
    return(gene_dict)

def buildTerminal(domain,gene,gene_dict,tree,newDB):
    print "building database for terminal edges"
    #build terminal nodes data 
    for clade in tree.get_terminals():
        edge = int(clade.confidence)
        org=gene_dict.get(clade.name)
        print clade.name, edge
        for gene in org:
            newDB.write(str(edge)+'\t'+gene[0]+'\t'+str(gene[1])+'\t1\t1\n')



def buildInternal(domain,gene,gene_dict,tree,newDB):
    print "building data for internal edges"
    i_clade=0
    for clade in tree.get_nonterminals():
        i_clade += 1
        print 'collecting data for internal node: '+ str(clade.confidence)+'. Clade '+str(i_clade)+' of '+str(len(tree.get_nonterminals()))
    
        ## Data on the clade that you want later.
    
        edge = int(clade.confidence)
        ## Iterate across all terminal nodes in clades to get the corresponding
        ## assemblies.
    
        ntip = len(clade.get_terminals())
        clade_members = []
        gene_list=[]
        count_list=[]
        for tip in clade.get_terminals():
        
            name = tip.name
            cladeNumber=int(tip.confidence)       
            clade_members.append(tip.name)
            org= gene_dict.get(tip.name)
            tmp_gene=[]
            tmp_count=[]
            for values in org:		
                tmp_gene.append(values[0])
                tmp_count.append(values[1])
                count_list.append(tmp_count)
                gene_list.append(tmp_gene)
        mutual_gene= set().union(*gene_list[:])
        mutual_gene_dict=dict.fromkeys(mutual_gene,0)
        mutual_obs_dict=dict.fromkeys(mutual_gene,0)
        mutual_count_dict=dict.fromkeys(mutual_gene,0)
        for i in range(0,len(gene_list)):
            for j in range(0,len(gene_list[i])):
                if gene_list[i][j] in mutual_gene_dict:
                    previous_gene=mutual_gene_dict.get(gene_list[i][j])
                    previous_count=mutual_count_dict.get(gene_list[i][j])
                    mutual_count_dict[gene_list[i][j]]=previous_count+1
                    mutual_gene_dict[gene_list[i][j]]=previous_gene+int(count_list[i][j])
            for gene,count in mutual_count_dict.items():
                mutual_obs_dict[gene]=float(count/float(len(gene_list)))
                previous_gene=float(mutual_gene_dict.get(gene))
                mutual_gene_dict[gene]=float(previous_gene/float(count))
                newDB.write(str(edge)+'\t'+gene+'\t'+str(mutual_gene_dict.get(gene))+'\t'+str(mutual_obs_dict.get(gene))+'\t'+str(count)+'\n')
def main():
    for domain in ["archaea","bacteria"]:
        for gene in ["EC","KO"]:
            cwd=os.getcwd()+'/'
            ref_dir_domain=ref_dir+domain+'/'
            newDB=open(ref_dir_domain+domain+gene+'DB.txt','w')
            newDB.write('edge_num\tgene\tmean_count\tmean_probs\tabs_occur\n')
            tree=getTree(domain)
            gene_dict=getTerminalData(domain,gene,tree)
            buildTerminal(domain,gene,gene_dict,tree,newDB)
            buildInternal(domain,gene,gene_dict,tree,newDB)

if __name__ =="__main__":
    main()
