import sys
import re
import argparse
import os
parser=argparse.ArgumentParser()
parser.add_argument("-i","--fasta")
parser.add_argument("-c","--counts")
parser.add_argument("-o","--outputFile")
parser.add_argument('-f', action='store_true')

args=parser.parse_args()

fasta=args.fasta
if os.path.isfile(args.outputFile)==False or args.f:
	sys.stdout = open(args.outputFile, 'w')
elif os.path.isfile(args.outputFile)==True and not args.f :
	print "File already exists"
	sys.exit()

#{fasta filename}.{domain}.unique.seqs
#this is the map of sequence identifiers (or OTU #) to edge
name= fasta.strip("fasta")
seq_edge="bacteria.unique_seqs.csv"
seq_edge_map= name+seq_edge
seq_edge_dict={}

#this creates a map between the sequence and the edge within paprica
with open(seq_edge_map,'rb') as csvfile:
	csvfile.readline()
	for row in csvfile:
		#col1.append(row.strip('"').split("\t",row.count(row))[0].strip('"'))
		#sys.stdout.write((row.strip('"').split("\t",row.count(row))
		temp_array=row.strip('\n').split(',')
		seq_edge_dict[temp_array[3]]=temp_array[2]
		#print float(temp_array[4])		


#for seq,edge in seq_edge_dict.items():
#	print seq,edge

#{fasta filename}.{domain}.ec e.g. seq.bacteria.ec uses seq.fasta for bacteria
#this is the map of mapped edges to enzyme commision numbers

edge_ec="bacteria.ec.csv"
edge_ec_map=name+edge_ec

#this is the map of total counts of EC's identified within the loaded sequence data as mapped per edge. As multiple sequences map to the same edge, you will eventually have to divide by the number of sequences that mapped to the edge to get EC per edge.
with open(edge_ec_map,'rb') as csvfile:
	#string=re.compile('"[A-Za-z0-9.,]*"')
	#for row in csvfile:
	#	print string.findall(row)
		#print row
	new_list=[row.strip('\n').split(',')for row in csvfile]

#print new_list

ec_edge=zip(*new_list)
#print ec_edge
#Create a dictionary of edges to their EC
ec_edge_dict={}
for i in range(0,len(ec_edge)):
	ec_edge_dict[ec_edge[i][0]]=i
count=ec_edge[ec_edge_dict.get('3230')][1:]

#print ec_edge[0][1]
#print count
#for edge,ec in ec_edge_dict.items():
#	print edge,ec
#print [float(i)*2 for i in count]

#{fasta filename}.{domain}.edge_data e.g. seq.bacteria.edge_data uses seq.fasta for bacteria
#this is the data about mapped edges to enzyme commision numbers
n_edge="bacteria.edge_data.csv"
n_edge_map=name+n_edge

#Get the number of sequences that mapped to an edge for re-normalization of EC/edge data.
n_edge_dict={}
with open(n_edge_map,'rb') as csvfile:
	csvfile.readline()
	for row in csvfile:
		values=row.strip('\n').split(',')
		#print values[0]
		
		n_edge_dict[values[0]]=values[1]
		
#for keys,values in n_edge_dict.items():
#	print keys,values

#get the count data of sequences from a tab-delimited table with rownames as identical identifiers created for the fasta input to paprica. 

with open(args.counts, 'rb') as seqtable:
	samples=seqtable.readline().strip('\n').split('\t')
	sample_ec=[[0]*len(new_list[1:]) for i in range(len(samples))]
	#print len(sample_ec[0])
	for row in seqtable:
		data=row.strip('\n').split('\t')

		edge=seq_edge_dict.get(data[0].strip('""'))
		#print edge

		#print (n_edge_dict.get(edge))
		n_edges=float(n_edge_dict.get(edge))
		#print n_edges.strip('""')
		#print N_sixteen
		edge_counts=[float(i) for i in ec_edge[ec_edge_dict.get(edge)][1:]]
		seq_counts=[float(i) for i in data[1:]]
		#print seq_counts
		for i in range(0,len(samples)):	
			sample_ec[i]=[x+y for x,y in zip(sample_ec[i],[seq_counts[i]* j /n_edges for j in edge_counts])]
			


ec_sample=zip(*sample_ec)
sys.stdout.write("\t".join(samples)+"\n")
#print ec_sample
i=1


for row in ec_sample:
	#print row
	
	sys.stdout.write(ec_edge[0][i]+"\t"+"\t".join(str(j) for j in row)+"\n")
	i=i+1

