
"""

@author: Avishek Dutta
@requires: Pandas

This script helps to extract taxonomic hierarchy and also compares the abundance of invidual taxa among different samples.
The input for this script is the output from combine_edge_results.py script.  
Copy the script in paprica output folder and run as:
python taxa_view_v2.py

"""


import pandas as pd
import os
import glob
import csv
import numpy as np

#Reading files having abundance

path = '.'
files_in_dir = [f for f in os.listdir(path) if f.endswith('.edge_tally.csv')]
for filenames in files_in_dir:
    edge_tally = pd.read_csv(filenames, header=None)

transpose = edge_tally.T
transpose.to_csv('transpose_file.csv', header=False, index=False)

#Merging files from different samples

path = '.'
files_in_dir = [f for f in os.listdir(path) if f.endswith('.taxon_map.csv')]
for filenames in files_in_dir:
    taxon_map = pd.read_csv(filenames)
    
#Extracting taxonomy

taxon_map.rename(columns = {list(taxon_map)[0]:'OTU'}, inplace=True)
taxon_map.to_csv("taxonomy.csv", index=False)

#Converting character

transpose_char=pd.read_csv("transpose_file.csv")
transpose_char.rename(columns={'Unnamed: 0': 'OTU'}, inplace=True)
transpose_char['OTU'] = transpose_char['OTU'].astype(int)
transpose_char.to_csv("transpose_file_int.csv", index=False)

#Taxa_abundance 

df1 = pd.read_csv('taxonomy.csv')
df2 = pd.read_csv('transpose_file_int.csv')
result = df2.merge(df1,on='OTU',how='left')
result.to_csv('taxa_abundance1.csv', index=False)

unique= pd.read_csv('taxa_abundance1.csv')
unique.drop_duplicates(subset=None, inplace=True)
unique.to_csv('taxa_abundance.csv', index=False)


group = pd.read_csv('taxa_abundance.csv')
group['phylum'].fillna('Unassigned', inplace=True)
group['class'].fillna('Unassigned', inplace=True)
group['order'].fillna('Unassigned', inplace=True)
group['family'].fillna('Unassigned', inplace=True)
group['genus'].fillna('Unassigned', inplace=True)
group['species'].fillna('Unassigned', inplace=True)
group.to_csv ('Unassigned.csv', index=False)

#Grouping phylum
phylum = group.groupby('phylum').sum()
phylum = phylum.drop('OTU', 1)
phylum.to_csv ('phylum.csv')


#Grouping class
class1 = group.groupby('class').sum()
class1 = class1.drop('OTU', 1)
class1.to_csv ('class.csv')


#Grouping order
order = group.groupby('order').sum()
order = order.drop('OTU', 1)
order.to_csv ('order.csv')


#Grouping family
family = group.groupby('family').sum()
family = family.drop('OTU', 1)
family.to_csv ('family.csv')

#Grouping genus
genus = group.groupby('genus').sum()
genus = genus.drop('OTU', 1)
genus.to_csv ('genus.csv')

#Grouping species
species = group.groupby('species').sum()
species = species.drop('OTU', 1)
species.to_csv ('species.csv')


#Deleting temporary files

os.remove ('transpose_file.csv')
os.remove ('taxonomy.csv')
os.remove ('transpose_file_int.csv')
os.remove ('taxa_abundance1.csv')
os.remove ('Unassigned.csv')
print ('The output of the file is present in taxa_abundance.csv')
