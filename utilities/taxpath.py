"""
@author: Avishek Dutta (avishekdutta14@gmail.com) 
@requires: python3 and Pandas

Requires all the .pathways.csv files from each samples and also require the
taxa_abundance.csv (output from taxa_view.py or taxa_view_v2.py).Copy all the .pathways.csv (output of paprica) file 
to a new folder in which only the pathways.csv files and taxa_abundance.csv 
is present and run the script "python taxpath.py".

N.B.: Please note that the pathway abundance for a particular taxon reported in the output is 
the cumulative abundance of the pathways present across all the samples used for the analysis.

"""
import pandas as pd
import os

path = '.'
files_in_dir = [f for f in os.listdir(path) if f.endswith('.pathways.csv')]
for filename in files_in_dir:
    prefix, ext = os.path.splitext(filename)
    if ext.lower() != '.csv':
        continue
    # Load the data into a dataframe
    df = pd.read_csv(filename, 
                               header=None, 
                               index_col=None, 
                               parse_dates=False)
    df_transposed = df.T
    # Save to a new file with an augmented name 
    df_transposed.to_csv(prefix+'_T'+ext, header=False, index=False)


#merging all the transposed file
path = '.'
files_in_dir = [f for f in os.listdir(path) if f.endswith('.pathways_T.csv')]
for filenames in files_in_dir:
    edge_data = pd.read_csv(filenames)
    edge_data.to_csv('merged_pathways_data.csv', mode='a')


#removal of all ec_T.csv file 

dir_name = "."
test = os.listdir(dir_name)

for item in test:
    if item.endswith("pathways_T.csv"):
        os.remove(os.path.join(dir_name, item))


#reading the file and renaming the column for finding the index later
df = pd.read_csv('merged_pathways_data.csv')

df1 = df.drop(['Unnamed: 0'], axis=1)

df1.rename(columns={'Unnamed: 0.1': 'OTU'}, inplace=True)

df2 = df1[df1.OTU != 'Unnamed: 0']


#df2['OTU'] = df2['OTU'].astype(float)

df3 = df2.astype(float)

#print(df3)

df4= df3.groupby(['OTU']).sum()

df5 = pd.read_csv('taxa_abundance.csv')

result = df5.merge(df4,on='OTU',how='left')

result.to_csv('taxpath_view.csv', index = False)

os.remove ('merged_pathways_data.csv')

print("The output is in taxpath_view.csv file")
