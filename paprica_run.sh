#### These are the critical steps for using paprica if you are use the provided database (a.k.a. ref_genome_database) or
#### have already built it using paprica_build.sh.  Used in this way paprica is nice and lightweight, but you won't have 
#### access to the PGDBs if you want to do something more sophisticated than just tally up the number of metabolic pathways
#### that have been inferred.

#### If you have a large number run this script in a loop (see the manual for an example).  Because the bottleneck is alignment,
#### and infernal is parallelized, it is best not to run samples in parallel.

#### Execute this script as ./paprica_run.sh [query] [domain].  If you are not in the install directory you will need to modify
#### the python commands so that python knows where to look for the scripts.  For example:
#### python /home/user/paprica/paprica_place_it.py -query $query -ref combined_16S.bacteria.tax -splits 1 -domain $2

#!/bin/bash

query=$1
domain=$2
paprica_path=/volumes/hd1/paprica_development_v0.30/paprica/
ref_dir=ref_genome_database

## 1. phylogenetic placement of query reads

python ${paprica_path}paprica_place_it.py -ref_dir $ref_dir -query $query -ref combined_16S.$domain.tax -splits 2 -domain $domain

## 2. find pathways and other information associated with edges.  if you subsampled in the previous step (i.e. with -n) your input
##    file is $query.sub.combined_16S.tax.clean.align.csv

python ${paprica_path}paprica_tally_pathways.py -ref_dir $ref_dir -i $query.combined_16S.$domain.tax.clean.align.csv -o $query.$domain -cutoff 0.5 -domain $domain