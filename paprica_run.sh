#### These are the critical steps for using paprica if you are use the provided database (a.k.a. ref_genome_database) or
#### have already built it using paprica_build.sh.  Used in this way paprica is nice and lightweight, but you won't have 
#### access to the PGDBs if you want to do something more sophisticated than just tally up the number of metabolic pathways
#### that have been inferred.

#### Be sure to check the beginning of the two python scripts for variables that you need to set.  If you have a large number
#### of samples to evaluate it is trivial to either 1) parallelize this script using parallel (watch for memory useage!) or 2) run
#### it in a loop, probably using ls *fasta to pass arguments to the loop.  Because the bottleneck is alignment, and infernal
#### is parallelized, it probably makes more sense to use a loop.

#!/bin/bash

query=$1

## 1. phylogenetic placement of query reads

python paprica_place_it.py -query $query -ref combined_16S.tax -splits 1

## 2. find pathways and other information associated with edges.  if you subsampled in the previous step (i.e. with -n) your input
##    file is $query.sub.combined_16S.tax.clean.align.csv

python paprica_tally_pathways.py -i $query.combined_16S.tax.clean.align.csv -o $query