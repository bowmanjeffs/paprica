#!/bin/bash

#### These are the critical steps for using paprica if you are use the provided database (a.k.a. ref_genome_database) or
#### have already built it using paprica_build.sh.  Used in this way paprica is nice and lightweight, but you won't have 
#### access to the PGDBs if you want to do something more sophisticated than just tally up the number of metabolic pathways
#### that have been inferred.

#### If you have a large number run this script in a loop (see the manual for an example).  Because the bottleneck is alignment,
#### and infernal is parallelized, it is best not to run samples in parallel.

#### Execute this script as ./paprica_run.sh [query] [domain].

query=$1
domain=$2

## 1. phylogenetic placement of query reads

paprica-place_it.py -ref_dir ref_genome_database -query $query -ref combined_16S.$domain.tax -splits 1 -domain $domain

## 2. find pathways and other information associated with edges.  if you subsampled in the previous step (i.e. with -n) your input
##    file is $query.sub.combined_16S.tax.clean.align.csv

paprica-tally_pathways.py -ref_dir ref_genome_database -i $query.combined_16S.$domain.tax.clean.align.csv -o $query.$domain -cutoff 0.5 -domain $domain