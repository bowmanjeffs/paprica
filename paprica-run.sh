#!/bin/bash

#### These are the critical steps for using paprica if you are use the provided database (a.k.a. ref_genome_database) or
#### have already built it using paprica_build.sh.  Used in this way paprica is nice and lightweight, but you won't have 
#### access to the PGDBs if you want to do something more sophisticated than just tally up the number of metabolic pathways
#### that have been inferred.

#### If you have a large number run this script in a loop (see the Wiki for an example).  Because the bottleneck is alignment,
#### and infernal is already parallelized, it is best not to run samples in parallel.

#### Execute this script as ./paprica-run.sh [query] [domain].

query=$1
domain=$2

## Select gene based on domain.

if [ $domain = "eukarya" ];then
	gene=18S
else
	gene=16S
fi

## 1. identify the domain (archaea, bacteria, eukarya) of your input reads

paprica-pick_domain.py -in $query

## 2. phylogenetic placement of query reads

paprica-place_it.py -ref_dir ref_genome_database -query $query.$domain -ref combined_$gene.$domain.tax -splits 1 -domain $domain &&

## 3. find pathways and other information associated with edges.  if you subsampled in the previous step (i.e. with -n) your input
##    file is $query.$domain.sub.combined_$gene.tax.clean.align.csv and your unique file is $query.$domain.sub.$domain.unique.seqs.csv.

paprica-tally_pathways.py -ref_dir ref_genome_database -edpl $query.$domain.combined_$gene.$domain.tax.clean.align.edpl.csv -i $query.$domain.combined_$gene.$domain.tax.clean.align.csv -o $query.$domain -cutoff 0.5 -domain $domain &&

echo "Thanks for using paprica!  Please be sure to read through the manual, and check out the tutorials at www.polarmicrobes.org"
