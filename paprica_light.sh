#!/bin/bash

#### these are the critical steps for using paprica if you don't need to build the database (a.k.a. the reference directory) ####
#### be sure to check the beginning of each script for user editable variables, particularly the location of the reference directory ####
#### this script requires that you have Pathway-Tools, Gnu parallel, pplacer (and Guppy), mothur (and the silva.seed alignment) ####

version=0.11

query=test

##-1. if you haven't alread done so, download the database with your preferred method, here using wget
wget ftp://ftp.ldeo.columbia.edu/archive/bowmanjs/paprica_database/*tgz

## 0. untar.  it will be fairly large (~30 Gb)

tar -xzvf ref_genome_database_a.tgz
tar -xzvf silva.seed_v119.tgz

### now the analysis starts ###

## 1. phylogenetic placement 
python paprica_place_it_v${version}.py $query combined_16S

## 2. get csv of placements
guppy to_csv --point-mass --pp -o $query.csv $query.combined_16S.pplacer.filter.jplace

## 3. get a fat tree with node numbers.  this is not absolutely necessary but is usually nice to have.
guppy fat --point-mass --pp --node-numbers -o $query.fat $query.combined_16S.pplacer.filter.jplace

## 4. find genomes represented in query
python paprica_get_genomes_v${version}.py $query.csv $query

## 5. build PGDB
chmod a+x generate_pgdbs.sh
parallel < generate_pgdbs.sh

## 6. tally pathways
python paprica_tally_pathways.py $query