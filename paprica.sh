#!/bin/bash

#### these are the critical steps for using paprica ####
#### be sure to check the beginning of each script for user editable variables ####

## currently set to run with provided test.fasta, after initial
## setup replace test with the name of your query file, excluding
## the extension
query=test

## you may need to update the paprica version number in this script
version=0.1

## make a reference directory, make sure it is the same as called
## for in the header of each script
mkdir ref_genome_database_v1

## 1. download genomes, combine elements, extract 16S
python paprica_make_ref_v${version}.py

## 2. make a reference package from 16S
python paprica_place_it_v${version}.py combined_16S

## 3. run mock analysis 
python paprica_place_it_v${version}.py $query combined_16S

## 4. get csv of placements
guppy to_csv --point-mass --pp -o $query.csv $query.combined_16S.pplacer.filter.jplace

## 5. get a fat tree with node numbers
guppy fat --point-mass --pp --node-numbers -o $query.fat $query.combined_16S.pplacer.filter.jplace

## 6. build the reference database - you need the node numbers from the previous step for this
python paprica_build_core_genomes_v${version}.py $query.fat

## 7. find genomes represented in query
python paprica_get_genomes_v${version}.py $query.csv $query

## 8. build PGDB
chmod a+x generate_pgdbs.sh
./generate_pgdbs.sh

## 9. tally pathways
python tally_pathways.py $query

#### Other Notes ####
## steps 1, 2, and 6 only need to be run the first time, or when it is desirable to update the database
## databases are NOT backwards compatible, i.e. node numbering will be different in each rebuild
## step 6 takes a long time, depending on how many cpus are available.  expect 24 hours for 24 cpus.
## if you have gnu parallel installed (recommended) run step 8 as parallel < generate_pgdbs.sh