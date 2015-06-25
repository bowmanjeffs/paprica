#!/bin/bash

#### these are the critical steps for using genome finder ####
#### be sure to check the beginning of each script for user editable variables ####
#### you will need to update the version numbers in this script ####

query=file-without-extension

## 1. download genomes, combine elements, extract 16S
python genome_finder_make_ref.py

## 2. make a reference package from 16S
python genome_finder_place_it.py combined_16S

## 3. run first analysis 
python genome_finder_place_it.py $query combined_16S

## 4. get csv of placements
guppy_to_csv --pp -o [query].csv $query.combined_16S.pplacer.filter.jplace

## 5. get a fat tree with node numbers
guppy fat --pp --node-numbers -o $query.fat $query.combined_16S.pplacer.filter.jplace

## 6. build the reference database - you need the node numbers from the previous step for this
python genome_finder_build_core_genomes.py $query.fat

## 7. find genomes represented in query
python genome_finder_get_genomes.py $query.csv $query

## 8. build PGDB
chmod a+x generate_pgdbs.sh
./generate_pgdbs.sh

## 9. tally pathways
python tally_pathways.py

#### Other Notes ####
## steps 1, 2, and 6 only need to be run the first time, or when it is desirable to update the database
## databases are NOT backwards compatible, i.e. node numbering will be different in each rebuild
## step 6 takes a long time, depending on how many cpus are available.  expect 24 hours for 24 cpus.
## if you have gnu parallel installed (recommended) run step 8 as parallel < generate_pgdbs.sh