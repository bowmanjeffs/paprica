#!/bin/bash

#### these are the critical steps for using genome finder ####

## 1. download genomes, combine elements, extract 16S
python genome_finder_make_ref.py

## 2. make a reference package from 16S
python genome_finder_place_it.py combined_16S

## 3. run first analysis 
python genome_finder_place_it.py [query] combined_16S

## 4. get csv of placements
guppy_to_csv --pp -o [query].csv [query].final_combined_16S.pplacer.filter.jplace

## 5. get a fat tree with node numbers
guppy fat --pp --node-numbers -o [query].fat [query].final_combined_16S.pplacer.filter.jplace

## 6. build the reference database - you need the node numbers from the previous step for this
python genome_finder_build_core_genomes.py [query].fat

## 7. build hypothetical metagegenome
python genome_finder_get_genomes.py [query].csv [query]

## 8. build PGDB
pathway-tools -lisp -patho /volumes/hd1/sea_ice_taxonomy/[query].query_genomes/ -disable-metadata-saving &> [query].log

## steps 1, 2, and 6 only need to be run the first time, or when it is desirable to update the database
## databases are NOT backwards compatible, i.e. node numbering will be different