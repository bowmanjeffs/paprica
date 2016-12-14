#!/bin/bash

#### These are the critical steps for building the PAPRICA database.  This is not necessary
#### to use PAPRICA, but will provide you with some added flexibility and the ability
#### to work directly with the PGDBs.  Be sure to check the beginning of each Python script
#### for user setable variables.  Depending on your system this script will take a substantial
#### amount of time to run and the PGDBs will take a substantial amount of space.  On a 24 core
#### box it takes roughly 8 hours to get all the genomes downloaded and the database built.
#### The PGDBs take up about 100 Gb of space.

pgdb_dir=~/ptools-local/pgdbs/user/
domain=$1
ref_dir=ref_genome_database

## Select gene based on domain.

if domain=eukarya;then
	gene=18S
else
	gene=16S
fi

## 1. download genomes, combine elements, extract 16S

paprica-make_ref.py -ref_dir $ref_dir -download T -domain $domain -cpus 2 &&

## 2. make a reference package from 16S or 18S

paprica-place_it.py -ref_dir $ref_dir -ref combined_$gene.$domain.tax -domain $domain -cpus 2 &&

## 3. run test.bacteria.fasta 

paprica-place_it.py -ref_dir $ref_dir -query test.$domain -ref combined_$gene.$domain.tax -domain $domain -splits 1 &&

## 4. build the reference database.

paprica-build_core_genomes.py -ref_dir $ref_dir -pgdb_dir $pgdb_dir -tree test.$domain.combined_$gene.$domain.tax.clean.align.phyloxml -domain $domain