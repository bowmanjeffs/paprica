#!/usr/bin/env Rscript

## This script implements the gRodon program as described in Weissman et al, 2021, PNAS.

library(gRodon)
library(Biostrings)

args <- commandArgs(trailingOnly=TRUE)
domain <- args[2]
ref_dir <- args[1]
path <- paste(ref_dir, domain, 'cds', sep = '/')

genomes <- sort(list.files(path, pattern = 'cds.fasta', full.names = F))

output <- as.data.frame((matrix(nrow = length(genomes), ncol = 8)))
row.names(output) <- matrix(unlist(strsplit(genomes, split = '.cds.fasta')), ncol = 1)[,1]
colnames(output) <- c("CUBHE", "ConsistencyHE", "CPB", "FilteredSequences", "nHE", "d", "LowerCI", "UpperCI")

for(genome in row.names(output)){
  print(paste('Estimating growth rate for', genome))

  try({
    genes <- readDNAStringSet(paste0(path, '/', genome, '.cds.fasta'))
    highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
    out <- as.data.frame(predictGrowth(genes, highly_expressed))
    output[genome,] <- out[1,]}, silent = T)
}

write.csv(output, file = paste(ref_dir, domain, 'gRodon_growth_estimates.csv', sep = '/'), quote = F, row.names = T)

