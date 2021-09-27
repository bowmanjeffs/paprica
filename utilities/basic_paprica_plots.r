## Define the prefix of your input files.  This is whatever you set the -o flag
## to when you ran paprica-combine_results.py.

prefix <- 'your.file.prefix'

## Define some functions to read in the paprica output.

read.edge <- function(prefix, domain){
        tally <- read.csv(paste0(prefix, '.', domain, '.edge_tally.csv'), header = T, row.names = 1) 
        tally <- tally[order(row.names(tally)),]
        tally[is.na(tally)] <- 0
        return(tally)
}

read.unique <- function(prefix, domain){
        unique <- read.csv(paste0(prefix, '.', domain, '.unique_tally.csv'), header = T, row.names = 1)
        unique <- unique[order(row.names(unique)),]
        unique[is.na(unique)] <- 0
        return(unique)
}

read.map <- function(prefix, domain){
        map <- read.csv(paste0(prefix, '.', domain, '.seq_edge_map.csv'), header = T, row.names = 1)
        return(map)
}

read.taxa <- function(prefix, domain){
        taxa <- read.csv(paste0(prefix, '.', domain, '.taxon_map.csv'), header = T, row.names = 1, sep = ',', as.is = T)
        return(taxa)
}

read.data <- function(prefix, domain){
        data <- read.csv(paste0(prefix, '.', domain, '.edge_data.csv'), header = T, row.names = 1)
        return(data)
}

####  ####

tally.bac <- read.edge(prefix, 'bacteria')
tally.arc <- read.edge(prefix, 'archaea')
tally.euk <- read.edge(prefix, 'eukarya')

unique.bac <- read.unique(prefix, 'bacteria')
unique.arc <- read.unique(prefix, 'archaea')
unique.euk <- read.unique(prefix, 'eukarya')

data.bac <- read.data(prefix, 'bacteria')
data.arc <- read.data(prefix, 'archaea')

taxa.bac <- read.taxa(prefix, 'bacteria')
taxa.arc <- read.taxa(prefix, 'archaea')
taxa.euk <- read.taxa(prefix, 'eukarya') 

map.bac <- read.map(prefix, 'bacteria')
map.arc <- read.map(prefix, 'archaea')
map.euk <- read.map(prefix, 'eukarya') 

## join bacteria and archaea

unique <- merge(unique.bac, unique.arc, by = 0, all = T)
tally <- merge(tally.bac, tally.arc, by = 0, all = T)

row.names(unique) <- unique$Row.names
row.names(tally) <- tally$Row.names

unique$Row.names <- NULL
tally$Row.names <- NULL

## convert NAs to 0

unique[is.na(unique)] <- 0
tally[is.na(tally)] <- 0

## eliminate libraries below size threshold

tally.select <- tally[rowSums(tally) > 5000,]
unique.select <- unique[rowSums(unique) > 5000,]

## drop a specific library or libraries

unique.select <- unique.select[grep('SRR14129902.16S.exp.', row.names(unique.select), invert = T),]
tally.select <- tally.select[grep('SRR14129902.16S.exp.', row.names(tally.select), invert = T),]

## eliminate ASVs below a certain abundance threshold

unique.select <- unique.select[,colSums(unique.select) > 10]

## normalize

unique.select.norm <- unique.select/rowSums(unique.select)
unique.select.log10 <- log10(unique.select)
unique.select.log10[unique.select.log10 < 0] <- 0
unique.select.log10 <- unique.select.log10/rowSums(unique.select.log10)

## get taxonomy for ASVs

get.names <- function(domain, map, taxa){
        unique.lab.Row <- map[colnames(unique.select), 'global_edge_num']
        unique.lab.Row <- taxa[unique.lab.Row, 'taxon']
        unique.lab.Row[unique.lab.Row == ""] <- domain
        unique.lab.Row[is.na(unique.lab.Row)] <- domain 
        return(unique.lab.Row)
}

lab.row.bac <- get.names('bacteria', map.bac, taxa.bac)
lab.row.arc <- get.names('archaea', map.arc, taxa.arc)
lab.row.euk <- get.names('eukarya', map.euk, taxa.euk)
lab.row <- cbind(lab.row.bac, lab.row.arc)

## make a heatmap of ASV abundance

## restrict to taxonomy if desired

#selected <- grep('polyedra', unique.lab.Row)

## alternatively restrict to top 50

selected <- order(colSums(unique.select.log10), decreasing = T)[1:50]

heat.col <- colorRampPalette(c('white', 'lightgoldenrod1', 'darkgreen'))(100)

heatmap(t(data.matrix(unique.select.log10[,selected])),
        #Colv = NA,
        scale = NULL,
        col = heat.col,
        labRow = lab.row,
        margins = c(10, 10))

## NMDS plot

library(vegan)
library(oce)

unique.mds <- metaMDS(unique.select.log10, k = 3)

mds.samples <- unique.mds$points

pdf('seasats_mds.pdf',
    width = 6,
    height = 6)

plot(mds.samples[,1], mds.samples[,2],
     ylab = 'Dim 2',
     xlab = 'Dim 1')


## which taxa drive variation in MDS plot?

mds.species <- unique.mds$species

target.asv <- row.names(mds.species)[order(abs(mds.species), decreasing = T)[1:10]]
target.clade <- map.bac[target.asv, 'global_edge_num']
taxa.bac[target.clade, 'taxon']
