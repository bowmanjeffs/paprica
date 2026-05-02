## Don't forget to set your working directory.

setwd('your/working/directory')

library(data.table)

## Define the prefix of your input files.  This is whatever you set the -o flag
## to when you ran paprica-combine_results.py.

prefix <- 'Combined_Datasets_Rarefied'

#### Define some functions to read in the paprica output ####

read.edge <- function(prefix, domain){
    tally <- fread(paste0(prefix, '.', domain, '.edge_tally.csv'), header = T, sep = ',', data.table = F)
    row.names(tally) <- tally$V1
    tally$V1 <- NULL
    tally <- tally[order(row.names(tally)),]
    tally[is.na(tally)] <- 0
    return(tally)
}

read.unique <- function(prefix, domain){
    unique <- fread(paste0(prefix, '.', domain, '.unique_tally.csv'), header = T, sep = ',', data.table = F)
    row.names(unique) <- unique$V1
    unique$V1 <- NULL
    unique <- unique[sort(row.names(unique)),]
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

#### read data files and prepare for analysis ####

tally.bac <- read.edge(prefix, 'bacteria')

data.bac <- read.data(prefix, 'bacteria')

taxa.bac <- read.taxa(prefix, 'bacteria')

map.bac <- read.map(prefix, 'bacteria')

taxa <- taxa.bac
map <- map.bac
tally <- tally.bac

metadata <- read.csv('BP combined metadata.csv')

metadata$file <- paste0(metadata$Run, '.exp.')

metadata <- metadata[which(duplicated(metadata$file) == F),]

row.names(metadata) <- metadata$file

setdiff(row.names(tally), metadata$file)

good.files <- sort(intersect(row.names(tally), row.names(metadata)))

metadata <- metadata[good.files,]

tally <- tally[good.files,]

## convert NAs to 0

tally[is.na(tally)] <- 0

## Eliminate libraries below size threshold.  You should
## pick something that works for your sample set.  In practice
## we try to avoid libraries with fewer than 5000 reads.

#### for each edge determine how often it appears in each amplified region ####

edge.ar <- as.data.frame(matrix(nrow = dim(tally)[2],
                                ncol = length(unique(metadata$Amplified_Region))))

rownames(edge.ar) <- colnames(tally)
colnames(edge.ar) <- unique(metadata$Amplified_Region)

for(edge in rownames(edge.ar)){
  print(edge)
  for(ar in colnames(edge.ar)){
    runs <- row.names(metadata)[which(metadata$Amplified_Region == ar)]
    t <- sum(tally[runs, edge])/sum(tally[runs,])
    edge.ar[edge, ar] <- t
  }
}

edge.ar.pca <- pca(edge.ar, scale = F)
plot(edge.ar.pca)

edge.points <- as.data.frame(edge.ar.pca$CA$u)
ar.points <- as.data.frame(edge.ar.pca$CA$v)

edge.cutoff <- 0.01

good.edges <- row.names(edge.points[which(abs(edge.points$PC1) < edge.cutoff &
                                            abs(edge.points$PC2) < edge.cutoff),])

bad.edges <- row.names(edge.points[which(abs(edge.points$PC1) >= edge.cutoff &
                                            abs(edge.points$PC2) >= edge.cutoff),])

map.bad.edges <- map[which(map$global_edge_num %in% bad.edges),]

taxa.bad.edges <- taxa[bad.edges,]

## instead of removing bad edges try collapsing those edges that belong to problematic taxa

genera.bad.table <- table(taxa.bad.edges$genus)
genera.bad <- row.names(genera.bad.table[genera.bad.table > 1])
genera.bad <- genera.bad[genera.bad != '']

dim(tally)

tally.select <- tally

for(genus in genera.bad){
  temp.taxa <- taxa[taxa$genus == genus,]
  temp.tally <- tally[,row.names(temp.taxa)]
  new.tally <- rowSums(temp.tally)
  tally.select[,genus] <- new.tally
  for(edge in row.names(temp.taxa)){
    tally.select[,edge] <- NULL
  }
  print(dim(tally.select))
}

range(rowSums(tally))

## re-run PCA

edge.ar.2 <- as.data.frame(matrix(nrow = dim(tally.select)[2],
                                ncol = length(unique(metadata$Amplified_Region))))

rownames(edge.ar.2) <- colnames(tally.select)
colnames(edge.ar.2) <- unique(metadata$Amplified_Region)

for(edge in rownames(edge.ar.2)){
  print(edge)
  for(ar in colnames(edge.ar.2)){
    runs <- row.names(metadata)[which(metadata$Amplified_Region == ar)]
    t <- sum(tally.select[runs, edge])/sum(tally.select[runs,])
    edge.ar.2[edge, ar] <- t
  }
}

edge.ar.2.pca <- pca(edge.ar.2, scale = F)
plot(edge.ar.2.pca)

edge.points.2 <- as.data.frame(edge.ar.2.pca$CA$u)

edge.cutoff.2 <- 0.2

bad.edges.2 <- row.names(edge.points.2[which(abs(edge.points.2$PC1) >= edge.cutoff.2 &
                                           abs(edge.points.2$PC2) >= edge.cutoff.2),])

map.bad.edges.2 <- map[which(map$global_edge_num %in% bad.edges.2),]

taxa.bad.edges.2 <- taxa[bad.edges.2,]


##

write.csv(tally.select, file = 'test_edge_tally.csv', quote = F)

#!!! probably need to remove small libraries here

## Subsample to the size of your smallest library.  The alternative is to move
## to relative abundance, but that gets weird if you want to apply a log
## normalization later on.

tally.select.ra <- tally.select/rowSums(tally.select) # relative abundance
tally.select.sub <- tally.select
sub.size <- floor(min(rowSums(tally.select.sub))) # figure out which samples were dominated by dropped edges and remove

for(row in 1:length(row.names(tally.select.ra))){
    temp <- table(sample(colnames(tally.select), size = sub.size, replace = T, prob = tally.select.ra[row,]))
    temp <- as.data.frame(temp)
    row.names(temp) <- temp$Var1
    temp$Var1 <- NULL
    tally.select.sub[row,row.names(temp)] <- temp[row.names(temp), 'Freq']
}

#tally.select.sub <- rrarefy(tally.select, 5000)

## Normalize. We find that a typical log10 normalization often works fine
## for community structure data, but of course this is data and analysis
## dependent.

tally.select.log10 <- log10(tally.select.sub)
tally.select.log10[tally.select.log10 < 0] <- 0
tally.select.log10 <- tally.select.log10/rowSums(tally.select.log10)

## get taxonomy for ASVs

get.names <- function(map, taxa){
  tally.lab.Row <- as.data.frame(map[colnames(tally.select), 'global_edge_num'])
    row.names(tally.lab.Row) <- colnames(tally.select)
    colnames(tally.lab.Row) <- c('global_edge_num')
    tally.lab.Row$taxon <- taxa[tally.lab.Row$global_edge_num, 'taxon']
    return(tally.lab.Row)
}

lab.row <- get.names(map, taxa)
lab.row$taxon[grep('phylum_reps', lab.row$global_edge_num)] <- 'Bacteria'

#### make a heatmap of ASV abundance ####

## OPTIONAL: Select a specific taxonomy

get.asvindex.by.taxa <- function(map, taxa, target_taxa, rank){
  target.edges <- row.names(taxa)[which(taxa[,rank] == target_taxa)]
  target.asvs <- row.names(map)[which(map$global_edge_num %in% target.edges)]
  selected <- which(colnames(unique.select) %in% target.asvs)
  return(selected)
}

selected <- get.asvindex.by.taxa(map.bac, taxa.bac, 'Cyanobacteria', 'phylum')

## Alternatively restrict to top 50:

selected <- order(colSums(tally.select.log10), decreasing = T)[1:50]

heat.col <- colorRampPalette(c('white', 'lightgoldenrod1', 'darkgreen'))(100)

side.colors <- metadata$Amplified_Region

library(paletteer)

for(i in 1:length(unique(metadata$Amplified_Region))){
  study = unique(metadata$Amplified_Region)[i]
  study.col = paletteer_d("colorBlindness::paletteMartin")[i]
  side.colors[which(side.colors == study)] <- study.col
}

tally.dist <- vegdist(tally.select.log10)
tally.dend <- as.dendrogram(hclust(tally.dist))

library(gplots)

heatmap.2(t(data.matrix(tally.select.log10[,selected])),
          trace = 'none',
          Colv = tally.dend,
          scale = NULL,
          col = heat.col,
          #labRow = lab.row[row.names(tally.select.log10)[selected], 'taxon'],
          margins = c(10, 10),
          ColSideColors = side.colors)

#### PCA plot ####

tally.pca <- pca(tally.select.log10)

plot(tally.pca)

tally.pca.points <- as.data.frame(tally.pca$CA$u)

plot(tally.pca.points[,1], tally.pca.points[,2],
     type = 'n')

points(tally.pca.points[which(metadata$Amplified_Region == 'V1_V2'), 1],
       tally.pca.points[which(metadata$Amplified_Region == 'V1_V2'), 2])

#### NMDS plot ####

library(vegan)
library(oce)

tally.mds <- metaMDS(tally.select.log10, k = 10)

mds.samples <- unique.mds$points

## Color points according to potential community growth rate.

point.cols <- terrain.colors(length(mds.samples[,1]))[as.numeric(cut(data.bac.select$gRodon.d.mean,breaks = length(mds.samples[,1])))]

plot(mds.samples[,1], mds.samples[,2],
     ylab = 'Dim 2',
     xlab = 'Dim 1',
     pch = 19,
     col = point.cols,
     cex = 2)

### which taxa drive variation in MDS plot? ###

scaling.factor <- 0.5 # data dependent

mds.species <- unique.mds$species

## Select top 50 contributors to first dimension.

target.asv <- row.names(mds.species)[order(abs(mds.species[,'MDS1']), decreasing = T)[1:50]]
target.clade <- map.bac[target.asv, 'global_edge_num']
target.taxa <- taxa.bac[target.clade, 'taxon']

arrows(0,0,
      (mds.species[target.asv,'MDS1'] * scaling.factor),
       (mds.species[target.asv,'MDS2'] * scaling.factor))
