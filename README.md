#PAPRICA
###PAthway PRediction by phylogenetIC plAcement

Pipeline to conduct a metabolic inference from 16S rRNA gene sequence libraries.

Check out PAPRICA_manual.pdf and paprica.sh to get started.

###Overview

PAPRICA conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach, however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.
PAPRICA uses genes shared between the members of all clades on a reference tree to determine what genes are likely to be associated with a phylogenetically placed read.  Because we are most interested in microbial function, PAPRICA presents this information in the form of metabolic pathways predicted for each point of placement.
PAPRICA was designed to use a significant amount of resources up front to construct a database (allow for 36 hours on a 12 core machine with hyperthreading enabled).  Successive analyses are comparatively lightweight.  Because metabolic pathways are predicted for each point of placement on the fly, but are then available for all subsequent analyses, sample run time decreases significantly for each new library.  If all pathways have already been predicted for the reads contained in a library runtime is just a few seconds. 
