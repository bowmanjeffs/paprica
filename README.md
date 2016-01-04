###Notice - If you've been redirected here from genome_finder you're in the right place!  We are currently working through a bug in paprica_build that results in empty dataframes.  The current version of paprica_build is not operational, paprica_run should execute just fine.

#PAPRICA
###PAthway PRediction by phylogenetIC plAcement

A pipeline to conduct a metabolic inference from 16S rRNA gene sequence libraries.  Check out paprica_manual.pdf and paprica_run.sh to get started.  Once you've downloaded the genome_finder directory the commands:

chmod a+x paprica_run.sh

./paprica_run.sh test

Should get you going and execute a run on the file test.fasta.

###Citation

Please cite PAPRICA as:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula." PloS one 10.8 (2015): e0135868.

###Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach, however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes (allow for 8 hours on a 12 core machine with hyperthreading enabled).  You can avoid this by using the provided database and the paprica_run.sh script. 
