###Notice - If you've been redirected here from genome_finder you're in the right place!

###Please go to [releases](https://github.com/bowmanjeffs/paprica/releases) to find the last stable release for download.

###If you'd like to try paprica without needing to download the dependencies we've created a Virtualbox appliance [here](https://dl.dropboxusercontent.com/u/7847544/paprica-demo.ova).  You will need to download the (free) Oracle Virtualbox software and then import the appliance.

#paprica
###PAthway PRediction by phylogenetIC plAcement

A pipeline to conduct a metabolic inference from 16S rRNA gene sequence libraries.  Check out paprica_manual.pdf and paprica_run.sh to get started.  The commands:

```
wget https://github.com/bowmanjeffs/paprica/archive/paprica_v0.XX.tar.gz
tar -xzvf paprica_v0.XX.tar.gz
mv paprica-paprica_v0.XX paprica
cd paprica
chmod a+x paprica_run.sh
./paprica_run.sh test
```

...where XX is the paprica version should get you going and execute a run on the file test.fasta.

###Mac OSX installation tutorial

Installing paprica and all its dependencies on OSX is a little more complicated than on Linux.  You can find a tutorial with detailed instructions [here](http://www.polarmicrobes.org/?p=1477).

###Citation

Please cite paprica as:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula." PloS one 10.8 (2015): e0135868.

###Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes (allow for 8 hours on a 12 core machine with hyperthreading enabled).  You can avoid this by using the provided database and the paprica_run.sh script. 
