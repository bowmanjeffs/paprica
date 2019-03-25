## Notice - If you've been redirected here from genome_finder you're in the right place!  Version 0.4 has a number of new features and improvements.  Please take a look at the [Wiki](https://github.com/bowmanjeffs/paprica/wiki) for further details.  There are also several tutorials available:
[Installing paprica on Mac OSX](http://www.polarmicrobes.org/installing-paprica-on-mac-osx/)  
[Basic analysis with paprica](http://www.polarmicrobes.org/analysis-with-paprica/)  
[Heatmaps and ordination with paprica output](https://www.polarmicrobes.org/tutorial-basic-heatmaps-and-ordination-with-paprica-output)  
[Annotating metagenomes with paprica-mg](http://www.polarmicrobes.org/tutorial-annotating-metagenomes-with-paprica-mg/)  
[Building the paprica database](http://www.polarmicrobes.org/building-the-paprica-database/)  

## Please go to [releases](https://github.com/bowmanjeffs/paprica/releases) to find the last stable release for download.

## If you'd like to try paprica without downloading the dependencies you have two options.
Option 1: We've created a Virtualbox appliance [here](http://www.polarmicrobes.org/extras/paprica-demo.ova) running Ubuntu.  You will need to download the (free) Oracle Virtualbox software and then import the appliance. The appliance isn't updated every time we make an improvement to paprica, so you'll probably want to re-clone the Github repository once you've got the VB up and running.  The log-in information for the VB machine is user: tester, pass: paprica (naturally).

Option 2: We've also created an Amazon Web Service machine instance.  Please read the tutorial located [here](http://www.polarmicrobes.org/paprica-on-the-cloud/).  NOTE: the AWS machine instance is no longer supported, but could be on request.

## paprica
PAthway PRediction by phylogenetIC plAcement

A pipeline to conduct a metabolic inference from 16S rRNA gene sequence libraries.  Check out the [Wiki](https://github.com/bowmanjeffs/paprica/wiki) and tutorials (listed above) to get started.  Once you've installed the depenencies the commands:

```
git clone https://github.com/bowmanjeffs/paprica.git
cd paprica
chmod a+x *py
chmod a+x *sh
./paprica-run.sh test bacteria
```
or
```
wget https://github.com/bowmanjeffs/paprica/archive/paprica_v0.XX.tar.gz
tar -xzvf paprica_v0.XX.tar.gz
mv paprica-paprica_v0.XX paprica
cd paprica
chmod a+x *py
chmod a+x *sh
./paprica-run.sh test bacteria
```
...where XX is the paprica version should get you going and execute a run on the file test.fasta for the domain Bacteria.

## Citation

Please cite paprica as:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula." PloS one 10.8 (2015): e0135868.

## Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes.  You can avoid this by using the provided database. 
