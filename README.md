![GitHub repo size](https://img.shields.io/github/repo-size/bowmanjeffs/paprica)
![GitHub top language](https://img.shields.io/github/languages/top/bowmanjeffs/paprica)
![GitHub stars](https://img.shields.io/github/stars/bowmanjeffs/paprica?style=social)
![Metabolic pathways in bacterial genomes in current paprica database](https://github.com/bowmanjeffs/paprica/blob/master/bacteria_terminal_path_distribution.png)
## Notice - If you've been redirected here from genome_finder you're in the right place!  Please take a look at the [Wiki](https://github.com/bowmanjeffs/paprica/wiki) for documentation.  There are also several tutorials available:
[Installing paprica on Mac OSX](http://www.polarmicrobes.org/installing-paprica-on-mac-osx/) (updated) 
[Basic analysis with paprica](http://www.polarmicrobes.org/analysis-with-paprica/) (updated)
[Heatmaps and ordination with paprica output](https://www.polarmicrobes.org/tutorial-basic-heatmaps-and-ordination-with-paprica-output) (updated)
[Annotating metagenomes with paprica-mg](http://www.polarmicrobes.org/tutorial-annotating-metagenomes-with-paprica-mg/) (not yet updated)
[Building the paprica database](http://www.polarmicrobes.org/building-the-paprica-database/) (not yet updated)

## We just completed a major update to paprica that provided a number of improvements.  The documentation and tutorials will lag the code development, but the basic commands and structure of the output are similar.  I encourage you to try the latest code by cloning the repository, and to create an issue if you encounter any problems.

## Alternatively you can go to [releases](https://github.com/bowmanjeffs/paprica/releases) to find the last stable release for download.

## To run paprica you will need:
* [epa-ng](https://github.com/Pbdas/epa-ng)
* [gappa](https://github.com/lczech/gappa)
* [infernal](http://eddylab.org/infernal/)
* [seqmagick](https://fhcrc.github.io/seqmagick/)
* python 3.6 or higher, with packages specified [here](https://github.com/bowmanjeffs/paprica/wiki/1.-Requirements-and-Installation)

Paprica can be run on OSX (see tutorial linked above), or preferably, on Linux or Windows using the Windows Subsystem for Linux (see the linux_install.sh script as a guide).

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
...where XX is the paprica version should get you going and execute a run on the provided file test.fasta for the domain Bacteria.

## Citation

Please cite paprica as:

Bowman, Jeff S., and Hugh W. Ducklow. "Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula." PloS one 10.8 (2015): e0135868.

## Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes.  You can avoid this by using the provided database. 
