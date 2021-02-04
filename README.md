![GitHub repo size](https://img.shields.io/github/repo-size/bowmanjeffs/paprica)
![GitHub top language](https://img.shields.io/github/languages/top/bowmanjeffs/paprica)
![GitHub stars](https://img.shields.io/github/stars/bowmanjeffs/paprica?style=social)
<img src="https://github.com/bowmanjeffs/paprica/blob/master/bacteria_terminal_path_distribution.png" alt="Metabolic pathways in bacterial genomes in current paprica database" width="400" align='right'>
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

## Citations

**Please cite paprica as**:

Bowman, Jeff S., and Hugh W. Ducklow, 2015. [Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0135868). *PloS one*, 10.8, e0135868.

**Paprica relies on a number of tools and databases developed by other groups.  If you cite paprica please be sure to also cite these publications:**

Barbera, P., Kozlov, A.M., Czech, L., Morel, B., Darriba, D., Flouri, T. and Stamatakis, A., 2019. [EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences](https://academic.oup.com/sysbio/article/68/2/365/5079844?login=true). *Systematic biology*, 68(2), pp.365-369.

Czech, L., Barbera, P. and Stamatakis, A., 2020. [Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data](https://academic.oup.com/bioinformatics/article/36/10/3263/5722201?login=true). *Bioinformatics*, 36(10), pp.3263-3265.

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C., Burgaud, G., De Vargas, C., Decelle, J. and Del Campo, J., 2012. [The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy](https://academic.oup.com/nar/article/41/D1/D597/1064851?login=true). *Nucleic acids research*, 41(D1), pp.D597-D604.

Haft, D.H., DiCuccio, M., Badretdin, A., Brover, V., Chetvernin, V., O’Neill, K., Li, W., Chitsaz, F., Derbyshire, M.K., Gonzales, N.R. and Gwadz, M., 2018. [RefSeq: an update on prokaryotic genome annotation and curation](https://academic.oup.com/nar/article/46/D1/D851/4588110?login=true). *Nucleic acids research*, 46(D1), pp.D851-D860.

Karp, P.D., Midford, P.E., Billington, R., Kothari, A., Krummenacker, M., Latendresse, M., Ong, W.K., Subhraveti, P., Caspi, R., Fulcher, C. and Keseler, I.M., 2021. [Pathway Tools version 23.0 update: software for pathway/genome informatics and systems biology](https://academic.oup.com/bib/article/22/1/109/5669859?login=true). *Briefings in bioinformatics*, 22(1), pp.109-126.

Nawrocki, E.P. and Eddy, S.R., 2013. [Infernal 1.1: 100-fold faster RNA homology searches](https://academic.oup.com/bioinformatics/article/29/22/2933/316439?login=true). *Bioinformatics*, 29(22), pp.2933-2935.

## Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes.  You can avoid this by using the provided database. 
