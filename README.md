![GitHub repo size](https://img.shields.io/github/repo-size/bowmanjeffs/paprica)
![GitHub top language](https://img.shields.io/github/languages/top/bowmanjeffs/paprica)
![GitHub stars](https://img.shields.io/github/stars/bowmanjeffs/paprica?style=social)
![docker build](https://img.shields.io/docker/cloud/build/jsbowman/paprica)
<img src="https://github.com/bowmanjeffs/paprica/blob/master/bacteria_terminal_path_distribution.png" alt="Metabolic pathways in bacterial genomes in current paprica database" width="400" align='right'>
## Notice - If you've been redirected here from genome_finder you're in the right place!  Please take a look at the [Wiki](https://github.com/bowmanjeffs/paprica/wiki) for documentation.  There are also several tutorials available:
[Installing paprica on Mac OSX](http://www.polarmicrobes.org/installing-paprica-on-mac-osx/) (updated) 
[Basic analysis with paprica](http://www.polarmicrobes.org/analysis-with-paprica/) (updated)
[Heatmaps and ordination with paprica output](https://www.polarmicrobes.org/tutorial-basic-heatmaps-and-ordination-with-paprica-output) (updated)
[Annotating metagenomes with paprica-mg](http://www.polarmicrobes.org/tutorial-annotating-metagenomes-with-paprica-mg/) (not yet updated)
[Building the paprica database](http://www.polarmicrobes.org/building-the-paprica-database/) (not yet updated)

## Announcements

### 9 April 2022 - The bacteria and archaea databases have been updated and now contain 11,305 and 346 genomes respectively. We have not yet fully validated the taxonomic assignments returned by the new bacterial trees, so please do due diligence with your own data and report any misassignments.

1 July 2021 - You can now download a [Docker image](https://hub.docker.com/repository/docker/jsbowman/paprica/tags?page=1&ordering=last_updated) for paprica.

If you're not using the Docker image please be sure to update to the most recent version of gappa, as paprica assumes the naming conventions of the most recent version.

## To run paprica you will need:
* [epa-ng](https://github.com/Pbdas/epa-ng)
* [gappa](https://github.com/lczech/gappa)
* [infernal](http://eddylab.org/infernal/)
* [seqmagick](https://fhcrc.github.io/seqmagick/)
* python 3.6 or higher, with packages specified [here](https://github.com/bowmanjeffs/paprica/wiki/1.-Requirements-and-Installation)

Paprica can be run on OSX (see tutorial linked above), or preferably, on Linux or Windows using the Windows Subsystem for Linux (see the [linux_install.sh](https://github.com/bowmanjeffs/paprica/blob/master/linux_install.sh) script as a guide).

If you use Docker or Singularity you can download a Docker image with `docker pull jsbowman/paprica:latest`

## paprica
PAthway PRediction by phylogenetIC plAcement

A pipeline to conduct a metabolic inference from 16S rRNA gene sequence libraries.  Check out the [Wiki](https://github.com/bowmanjeffs/paprica/wiki) and tutorials (listed above) to get started.  Once you've installed the depenencies the following commands will get you started:

```
git clone https://github.com/bowmanjeffs/paprica.git
cd paprica
chmod a+x *py
chmod a+x *sh
./paprica-run.sh test bacteria
```
Alternatively, if you're using the Docker image you can skip installing the dependencies and just do:
```
docker pull jsbowman/paprica:latest
docker run -it jsbowman/paprica
cd /paprica
./paprica-run.sh test bacteria
```

## Citations

**Please cite paprica as:**

Bowman, Jeff S., and Hugh W. Ducklow, 2015. [Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0135868). *PloS one*, 10.8, e0135868.

**Paprica relies on a number of tools and databases developed by other groups.  If you cite paprica please be sure to also cite these publications:**

Barbera, P., Kozlov, A.M., Czech, L., Morel, B., Darriba, D., Flouri, T. and Stamatakis, A., 2019. [EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences](https://academic.oup.com/sysbio/article/68/2/365/5079844?login=true). *Systematic biology*, 68(2), pp.365-369.

Czech, L., Barbera, P. and Stamatakis, A., 2020. [Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data](https://academic.oup.com/bioinformatics/article/36/10/3263/5722201?login=true). *Bioinformatics*, 36(10), pp.3263-3265.

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C., Burgaud, G., De Vargas, C., Decelle, J. and Del Campo, J., 2012. [The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy](https://academic.oup.com/nar/article/41/D1/D597/1064851?login=true). *Nucleic acids research*, 41(D1), pp.D597-D604.

Haft, D.H., DiCuccio, M., Badretdin, A., Brover, V., Chetvernin, V., O’Neill, K., Li, W., Chitsaz, F., Derbyshire, M.K., Gonzales, N.R. and Gwadz, M., 2018. [RefSeq: an update on prokaryotic genome annotation and curation](https://academic.oup.com/nar/article/46/D1/D851/4588110?login=true). *Nucleic acids research*, 46(D1), pp.D851-D860.

Karp, P.D., Midford, P.E., Billington, R., Kothari, A., Krummenacker, M., Latendresse, M., Ong, W.K., Subhraveti, P., Caspi, R., Fulcher, C. and Keseler, I.M., 2021. [Pathway Tools version 23.0 update: software for pathway/genome informatics and systems biology](https://academic.oup.com/bib/article/22/1/109/5669859?login=true). *Briefings in bioinformatics*, 22(1), pp.109-126.

Nawrocki, E.P. and Eddy, S.R., 2013. [Infernal 1.1: 100-fold faster RNA homology searches](https://academic.oup.com/bioinformatics/article/29/22/2933/316439?login=true). *Bioinformatics*, 29(22), pp.2933-2935.

**If you use the gRodon.d estimate of doubling time please also cite:**

Weissman, J.L., Hou, S. and Fuhrman, J.A., 2021. [Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns](https://www.pnas.org/doi/10.1073/pnas.2016810118). *Proceedings of the National Academy of Sciences*, 118(12).

**If you use the paprica-mt module please also cite:**

Bowman, J.S., Van Mooy, B.A., Lowenstein, D.P., Fredricks, H.F., Hansel, C.M., Gast, R., Collins, J.R., Couto, N. and Ducklow, H.W., 2021. [Whole Community Metatranscriptomes and Lipidomes Reveal Diverse Responses Among Antarctic Phytoplankton to Changing Ice Conditions](https://www.frontiersin.org/articles/10.3389/fmars.2021.593566/full). *Frontiers in Marine Science*, 8, p.119.

Li, H. and Durbin, R., 2009. [Fast and accurate short read alignment with Burrows–Wheeler transform](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true). *Bioinformatics*, 25(14), pp.1754-1760.

**If you use the paprica-mg module please also cite:**

Buchfink, B., Xie, C. and Huson, D.H., 2015. [Fast and sensitive protein alignment using DIAMOND](https://www.nature.com/articles/nmeth.3176). *Nature methods*, 12(1), pp.59-60.

## Overview

Paprica conducts metabolic inference on (preferably, but not exclusively, NGS) 16S rRNA gene libraries.  Instead of using an OTU based approach however, it uses a phylogenetic placement approach.  This provides a more intuitive connection between its “hidden state prediction” and library analysis components, and allows resolution at the strain and species level for some spots on the prokaryotic phylogenetic tree.

Paprica uses pathways shared between the members of all clades on a reference tree to determine what pathways are likely to be associated with a phylogenetically placed read.

Paprica was designed to use a significant amount of resources up front to construct a database and draft metabolic models for all available completed genomes.  You can avoid this by using the provided database. 
