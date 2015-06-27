###Genome Finder
**Introduction**
Genome Finder is a pipeline to conduct a metabolic inference on a collection of 16S rRNA gene sequences.  Given a set of 16S rRNA gene sequences the pipeline returns a collection of metabolic pathways that are expected.  Genome Finder as designed with the analysis of large 16S rRNA gene sequence libraries in mind, such as those generated by 454 or Illumina sequencing, but is also appropriate for small datasets.
**Requirements**
Genome Finder requires a Linux-like operating environment.  It was developed on Ubuntu.  Genome Finder uses several open source tools in addition to some novel components.  Required tools that should be located in your path are:
1.  Python 2.7
2. Biopython
3. The Python joblib module
4. Blast+ (and the 16SMicrobial database)
5. Mothur (and one of the Mothur formatted alignment databases)
6. FastTreeMP
7. pplacer, Guppy, and taxit from the Matsen group
**Components**
Eight Python scripts are provided to create the database and execute the metabolic inference.  The genome_finder.sh script is a template to guide database construction and analysis.  An addition text file, mean_phi_values.txt, creates some data on the available completed genomes that is essential to Genome Finder.
**Basic operation**
The first time genome_finder.sh executes it must construct the database necessary to conduct the metabolic inference.  Counterintuitively it needs to run part of an analysis to construct the database.  Database construction is very time intensive and can take several days, however, subsequent executions are fairly fast.
Each script has, near the top, a set of variables that should be modified by the user to reflect directory pathways and other aspects of the operating environment.
*genome_finder_make_ref.py* (only necessary to construct database)
1.  Downloads all completed genomes from Genbank
2. Finds the 16S rRNA genes in these genomes
3. Compiles 16S rRNA gene copy number with the provided phi values
*genome_finder_place_it.py* (necessary to construct database and for analysis)
1.  Uses the 16S rRNA genes from the completed genomes to build a reference tree with FastTreeMP after alignment with Mothur
2. Uses taxit to convert the alignment and tree into a reference package
3. Executes pplacer to place query reads on the reference tree
*genome_finder_build_core_genomes.py* (necessary to construct database)
1.  Uses the reference tree and blastn to build a consensus genome for each node on the reference tree.
2. Generates statistics for each consensus genome, including size, expected degree of genomic plasticity, and expected 16S rRNA gene copy number
### IN PROGRESS
