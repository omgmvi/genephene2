# genephene2
TO DO: WRITE THE README
        copy the folder - this time proper structure and checks
        change the GENOME folder to the genome_metadata
        change Taxonomy folder to the base
        create Genome folder and all its structure


Introduction : GenePhene2 is meant to be the continuation of the project GenePhene from Prof. Christopher Quince and its collaborators.

The objective is to elaborate a machine learning model able to predict bacteria and archaea metabolic phenotypes from their genomic functional annotation.
More specifically, many bacterias and archaeas cannot be cultivate using standard laboratory technique and therefore knowledge of the phenotypic characteristics are mostly unknown.
On the other hand modern sequencing technologies and the bioinformatic tools associated has facilitated the recovery and reconstruction of 'metagenomic assembled genomes". 
That would mean that we are able to observe the genomes of unknown microorganims but our knowledge of their activities is only provided by comparison with other siblings taxa through RNA16S or similar barcoding.

The approach to be taken here is to link the functionally annotated genomes, through ORF finding and Ortholog detection algorithms, to more or less specific phenotypes such like 'aerobe' or 'catalase activity'. This linkage is being built on machine learning algorithms for supervised classification where the phenotypes play the role of predicted variables and the Ortholog gene name the role of predictors. To advance the current tools we plant to transform the annotated gene names in a continuous numberical variable using NeuroLinguistic programming tools like Word2Vec and Doc2Vec instead of the so called bag of words approximation.


# CODE ORGANIZATION

The project need the building of two types of databases, the phenotypic and the genomic ones. The phenotypic databases are set up as a single table with phenotypic traits in columns and microorganism names in rows. The cells are filled with either numerical values for temperature or pH or categorical values like yes/no/NA. The genomic databases require some middle tables to be built but the final tables for using at modelling stages are single tables with the name of the gene on columns and microorganisms on rows. The cells are filled with gene count (absolute frequency). *That information may change when the Word2Vec models are fitted*. Three different genomes are being constructed for each type of ortholog genes, KEGG orthologs, Cluster of Ortholog Genes (and EGGNOG) and the protein family PFAM.

The database and their manipulation code are split in folders where for each database all the code to build them up is within the folder (with the exception of the tax2ID metadata see below).

Description of folders / databases

./GenePhene2/
    /Databases : It refer to the phenotypic databases, each one correspond to a publication
        /FAPROTAX
        /Hiding
        /GenePhene2
        /MDB
        /Midas3
    /Genomes : It contain  genomic information metadata, at the base level there are the TaxID2Genome.txt files that is metadata to connect the microorg. name to its genome
        /FAPROTAX
        /Hiding
        /GenePhene2
        /MDB
    /Taxonomy : it contains the taxonomic metadata, it really only is a translation between microorganism name in the phenotype database and the Taxon ID in NCBI to find its Reference Genome
        /FAPROTAX
        /Hiding
        /GenePhene2
        /MDB
        /scripts
    /Models : It contain the information to make the models particularly the scripts to generate the JSON config files
        /scripts
        /config.files
        /GenePhene2/
            data.files
            model.results/

## Pipeline
Shortly, the pipeline has to start with the phenotype database and 'join' the microorganism name with the name2taxid and this one with the taxid2genome and lastly with the KEGG_summary.tsv or cogs_summary.tsv or pfamA_summary.tsv. Graphically,

Databases/FAPROTAX/FAPROTAX_1.2.4/FAPROTAX.tsv <--> Taxonomy/FAPROTAX/SpeciesTaxID.tsv <--> Genomes/FAPROTAX/TaxID2Genome_metadata.tsv <--> (GENOME FILE STILL NOT CREATED)
Databases/GenePhene2/Metabolic_features.tsv <--> Taxonomy/GenePhene2/SpeciesTaxID.tsv <--> Genomes/GenePhene2/TaxID2Genome_metadata.tsv <--> Databases/GenePhene2/Odin_data/cogs_summary.tsv (Still not Structured)


