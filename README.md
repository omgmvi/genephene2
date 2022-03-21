# **Genephene2** : *A set of models to predict microbe metabolic phenotypes from their annotated genomes*


## Introduction
GenePhene2 is a continuation of the GenePhene project from Prof. Christopher Quince and his collaborators.

The aim of this project is to find appropiate  machine learning model able to predict bacteria and archaea metabolic phenotypes from their genomic functional annotation having having in mind metagenomics experiments. More specifically, it is now possible to take environmental DNA samples from many environments, e.g. soil, water, sludge, air, that can now be sequenced either in long reads, i.e. Oxford Nanopore, or short reads, e.g. Illumina. Those samples contain many species of  bacterias and archaeas that cannot be cultivated using standard laboratory technique and therefore little knowledge of their phenotypic characteristics are known. On the other, hand modern sequencing technologies and associated bioinformatic tools has facilitated the recovery and reconstruction of 'metagenomic assembled genomes",~~i.e. single species genomes reconstructed by assembly contiguous sequences and further filtering of quimeric sequences~~. Therefore, we have access now to new species or strains whosemetabolic activities are only deduced by comparison with their siblings taxa, calculated by RNA16S or similar barcoding.~~The approach of phenotype imputation by phylogenetic association is troublesome in the sense that many strains or new microbes species are actually defined by their differential phenotypic capabilities from their closest siblings. In short, new taxa are defined as those that do something their sister taxa don't.~~ The phylogenetic based phenotypic imputation is the chosen technique used by some popular database-algorithm combination, e.g. FAPROTAX, Tax4Fun, MIDAS, that maps RNA 16S profiles to their imputed metabolic activities.
The approach taken here is linking the functionally annotated genomes, after Open Reading Frame and Ortholog gene imputation, to more or less specific phenotypes such like 'aerobe' or 'catalase activity'. Each phenotype is associated after fitting supervised classification algorithms where  phenotypes play the role of predicted variables and ortholog genes the one of predictors.

To solve several difficulties regarding undetermined matrix (many more genes than individuals and therefore columns are not independent), phylogenetic signal in the gene space and other correlations, we test to transform annotated genomes in continuous numerical variables using the NLP (for Natural Language Processing) algorithm  Doc2Vec instead of the so called Bag of Words approximation. If possible, other _feature engineering_ methods would be tried.

At this stage, several phenotypic databases will be used for training the models and check the performance of them, being them GenePhene ( an small custom database of many metabolic traits from acidogenic taxa for testing purposes), FAPROTAX, Hiding in Plain Sight and MDB (Methanogens Database). The database MIDAS was dismissed due to the low number of species, being most entries genera or families. These databases have been curated to keep only entries corresponding to bacterial with full names and their metabolic related phenotypes. We dismiss taxa levels over species, i.e. genera, families and above, and numerical phenotypes like temperature range or non-metabolic ones like sporulation. Therefore, the traits kept were of type _Consume resource X_, _Produce resource Y_ or more general trtraits like _Catalase activity_. Yet, each database encode a quite different set of traits and a quite non-overlapping number of bacterias.

# CODE ORGANIZATION

The code, at the moment, revolt around a single script, model\_setup.R that perform all the next steps in the analysis. The script and associate 'sender' scripts are found at ~/Genephene2/Models/scripts/

1. Generate Input data structure - single table .
2. Clean data according to rules
3. Stratify sample (if desired).
4. Split data set in validation and training dataset
5. Perform desired model
6. Evaluate model performance

The script was designed with two ideas in mind. On one hand, the number of machine learning models, i.e. logistic regression, neural networks, support vector machines, to try may be big or at least changing over time, so a flexible way to apply different models to each data set, e.g. phenotype, genome, test dataset. The second one is that the input data was going to be changing often and the combinations were going to be big; three orthologs times two types of genome transforms expandable to more transformes times several databases. 
In addition, the script was thought to be running in paralell in a HPC computer and therefore made sense to sacrifice a bit of overhead in setting up the input data structure flexibly (see below) to have a quick and flexible way to try as many phenotypes and models as were possible. The counterpoint is that the model_setup.R script couples all the steps when decoupling them onto a workflow provide more sense and flexibility

The model_setup.R script have a companion script at the folder ~/GenePhene2/Models/templates to show how to send many models for fitting using bash as if it were a job queue ( a for loop). The 'too long;don't read' is that a set of JSON files are need to make it run and how to generate them can be found at this **link** and how to run a model in this **link**.

## Data Input Structure

The machine learning models for phenotype imputation requires a table with one column representing the phenotype to model, typically as a yes/no, and many columns with the 'genes', either orthologs or Doc2Vec variables. At this moment, we have train models for four phenotypic database, named FAPROTAX, Hiding in plain sight, Methanogenes DataBase and a custom one for practice that we have called GenePhene. Each of them come from a different paper with arbitrary formats, so we have translated them to a common text format with columns __|_Genus_|_Species_|_Name_|Phenotype 1|Phenotype 2|...|Phenotype N|:__.
In addition, for each phenotypic data base there are 3 other corresponding tables.
- Taxonomy
    - **|_Species_|TaxID|**
    - Translates microbe names (Genus + Species separated by a single space) to a TaxID (The number has been searched semi-manually at NCBI database but other DBs can be used if decide).
- Metadata
    - **|assembly\_accession|taxid|species\_taxid|_organism\_name_|**)
    - Connects NCBI taxid (notice that there are two columns with TaxIDs at NCBI DB one for species and other for other taxa levels).
- Genome, which come in two flavors, at the moment: 
    - **BagOfWords** genomes. A different file for each ortholog annotation, COG,KEGG and PFam, then the header is:
        - **|_GenomeID_|COG0001|...|COGN|**  OR **|_GenomeID_|K0001|...|KN|** OR **|_GenomeID_|PF0001|...|PFN|** acording to the ortholog database type.
        - NOTE 1: In _a posteriori_ thought, make sense unify all these on GeneN independently the  type of ortholog.
        - NOTE 2: The column GenomeID cannot be re-named in any other way since some part of the code are hardwired to it and this need to be changed soon.
        - Bag of Words genomes  contain each ortholog gene as column and each 'cell' contain the frequency of each ortholog in the genome. Usually there are 0,1 or 2 but there might be up to 22 copies of  the same ortholog).
    - **Doc2Vec** genomes. Again one file  per ortholog type. See the details to compute them below.
        - **|_GenomeID_|D2V0001|...|D2VN|** (independently of which ortholog genome is).
    
    Given this table structure (schema in relational DB terminology) an script can merge (inner join or other joins in SQL terminology) all of them to get a single table with the phenotypes and genomes removing the need for the intermediate databases. Yet, in our script model_setup.R do perform each time this join unnecesarily (due to the history of how this development was done by tinkering rather than proper planning). This join operation was usually not very problematic in performance, since the overhead time was seconds for the databases we were using, with the big exception of a clean operation for multiple genomes.
    In principle, uncoupling the join from the rest of the code is not a problem but it is worthy to notice few advantages of keeping the database as separate tables. On one side, many genome type will be tried, and more may come in the future, so we need to keep the original tables for future joins. In addition, manipulation the current tables provide advantages, for example, the phenotypic databases can be split by phenotype and subset randomly or given any heuristic to have different dataset to perform repetitions or subsets. That is not very relevant as we will see in the next section. Yet, we can subset the taxonomy database to perform the models in subgroups, like just in Gram + bacteria, or just archaea or just using a number of representatives of every group and avoid phylogenetic bioas. Similarly, there are some genome repetition, meaning some microbes have several references genomes or multiple strains, that also bias the database. That repetition can be manipulated using the metadata table (since we will mention later than the data clean is part of the most problematic and slow part of the script).

    So, at the moment, the input for the script is a JSON file with the required field to find, clean and join the data base plus other configurations we'll show in the next paragraphs.
    Link to the methods to compute the join, the cleaning of data and the JSON.

## Data Cleaning
Before joining the table it is neccesary and recomendable to clean data to reduce memory consuming (thinking on possible big data problems). Independently of when it was done, the data cleaning involvedi (NOTE: there some information on the option to set for some of the cleaning, full explanation at the JSON files link below):

- Select phenotype column (remainder: as the script is now, the goal was to make the join with a single phenotype - when this part of the algorithm is decoupled there were no need to filter the column)
- Choose what to do with NA values. Most of the databases we chose for this work were unstructure, having only the name of the bacteria and the list of compound they degrade or produce among other phenotypes. When translating them to tables, only the 'yes' values were present, and the rest where left as missing values. Now, we need to choose what to do with them. We have then the option 'NA_substitution' (see JSON files below) that can be setup either as'None' for doing nothing (but then we only have one class for the FAPROTAX and Hiding databases) or translate to a different symbol with the 'NA2level' option. In the latter case another option 'NA_substitution_value' must be set up with the value we desire for the NA. For the databases that already have phenotypes as yes/no, this value can be a 'no' or a third class.
- Genome cleaning. We found that many microbes have several reference genomes or several straing sequenced. We had to implement an option 'Fix_multiple_genomes' with several possible values like 'identity' to keep the genome as is, 'mean', 'median', 'var','sum' to aggregate the values for each gene, 'pickone' to get one random genome (but at the moment pick just the first). It should be a 1/0 function (to be added if needed) to eliminate the frequency and leave it as presence/abscense

## Trait stratified sampling 
Once the usable data table has been generated/load we need to take the decision of wether balance the sample, make an stratified sampling, or keep the data as is. In the case of choosing to stratify the sample, we use the package caret (as for most of the rest of machine learning operations) that contain the functions downsample and upsample. The first, reduce the number of cases from every class to the that of the minimum, while in the upsample repeat rows from the low number classes to equate to the largest one. There is another option in the package called SMOTE that seems to mix columns from the lower class to make it bigger. We did not include that as an option since it looks confusing to mix genomes (although this can be even helpful to remove phylogenetic bias in the genome)

The option is called 'fix_balance_classes' with options \<downsample|upsample|None\>

## Validation / Training dataset split
    From here onwards, we are using the options at the modelling code (see JSON section below).
    
    NOTE: As we will see in the next section, the script was made for calculating several models for the same dataset. Yet, every call to the model_setup.R will make a different validation dataset and there's no options at the moment to make a list of validation dataset (it is possible with the function but our datasets where too small for this). Then, to repeat a model with different validation/training split one have to call another instance of the split and conversely, one cannot repeat the same validation/trainin split in another instance since we do not save them nor accept them as an option (a future idea to implement). What this script can do is to fit several types of models to the same input data, and all this models will share the validation/trainin data. More clear, every run of the script will have a random split, but if you use a multimodel run, within this run you will be using a single split.

## Machine learning models and their evaluation.

The script model_setup.R was designed with the purpose to apply 1 or several machine learning models to the same dataset. That was motivated by the idea that in the case of big data the overhead of reading the data may be large, although eventually we did not tray such large databases and the overhead was govern by the models itsef and the process of the genome. Therefore, this script is able to perform several models to a single training dataset and evaluate the performance in a common validation dataset. This provide the chance to compare several models using the exact similar data.

The script was designed to be able to add new machine learning models. At the moment, it only has three and coming from the package caret, glmnet - for a logistic regression with elastic net regularizaation - Naive Bayes and XGBoost. Adding new models it is relatively easy by adding few functions in specific places of the script but a modification is needed to flexibilize the model addition without tinkering the code.

Each model is trained on the training dataset using any type of cross-validation or without. To achieve that, input options (see JSON section below) has been set up to account for the type of cross-validation and the number of folds, e.g. leave one out (LOOCV), K-fold with 10-fold, ... and it is responsability of the modelling function to run the cross-validation. Being the script programmed with the package caret in mind, the function _train()_ together with _trainControl()_ take care of this.

The output of the machine learning models from the caret package has not been managed so far, so that in case of using any other package it would be needed to homogenize or make some sort of polymorphism. The current output contain the most fit model parameters together with a full description of model fitting, steps, parameters tried and corresponding metrics, etc.This output is then used by the **evaluation** function to calculate a confusion matrix (since they are all classification models)  and with that set of performance metrics.

The evaluation and the models are saved  (well, at the moment the models are not being saved to avoid hard drive space clutter) at the output folders. A modification required though is to save only the coefficients as a matrix or any recoverable structure that can be used in other packages. Using serialization off-the-shelve solutions in R,  like pickle package in python, are possible but undesired for sharing with other software and also they are bigger than they are desired, since the object output from caret's _train_ funcion contain unnecesary and massive information for its own internal work.

# Data base organization

As explained above, the database are organize around four tables. Albeit being a case of 'technological debt' the tables are not encoded in a relational (SQ)  database manager which would help and accelerate the whole system. Thus, they text files, tab separated (tsv) and with \n as EOL and UTF-8 encoding, were decide to be saved in separated folders per type, e.g. Genome, Metadata, ... and per original phenotype database. Names and folder can be easily changed and rearranged since the information is passed to the programs through JSON files (see below). Yet a short description follows with the files and scripts elaborated for each type.

## Description of folders / databases

../GenePhene2/
    /Databases : _The name would be better **Phenotypes** since it contains each metabolic phenotype databases._
                 Each subfolder refer to the phenotypic database, corresponding to a published one.
        /FAPROTAX
            /FAPROTAX_1.2.4_zip : The original FAPROTAX database. contain the files FAPROTAX.txt with the whole phenotypes in a unstructure format
            /FAPROTAX.tsv : The phenotypic database
            /see.tsv.sh : script to visualize the database in the bash command line
            /get_list_of_phenotype.sh : script to extract the phenotypes
            /get_list_of_speces.sh : script to extract the names of the species
            /parse_faprotax.R : script to parse the original files from FAPROTAX
            /summarize.R : script to calculate the balance of phenotype columns
            /graphicalsummary.R : script to make an ASCII plot for the balance of phenotypic columns
        /Hiding : database from the paper Hiding in plain sight - one of the databases backing BacDive
            /instructions : description of steps to convert from the original database IJSEM_pheno_db_v1.0.txt to ijsem.txt
            /Substrates.tsv : Database of Substrates (Consume X)
            /Metabolic.tsv : Database of several metabolic activities
            /extract_phenotypes.R : script to pass from a manually corrected DB (see instructions) wiht UNICHAR and \n EOL to ijsem.tsv with ranges split as columns and others
            /get_phenotypes.R: Split ijsem.tsv onto two tables one for Substrates and the other for Metabolic activities.
            /extract_species.R : scripts to see the list of species
            /extract_phenotypes.sh : (avoid it) calls extract_phenotypes and summarize.R
            /get_list_phenotypes.sh :script to see the full list of phenotypes as list_of_phenotypes.txt
            /get_list_species.sh : script to see the full list of species calling extract_species.R
            /see_tsv.sh : visualize the tsvs' in command line
            /see_Metabolic_tsv.sh
            /see_Substrates_tsv.sh
            /summarize.R : calculate balance of every column
            /graphicalsummary.R : calculate balance and make an ASCII plot
        /GenePhene2: A manually compiled set of 41 acidogenic bacterias and hundreds of phenotypes for testing purposes. The 'behind the scenes' for making it is off this folder
            /Metabolic_features.tsv : Database of metabolic activities such as consume X, produce Y
            /Metabolic_traits : Database of metabolic traits such as acidogenics, catalase positive, anaerobic ...
            /get_list_phenotypes.R : script to extract the list of phenotypes
            /get_list_phenotypes.sh : run the R script and keep the outut in a file list_of_phenotypes.txt
            /get_list_species.R : script to extract the full list of species
            /get_list_species.sh : run the R script and keep the output in a file list_of_species.txt
            /summarize.R : script to calculate the balance and coverage (yes/no) of every phenotype
            /graphicalsummary.R : script to get an ASCII plot for balance and coverage
            /see_Metabolic_traits.sh : visulize the database in the command line
        /MDB: Metanogen DataBase : The Phy2Met2 database at http://phymet2.biotech.uni.wroc.pl
            /MDB.csv : the original database
            /MDBcsv2tsv.R : small script to pass it to tab file and compute all the rest of files
            /instructions : some description of how the database was clean.
            /Substrates.tsv : list of phenotypes as input substrates
            /list_of_phenotypes.txt
            /summary_phenotypes_coverage.txt
            /summary_phenotypes_balance.txt
            /graphicalSummary.R : makes an ascii plot
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


