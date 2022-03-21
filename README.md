# **Genephene2** : *A set of models to predict microbe metabolic phenotypes from their annotated genomes*

**TL;DR**

>At the moment this is an script to fit machine learning models for classification of phenotypic traits according to presence/abscense of genes but flexibly build to accept more types of genomic information. At the moment, the script fit the models and return files with their performance in classification. 

>For genomic information based on gene frequency, the accuracy range around 0.80-0.90 (roughly) but for genomes transformed in continuous embedding (Doc2Vec) the results drop dramatically

**TODO:** 
- [ ] Find new ways to do genomic embeddings.
    - [ ] Try Kohonen self-organizative networks and 
    - [ ] Try Dictionary learning
    - [ ] Try Scatter Search and GA based logistic regression 
- [ ] Describe the reason for low performance in the current models
- [x] Check if Doc2Vec is failing due to coding mistakes
- [ ] Add more machine learning models
- [ ] Find an OOP solution to reduce coding and flexibilize user input of functions
- [ ] Move the prototype to python as separated scripts
- [ ]  Use a CMake kind of script to unroll the databases - as it is now sounds super complicated and easy to forget

**QUICK LINKS**
- A description of how data has been and need to be input
- The Doc2Vec genome embedding
- How to run a model
- How to elaborate a larger number of models 'as in a HPC' or SnakeMake

## Introduction
GenePhene2 is a continuation of the GenePhene project from Prof. Christopher Quince and his collaborators.

The aim of this project is to find appropiate  machine learning model able to predict bacteria and archaea metabolic phenotypes from their genomic functional annotation having having in mind metagenomics experiments. More specifically, it is now possible to take environmental DNA samples from many environments, e.g. soil, water, sludge, air, that can now be sequenced either in long reads, i.e. Oxford Nanopore, or short reads, e.g. Illumina. Those samples contain many species of  bacterias and archaeas that cannot be cultivated using standard laboratory technique and therefore little knowledge of their phenotypic characteristics are known. On the other, hand modern sequencing technologies and associated bioinformatic tools has facilitated the recovery and reconstruction of 'metagenomic assembled genomes",~~i.e. single species genomes reconstructed by assembly contiguous sequences and further filtering of quimeric sequences~~. Therefore, we have access now to new species or strains whosemetabolic activities are only deduced by comparison with their siblings taxa, calculated by RNA16S or similar barcoding.~~The approach of phenotype imputation by phylogenetic association is troublesome in the sense that many strains or new microbes species are actually defined by their differential phenotypic capabilities from their closest siblings. In short, new taxa are defined as those that do something their sister taxa don't.~~ The phylogenetic based phenotypic imputation is the chosen technique used by some popular database-algorithm combination, e.g. FAPROTAX, Tax4Fun, MIDAS, that maps RNA 16S profiles to their imputed metabolic activities.
The approach taken here is linking the functionally annotated genomes, after Open Reading Frame and Ortholog gene imputation, to more or less specific phenotypes such like 'aerobe' or 'catalase activity'. Each phenotype is associated after fitting supervised classification algorithms where  phenotypes play the role of predicted variables and ortholog genes the one of predictors.

To solve several difficulties regarding undetermined matrix (many more genes than individuals and therefore columns are not independent), phylogenetic signal in the gene space and other correlations, we test to transform annotated genomes in continuous numerical variables using the NLP (for Natural Language Processing) algorithm  Doc2Vec instead of the so called Bag of Words approximation. If possible, other _feature engineering_ methods would be tried.

At this stage, several phenotypic databases will be used for training the models and check the performance of them, being them GenePhene ( an small custom database of many metabolic traits from acidogenic taxa for testing purposes), FAPROTAX, Hiding in Plain Sight and MDB (Methanogens Database). The database MIDAS was dismissed due to the low number of species, being most entries genera or families. These databases have been curated to keep only entries corresponding to bacterial with full names and their metabolic related phenotypes. We dismiss taxa levels over species, i.e. genera, families and above, and numerical phenotypes like temperature range or non-metabolic ones like sporulation. Therefore, the traits kept were of type _Consume resource X_, _Produce resource Y_ or more general trtraits like _Catalase activity_. Yet, each database encode a quite different set of traits and a quite non-overlapping number of bacterias.

# CODE ORGANIZATION

The code, at the moment, revolt around a single script, model\_setup.R that perform all the next steps in the analysis. The script and associate 'sender' scripts are found at ~/Genephene2/Models/scripts/. All this works as a homemade job queue where the model\_setup.R is called using a JSON file as input parameter(plus some others). The sender is tipically an script to apply the model\_setup.R to all the files in a folder or similar. In addition, an script to generate the JSONs has been prepared, it uses another JSON file with all the information required to generate the full set of final JSONfiles.

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

# Database organization

As explained above, the database are organize around four tables. Albeit being a case of 'technological debt' the tables are not encoded in a relational (SQ)  database manager which would help and accelerate the whole system. Thus, they text files, tab separated (tsv) and with \n as EOL and UTF-8 encoding, were decide to be saved in separated folders per type, e.g. Genome, Metadata, ... and per original phenotype database. Names and folder can be easily changed and rearranged since the information is passed to the programs through JSON files (see below). Yet a short description follows with the files and scripts elaborated for each type.

## Description of folders / databases
This is a short description of the folder and file structure used in this project. The most relevant files are explained shortly
The databases are organized so each type of table is in a different folder, and within those a folder for every dataset, e.g. FAPROTAX, Hiding ...
Small scripts to perform small tasks, sort of untar files, add headers, etc and that was worthless to make it complex to check assumptions or add parameters are left in the same folder where they are used.
Scripts prompt to be reused are usually moved to a base folder named scripts and accept parameters.

- ../GenePhene2/
   -  /Databases : _The name would be better **Phenotypes** since it contains each metabolic phenotype databases._
                 Each subfolder refer to the phenotypic database, corresponding to a published one.
        - /FAPROTAX
            - /FAPROTAX_1.2.4_zip : The original FAPROTAX database. contain the files FAPROTAX.txt with the whole phenotypes in a unstructure format
            - /FAPROTAX.tsv : The phenotypic database
            - /see.tsv.sh : script to visualize the database in the bash command line
            - /get_list_of_phenotype.sh : script to extract the phenotypes
            - /get_list_of_speces.sh : script to extract the names of the species
            - /parse_faprotax.R : script to parse the original files from FAPROTAX
            - /summarize.R : script to calculate the balance of phenotype columns
            - /graphicalsummary.R : script to make an ASCII plot for the balance of phenotypic columns
        - /Hiding : database from the paper Hiding in plain sight - one of the databases backing BacDive
            - /instructions : description of steps to convert from the original database IJSEM_pheno_db_v1.0.txt to ijsem.txt
            - /Substrates.tsv : Database of Substrates (Consume X)
            - /Metabolic.tsv : Database of several metabolic activities
            - /extract_phenotypes.R : script to pass from a manually corrected DB (see instructions) wiht UNICHAR and \n EOL to ijsem.tsv with ranges split as columns and others
            - /get_phenotypes.R: Split ijsem.tsv onto two tables one for Substrates and the other for Metabolic activities.
            - /extract_species.R : scripts to see the list of species
            - /extract_phenotypes.sh : (avoid it) calls extract_phenotypes and summarize.R
            - /get_list_phenotypes.sh :script to see the full list of phenotypes as list_of_phenotypes.txt
            - /get_list_species.sh : script to see the full list of species calling extract_species.R
            - /see_tsv.sh : visualize the tsvs' in command line
            - /see_Metabolic_tsv.sh
            - /see_Substrates_tsv.sh
            - /summarize.R : calculate balance of every column
            - /graphicalsummary.R : calculate balance and make an ASCII plot
        - /GenePhene2: A manually compiled set of 41 acidogenic bacterias and hundreds of phenotypes for testing purposes. The 'behind the scenes' for making it is off this folder
            - /Metabolic_features.tsv : Database of metabolic activities such as consume X, produce Y
            - /Metabolic_traits : Database of metabolic traits such as acidogenics, catalase positive, anaerobic ...
            - /get_list_phenotypes.R : script to extract the list of phenotypes
            - /get_list_phenotypes.sh : run the R script and keep the outut in a file list_of_phenotypes.txt
            - /get_list_species.R : script to extract the full list of species
            - /get_list_species.sh : run the R script and keep the output in a file list_of_species.txt
            - /summarize.R : script to calculate the balance and coverage (yes/no) of every phenotype
            - /graphicalsummary.R : script to get an ASCII plot for balance and coverage
            - /see_Metabolic_traits.sh : visulize the database in the command line
        - /MDB: Metanogen DataBase : The Phy2Met2 database at http://phymet2.biotech.uni.wroc.pl
            - /MDB.csv : the original database
            - /MDBcsv2tsv.R : small script to pass it to tab file and compute all the rest of files
            - /instructions : some description of how the database was clean.
            - /Substrates.tsv : list of phenotypes as input substrates
            - /list_of_phenotypes.txt
            - /summary_phenotypes_coverage.txt
            - /summary_phenotypes_balance.txt
            - /graphicalSummary.R : makes an ascii plot
    - */Genomes* : 
        - /FAPROTAX
            -   */BagOfWords_Genome* : The folder containing the gene frequencies genome
                - **/FAPROTAX_BoW.tar.gz** : The compressed file with gene frequencies for COG, KEGG and pFAM
                - /extract_db.sh : Companion script that uncompress the tar file and correct the name of the column Genome_ID <--> GenomeID
            -   */Doc2Vec_Genome* : The folder containing the required file for generating the Doc2Vec embeddings - note they are not the embedding themselves
                -   **/FAPROTAX_D2V.tar.gz** : The files within contain the list of genes annotated as a list of contigs as genome_X geneA geneB ... ordered as found in the genome
                                                These files have to be cleaned for too small contigs and also pass through the script that generate the embedding - plus adding the headers
        - /Hiding
            - */BagOfWords_Genome*: Similar to the FAPROTAX one
            - */Doc2Vec_Genome* :  Like the one in the FAPROTAX folder but in this case the file had to be split to small sized (due to github) and have to be reconstructed
                - HIDING_D2V.tar.gz.part_* : These files ending in part_aa part_ab and part_ac are recomposed by cat HIDING_D2V.tar.gz.part_a* > HIDING_D2V.tar.gz 
                                                then the same process as described in the FAPROTAX 
        - /GenePhene2: Same as FAPROTAX
        - /MDB: Same as FAPROTAX
        - **/scripts**: several helper scripts - including the one to build the Doc2Vec embeddings
            -   Add_D2V_header.sh : Just add the header to the final D2V file (so the last script to run) - work in a list of files when using globs (within "")
            -   remove_contig.sh : take a filename and a contig lenght and filter out those rows that have less than the length
            -   remove_contig_files.sh : Apply remove_contig.sh to a list of files using a wildcard (within " ") - needed to solve the glob expansion
            -   /D2V_scripts : The script in python to build the D2V embeddings - originally several, they have been grouped into a single one that need to be modified to select the database
                    - **/Doc2Vec_Genomes.py**: This is the script that convert contigs in embeddings - It requires to be told the number of dimensions desired plus the folders and files
                                                At the moment, the parameters are hardcoded and need to be changed to accept parameters as inpput
    - /Taxonomy : it contains the taxonomic metadata, it really only is a translation between microorganism name in the phenotype database and the Taxon ID in NCBI to find its Reference Genome
        - /FAPROTAX
            - /SpeciesTaxID.tsv : The database for taxonomy
            - /SpeciesName2TaxID.R : a helper script to search for the taxid in the NCBI taxonomy database
            - /known_mistakes_or_variants.txt : list of wrong names and their corrections
        - /Hiding:
            - /SpeciesTaxID.tsv : The database for taxonomy
            - /SpeciesName2TaxID.R : a helper script to search for the taxid in the NCBI taxonomy database
            - /known_mistakes_or_variants.txt : list of wrong names and their corrections
            - /missingSp.txt : another list of wrong names and corrections
        - /GenePhene2
            - /SpeciesTaxID.tsv : The database for taxonomy
            - /SpeciesName2TaxID.R : a helper script to search for the taxid in the NCBI taxonomy database
        - /MDB
            - /SpeciesTaxID.tsv : The database for taxonomy
            - /SpeciesName2TaxID.R : a helper script to search for the taxid in the NCBI taxonomy database
        - /instructions : How to download and untar the NCBI database
    - /Metadata : It contain  genomic information metadata, at the base level there are the TaxID2Genome.txt files that is metadata to connect the microorg. name to its genome
        - /FAPROTAX : 
            - /TaxID2Genome_metadata.tsv : the table to connect the assembly accession to the taxid and species taxid
            - /genome_metadata.tsv : Information to let Seb Raguidau to download and annotate the genomes
        - /Hiding
            - /TaxID2Genome_metadata.tsv : the table to connect the assembly accession to the taxid and species taxid
            - /genome_metadata.tsv : Information to let Seb Raguidau to download and annotate the genomes
        - /GenePhene2
            - /TaxID2Genome_metadata.tsv : the table to connect the assembly accession to the taxid and species taxid
            - /genome_metadata.tsv : Information to let Seb Raguidau to download and annotate the genomes
        - /MDB
            - /TaxID2Genome_metadata.tsv : the table to connect the assembly accession to the taxid and species taxid
            - /genome_metadata.tsv : Information to let Seb Raguidau to download and annotate the genomes
        - **/scripts**_: Here is the method to get the metadata from the taxonomy files and the refseq (from NCBI)
            -   /TaxID2Refseq.R : The actual script to generate the metadata from the assembly_summary_refseq.txt and the corresponding SpeciesTaxID.tsv from Taxonomy folders
    - **/Models** : It contain the information to make the models particularly the scripts to generate the JSON config files
        -   /scripts
            -   **/model_setup.R** : THE script/program that run the model [look at the how to section]
            -   /expand_json_model_configuration.R : a companion script that helps to create a full list of JSON files each corresponding a run fo the model_setup.R
            -   /templates
                -   /expansion_GenePhene2_test.sh : an example of how to expand a configuration file using expand_jso_model (the script is useful when needed to repeat during the preparation)
                -   /send_all_folder_model.sh : An example of how to send all jobs to the bash as if it where a job queue in a HPC (unfortunately they arent)
        - **/config.files** :  Here are the configuration files used for generate the input data
            -   /FAPROTAX : A config.file that includes all phenotypes from FAPROTAX database and several types of BagOfWords genomes to be expanded
            -   /**GenePhene**  : A config.file with all phenotypes (in two separate tables) and several genomes of type  BagOfWords and Doc2Vec
            -   / ... : several other examples for Hiding and MDB and test with less phenotypes
            -   **/model.glmnet_elasticnet.json** : a configuration file with the options needed to run a logistic regression with elastic net regularization
            _   /model.Naive_Bayes.json : a config. file to run a Naive Bayes
            -   /model.XGBoost.json
            -   **/models.json** : a multimodel example, here the example runs together (in one go of model_setup.R) the glmnet, Naive Bayes and XGBoost models in one phenotypic dataset.

## BUILD DATABASES
### PHENOTYPES
    The phenotype databases are not compressed and can be used as is. The scripts assciated to them were useful to reduce the lenght of commands during the developing phase.
### GENOMES
#### Bag of Words genomes.
These genomes are just tables with a GenomeID column that contain the 'name' of the genome to be used as primary key to execute the join with the other tables. The rest of columns are the frequency tables for each gene of every ortholog type. The ortholog genes has not been completed between dataset (FAPROTAX, HIDING,....) so, each one can contain a different number. Also, different contigs are in different rows so the GenomeIDs are not unique but Doc2Vec takes care of that.

To build the BoW genomes just untar the compressed file and use awk to change the name of the column Genome_ID to GenomeID.
> tar -zxvf dataset_BoW_ortholog-type.tar.gz
> awk 'BEGIN {OFS = "\t"} NR == 1 {if ($1=="Genome_ID"){$1 = "GenomeID"}} 1' FAPROTAX_BoW_COG.tsv
to simplify the process there's a joint script doing that for the three ortholog types of genomes
> extract_db.sh

#### Doc2Vec
The files named _dataset_D2V_Ortholog-type.tsv are a row-based list of contigs whose first element is the GenomeID (notice there's no header since this is not a table but the lenght of the rows are very different) and the rest of the elements are the ortholog genes annotated in a contig, kept in order. These lines can be used for Doc2Vec as follows
- Gensim Python library can execute a Doc2Vec neural network for embeding provided a "sentence" and a "tag". In this case the sentence would be the contig and the tag the genome.
- Gensim D2V is not able to update the model in batches (that does not mean that it does not do batch gradient descent - I have no clue if it does) so all the genome has to be given in one go, although the library have a tool to map the files to the memory so big data are not a problem for fitting the model.
- Gensim D2V needs a the number of dimensions of the embedding (I've been using 50, 100, 500, 1000, 5000) and the tag - I have modified the code so the tag is the first element of the row rather the name of the file like in the original one - also if the internal model is fitting as a skip-gram network or a Continuous BagOfWords (sorry for the misunderstanding but this is not the same as in our BagOfWords genomes, it is reusing the terminology for related but different behaviour of the internal model) that translate to two types of possible Doc2Vec, the Distributed Memory and the Distributed Bag of Words .
- The number of epoch is kept as default in 40 and the default window size seems to be 5.


The steps then are
1. untar the tar files to get the 3 orthologs genomes positionally ordered as described before.
2. Filter out small contigs using, for example for FAPROTAX
> ~/GenePhene2/Genomes/scripts/remove_contig_files.sh "../FAPROTAX/Doc2Vec_Genome/FAPROTAX_SINGLE_D2V\_\*.tsv" 10
3. run Doc2Vec_Genomes.py after modifying the hard-code parameters (just Ndim = 100, Ortho = "KEGG", Geno = "FAPROTAX") - This must change to accept the parameters in the CLI and should take a minute
>python3  /GenePhene2/Genomes/scripts/D2V_scripts/Doc2Vec_Genomes.py
4. Add headers to the new files in the folder
>bash  ~/GenePhene2/Genomes/scripts/Add_D2V_header.sh "../FAPROTAX/Doc2Vec_Genome/Faprotax_Doc2Vec_*_KEGG.tsv"

### THE JSON FILE
### BIG JOIN
# LEFTOVERS (TO REMOVE)
## Pipeline
Shortly, the pipeline has to start with the phenotype database and 'join' the microorganism name with the name2taxid and this one with the taxid2genome and lastly with the KEGG_summary.tsv or cogs_summary.tsv or pfamA_summary.tsv. Graphically,

Databases/FAPROTAX/FAPROTAX_1.2.4/FAPROTAX.tsv <--> Taxonomy/FAPROTAX/SpeciesTaxID.tsv <--> Genomes/FAPROTAX/TaxID2Genome_metadata.tsv <--> (GENOME FILE STILL NOT CREATED)
Databases/GenePhene2/Metabolic_features.tsv <--> Taxonomy/GenePhene2/SpeciesTaxID.tsv <--> Genomes/GenePhene2/TaxID2Genome_metadata.tsv <--> Databases/GenePhene2/Odin_data/cogs_summary.tsv (Still not Structured)


