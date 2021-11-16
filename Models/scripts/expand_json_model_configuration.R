#!/usr/bin/env Rscript
### Odin Moron-Garcia 

## Script to expand a general JSON file containing the folders, file address and some parameter of phenotypic traits and genomic database
# The reason is to write a configuration file that contain the phenotypic database and the desired phenotypes for further modelling
# After the expansion multiple files with a single phenotype will be ready for passing to the modelling pipeline - so it is ready for HPC array jobs


## The files are meant to be in the folder ~/Models.config.files and having an structure like :
#   {   "Database" : <DATABASE_NAME>,
#       "Phenotype" : 
#           {
#               "Folder_Phenotype" : <Folder where phenotype can be found>  # e.g. "/home/ubuntu/Databases/GenePhene2/",
#               "Phenotypes" : [                                            # note: The phenotypes can be in more than one file - particularly when they are separated in Substrates or Activities
#                                                                           #       Therefore, different files and subsets of phenotypic columns can be chosens (is a list  of tuples [{},{},...,{}])
#                       { "PhenotypeDB" : <PHENOTYPIC DATABASE FILENAME>,"Phenotypes"  : [<Pheno1>, <Pheno2>,..., <PhenoN>]},
#                       { "PhenotypeDB" : <PHENOTYPIC DATABASE FILENAME>,"Phenotypes"  : [<Pheno1>, <Pheno2>,..., <PhenoM>]},...
#                   ],
#       "Bacterial_metadata_columns" : [<col1>,<col2>,...,<colN>],  ##At the moment: ["Name","Genus","Species","Strain"]
#       "NA_substitution" : <NA2level>,
#       "NA_substitution_value" : <Unknown|NA|0|...>
#           }
#
#       "Taxonomy" :
#           {
#               "Folder_Taxonomy" : <Foldere where the Taxonomy files can be found>,    # e.g. "/Home./ubuntu/Genomes/Taxonomy/GenePhene2/"
#               "TaxonomyDB" : <Species Name 2 Tax ID files>                            # Note: At the moment all these files are called "SpeciesTaxID.tsv" with 2 columns Name and TaxID
#           }
#       "Metadata" :
#           {
#               "Folder_Metadata" : <Folder with Metadata file>,                        # e.g. "/home/ubuntu/Genomes/GenePhene2" since at the moment the file is kept under the genome folder
#               "Taxa2assemblyDB" : <Taxi ID 2 Genome file>                             # Note: At the moment all these filas are called "TaxID2Genome_metadata.tsv"
#           }                                                                           #       with columns assembly_accession      taxid   species_taxid   organism_name
#
#       "Genomes": [                                                                    # Note: Our Genomes come in 3 types of orthologs and will be soon as Doc2Vec embbedings. 
#                                                                                       #       We put here all the available genomes in a list of tuples
#               {   "Genome" : <Ortholog type>,                                         # e.g. "KEGG"|"COG"|"pFam"|"D2V_KEGG"|"D2V_COG"|"D2V_pFam"
#                   "Genome_type" : <a string for the type of  genome>                  # e.g. "BoW"|"D2V"|"D2V50"|"D2V100"
#                   "file_ending" : <an open field to make distinction between files>
#                   "Folder_Genome" : <Folder to locate the Genomic DB>,                # e.g. "/home/ubuntu/Genomes/GenePhene2/odin_data",
#                   "GenomeDB" : <File containing the genome>,                          # e.g. "KEGG_summary.tsv"|"COGS_summary.tsv"|"pFam_summary.tsv"
#                   "Pattern_gene_columns" : <Perl like regexp>,                        # Note: How to find the name of the gene columns "K\\d+"|"COG\\d+"|"PF\\d+\\.\\d+"
#                   "Pattern_genomic_table_columns" : <Perl like regexp                 # Note: How to find all non-metadata columms "K\\d+|\/"|"COG\\d+|\/"|"PF\\d+\\.\\d+|\/"
#               },
#               { "Genome" : <>, "Folder_Genome" :<>,"GenomeDB" : <>,"Pattern_gene_columns" : <>, "Pattern_genomic_table_column":<>},...
#           ]                                                                   {
#   }

# This script 'unpack' the several phenotypes into a single phenotype DB and phenotype column and the several Genomes into a single genome generating many JSON formats for a single phenotype and genom

## NOTE: The script has the file to unpack harcoded in the global sections - TO DO: accept it as parameter
#############################
#                           #
#           PREAMBLE        #
#                           #
#############################

##### EXTERNAL LIBRARIES #####

require(rjson)


##### FUNCTION LIBRARY  ######


### Unpacking ###


## Phenotypes list

## These functions are responsible of unpacking phenotype databases and their columns.

# The general idea behing these functions  is to use an lapply (an for loop likeish iteration over elements)
# So, to each element of list-like structure, copy the original structure as having a single element and get rid of the list structure
# e.g.  Pheno : [ p1, p2, p3 , ... , pn] -> [Pheno : p1 , Pheno: p2, Pheno : p3, ..., Pheno : pn]


substitute_PhenotypeDB <- function(Phenotypes_list,Full_list){
        # If There are more than one Phenotype Database - i.e. more than one file with phenotipic traits columns
        # Separate the structure in a list of structures each one with a single phenotypic file - We want eventually each analysis done with a single phenotype
        # Remove the multiple PhenoDB list to make a list of single PhenoDB (which still contains a list of many phenotipic traits)
        # NOTE : The output structure is now different in the names since the list Phenotypes pass to be called Phenotype (to represent is not a list anymore)
        # and the sublist Phenotypes pass to be called Phenotype_traits - This may not be the most ideal situation but they are temporal and internal structures only.
        # NOTE2 : Notice that the method use the fact that R pass function arguments by values (copy on modify really that acts as pass by value).
    Full_list$Phenotype$PhenotypeDB <- Phenotypes_list$PhenotypeDB
    Full_list$Phenotype$Phenotype_traits <- Phenotypes_list$Phenotypes
    Full_list$Phenotype$Phenotypes <- NULL
    Full_list
}


substitute_Phenotype_traits <- function(Full_list){
    ## This unpack the list of phenotypic traits that is implemented as columns of a table, to be a list of structures with a single phenotype
    ## TO DO - does it work fine when there are more than one phenoDB in the structure? I don't think so, if not - some error or warning should be given

    # The expected input -> output is Phenotypes : {PhenoDB : <Filename>, Phenotypes : [pheno1,pheno2, ..., phenoN] } -> The  list of copies from the original list
    # for Database, metadata, taxonomy, genomes and Phenotype but 
    # with Phenotypes : {PhenoDB : <Filename>, Phenotypes : phenoK}
    # NOTE: The new Phenotypes structure change again names  to be Phenotypic_trait rather than Phenotypic_traits reflecting that there's now a single pheno trait in the whole structure
    lapply(Full_list$Phenotype$Phenotype_traits,function(Ptl){
        Full_list$Phenotype$Phenotypic_trait <- Ptl
        Full_list$Phenotype$Phenotype_traits <- NULL
        Full_list
    })
}

## Genome list
# The original document can contain several types of genomic files - according to the type of ORF functional classification COG | KEGG | pFam or either using word embbeddings or not
# In this case the unpacking is similar to the Phenotype object unpacking, the Genomes is a list of single genomes that gets unpacked onto a list of copies of the complete document
# now with a single genome . NOTE: To reflect the change from many genome files to a single one the structure is renamed from Genomes to Genome

substitute_Genome <- function(Full_list){

    lapply(Full_list$Genomes,function(Gen){
        Full_list$Genome <- Gen
        Full_list$Genomes <- NULL
        Full_list
    })
}

### Other helper functions ###

# These just assist the construction of filenames and homogeneous writing of files (the export of the whole unpacked list documents)

create_output_filenames <- function(x){
    # Assumming x as a instance of a experiment document containing the elements Database, Phenotypic trait and Genome type build an underscored separate name with .dat extension
    gsub(x = paste(paste(x$Database,x$Phenotype$Phenotypic_trait,x$Genome$Genome,x$Genome$Genome_type,x$Genome$file_ending,sep = "_"),".dat",sep = ""),pattern="\\s",replacement = "_")
}

## Small function that handle the export of single experiment documents to JSON and then to a file - it contain an print for confirming the export but lacks any test or error handling

write_files <- function(Full_list,file,folder){
    f_handler <- file(description=file.path(folder,file),open = 'w')
    writeLines(toJSON(Full_list),f_handler)
    close(f_handler)
    print(paste0("File saved:",file.path(folder,file)))
}


#########
# SETUP # 
#########

# Variables that should be inputed
args <- commandArgs(trailingOnly=T)

if(length(args) != 3){
    ## GENEPHENE2
    # Where the original multi-phenotype multi-genome list is?
    folder_input <- "../config.files/"
    file_input <-"GenePhene2_test"
    #file_input <-"GenePhene2"
    #Where do you want to put the full list of single-phenotypeDB, single-phenotypic trait,  single-genome experiment document?
    folder_output <- "/home/ubuntu/Models/GenePhene2/test.files/"
    #folder_output <- "../GenePhene2/data.files/"


    #FAPROTAX
    #folder_input <- "../config.files/"
    #file_input <-"FAPROTAX_test"
    #folder_output <- "/home/ubuntu/Models/FAPROTAX/test.files/"
}else{
    folder_input  <- args[1]
    file_input    <- args[2]
    folder_output <- args[3]

}
##########
# SCRIPT #
##########

### The script is a single call to the many functions written above

# READ THE multi-phenotypeDB, multi-phenotypic traits, multi-genome document assuming it is JSON and assuming it conforms with the structure
## TODO: Test the structure for conformance

fromJSON(file = file.path(folder_input,file_input))->JSON_list

### Start unpacking the multi-phenoDB into a list of documents single PhenoDB

lapply(JSON_list$Phenotype$Phenotypes,substitute_PhenotypeDB,Full_list = JSON_list)->JSON_list

# Unpack each PhenoDB with a list of Phenotypic traits into a list of documents with a single trait
lapply(JSON_list,substitute_Phenotype_traits)->JSON_list
JSON_list <- unlist(JSON_list,recursive = F) #lapply-ed functions need some cleaning afterward - in this case there's a list of lists that needs to be one level up.

## Unpack the multi genome list into a list of single genome documents.
lapply(JSON_list,substitute_Genome)->JSON_list
JSON_list <- unlist(JSON_list,recursive = F)


### EXPORT THE LIST

# Bulk-generate the filenames for every element of the list

lapply(JSON_list,create_output_filenames)->output_filenames
## Send the list of documents to a writer - notice that with the multiple apply (mapply) we send a document and a file name to an instance of the function.
mapply(FUN = write_files,JSON_list,output_filenames,folder = folder_output)
