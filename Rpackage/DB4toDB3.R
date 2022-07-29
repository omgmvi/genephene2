#!/usr/bin/env Rscript

# Odin Moron-Garcia
# Developed from
# 23 July 2021 - continuation on 20-30 August
# 30 Sept. 2021 : Update documentation 
# 18 July 2022 - Separation in two or more scripts
#       1. reading JSON onto a 3 tables files 
#       2. fitting model


#### TODO ###



#################
## DESCRIPTION ##
#################
## This script is meant to read 4 files containing : 1. phenotypic data + bacterial names 2. Bacterial names + NCBI taxid 3. NCBI taxid + NCBI species taxid + genome_id 4. Genome_id + gene data
## And perform some SQL-like joins until make three files : Response_data (phenotypes), Feature_data (genotype) and Metadata (other type of data - bacterial names) - They will have a common index column with no column name that will be read as names by R.

# This is a change from previous version that actually make a single table with 1 phenotype, bacterial name and gene data. This 'new' way, I don't have to  handle variable names in the rest of the scripts. For the sake of brevity, I will use the genome_ID as common index between the three tables.

# The configuration file contain the address of the Phenotype DB, the phenotypic trait within the columns of the PhenoDB plus
# taxonomyDB (Name2TaxID), metadata (TaxID2GenomeID)  and genome (GenomeID plus dependant variables for the model ) databases plus some extra variables for file treatment


### How to build datasets are built.

# 1. Read the phenotype DB, extract the bacterial name related info ("Species","Genus","Name" and "Strain") and the phenotypic trait
# 2. Read the TaxonomyDB and match (inner join) microoganism names  (column in phenotypeBD) with its TaxID
# 3. Read the MetadataDB to match the TaxID with the reference genome ID
# 4. Read the GenomeDB that contain genes per column - these coming from COGs, KEGG Orthologs and pFam databases

# The result should be by columns |Genus|Species|Name|Strain|Gene1|...|GeneN

# Merge - inner join sequentially
# PhenoDB -> tax_id -> assembly_accession -> COG|KO|pFam
# Each phase of this join is coded in a function


## S3 classes definitions

# Class Document : list of pairs (folder,file) with the address of a DB in the hard drive

# Class Documents:
    # list of 4 Document object for Phenotype, Taxonomy, Metadata and Genome

# Class Data_Setup: a list of two objects Phenotype_config and Genome_Config


# Class Phenotype_config : list with 4 elements (Phenotypic_trait, Bacterial_metadata_columns, NA_substitution, NA_substitution_value)    
# Class Genome_config : list of 4 elements (Pattern_gene_columns, Pattern_genomic_table_columns, Fix,multiple_genomes, Genome_type)


##### LIBRARIES #####

require(rjson)
require(plyr) # use of dlply and ddply to summarize values
require(uuid)

##### FUNCTIONS ####

## Helper functions ##

read_config_file <- function(folder,file){
    # Given a configuration file and its address in the hard-drive, build the pathname
    # and assuming it is a well-formed JSON file
    # read it all and return it.

    fil<- file.path(folder,file)
    rjson::fromJSON(file = fil)
}


Gene_names <- function(col_names,search_str){
    # The function has a name given by the use in this script
    # it accepts a vector of strings and a regular expression
    # it returns a logical vector with TRUE at the position accepte by the regexp
    # The use:
    # Search the column names from a dataframe for those obbeying a regular expression.
    # This fuction will search a regexp pattern given by the "genome" substructure to help finding the gene names like K00001

	grepl(col_names,pattern =search_str,ignore.case=T,perl = T)
}

# Some "getters" making use of the Gene_names
# They are only wrappers to provide easy access to some columns - and naming the functions to facilitate the script reading

get_genome_columns_function <- function(pattern_gene){
    function(model_data){
       model_data[,Gene_names(names(model_data),pattern_gene)] 
    }
}

get_metadata_columns_function <- function(pattern_gene){
    function(model_data){
        model_data[,!Gene_names(names(model_data),pattern_gene)] 
    }
}

get_phenotype_function <- function(phenotype){
    if(length(phenotype)>1){
        function(model_data){
            model_data[phenotype]
        }

    }else{
        function(model_data){
            model_data[[phenotype]]
        }
    }
}

#Little logging function

log <- function(txt){
    print(paste("[",Sys.time(),":",txt,"]"))
}


##### FUNCTIONS IN TIERS
##### 1. Read configuration files
#####   1.a Checks on the configuration files - 4 types: [Database,Phenotype, (Metadata),Genome] with path and filename + specific option (to be changed for JSON schema validators)
#####   1.b Check on the existence of files, reading the databases - being each one a dataframe (with an extra S3 class)
##### 2. 
##### 
##### 



#####################################################################
# Tier 1a: configuration data loading and manipulation - functions  #
#####################################################################

extract_DB_from_config_files<- function(config_list){
    #This function is not really neccesary - yet ensure the original configuration file conforms with the specification and at least the substructures exists
    # Notice that the name of every list must be folder and file - they are used later!
    # This list facilitate to read in series all the files at once - and upload them in memory at once as well.

    structure(list(
        Phenotype = structure(list(folder = config_list$Phenotype$Folder_Phenotype,   file = config_list$Phenotype$PhenotypeDB),class = "Document"),
        Taxonomy  = structure(list(folder = config_list$Taxonomy$Folder_Taxonomy,     file = config_list$Taxonomy$TaxonomyDB),class = "Document"),
        Metadata  = structure(list(folder = config_list$Metadata$Folder_Metadata,     file = config_list$Metadata$Taxa2assemblyDB),class = "Document"),
        Genome    = structure(list(folder = config_list$Genome$Folder_Genome,         file = config_list$Genome$GenomeDB),class = "Document")
    ),class = "Documents")
}


extract_specific_config <- function(config_list){
# This function serve the purpose of check that the configuration file conforms, at least by names, with the expected nested structure 
# That is, the contents are to checked
    structure(list(
            Database  = structure(list( Database                    = config_list$Database),class = "Database_config"),
            Phenotype = structure(list( Phenotypic_trait            = config_list$Phenotype$Phenotypic_trait,
                                        Bacterial_metadata_columns  = config_list$Phenotype$Bacterial_metadata_columns,
                                        NA_substitution             = config_list$Phenotype$NA_substitution,
                                        NA_substitution_value       = config_list$Phenotype$NA_substitution_value),class = "Phenotype_config"),

            Genome    = structure(list( Pattern_gene_columns            = config_list$Genome$Pattern_gene_columns,
                                        Pattern_genomic_table_columns   = config_list$Genome$Pattern_genomic_table_columns,
                                        Fix_multiple_genomes            = config_list$Genome$Fix_multiple_genomes,
                                        Genome_ortholog_type            = config_list$Genome$Genome,
                                        Genome_type                     = config_list$Genome$Genome_type),class = "Genome_config")
    ),class = "Data_setup")
}


##############################################
# Tier 1b: actual data loading - functions   #
##############################################

# More like helper functions - wrappers to make partial application easily

check_exist <- function(list_of_documents){
        # Check that the folders and least (results of extract_DB_from_config_files function) really exists - just it does it comfy
        # It should actually give an error but I rely on the print functions
    existence<-lapply(list_of_documents,function(ll){
        list(   folder  = file.exists(ll$folder),
                file    = file.exists(file.path(ll$folder,ll$file)))})
    unlist(existence)
}

read_data <-function(folder,file,header = T,sep = "\t",comment.char = "",quote = "",check.names = F,...){
        # A partial application through the function arguments defaults and the ellipsis ...
        # In fact, this function ensure that the files will be tsv with headers and no comment char like # or quotes are expected
        # Yet, it is important to notice that the check.names = F will let spaces separated column names at the phenotypes remaing like such (problematic in some pipelines)
    read.table(file.path(folder,file),header = header,sep = sep,comment.char = comment.char,quote = quote,check.names=check.names ,...)
}

collect_data <- function(list_of_documents){
        # This funcion puts to work the read data function previously by making a "map" on the list of databases - assume that the input is a list with at least elements of name folder and file
        # The check exist ensure this function does not provide and error - but the files content can still be incorrect
    names(list_of_documents)->DB
    data<-lapply(list_of_documents,function(x){read_data(x$folder,x$file)})

}


#####################################################
# Tier 2a: Database processing & joins - functions  #
#####################################################

## Some dedicate treatment of the databases: After reading the databases in tier 1

        ## Phenotype DB:

        # Select the colum(s) of interest (microorganism names and the Phenotypic trait) and get rid of everything else.
        # Choose what to do with the NA at the phenotypic trait - enconded in the NA_substitution <NA2level|NonNA2level|None> and NA_substitution_value parameters of Phenotype structure

process_PhenotypeDB <- function(datos,config){
            # Remove the Phenotype columns we are not using here
    datos$Phenotype <- subset(datos$Phenotype,select=c(config$Phenotype$Bacterial_metadata_columns,config$Phenotype$Phenotypic_trait))
    if(length(config$Phenotype$Phenotypic_trait)==1){
        log("1 var")
        datos <- fix_Phenotype_NA(datos,phenotypic_trait = config$Phenotype$Phenotypic_trait,method = config$Phenotype$NA_substitution,level = config$Phenotype$NA_substitution_value)
    }else{
        log("many var")
        for(Phenotypic_trait in config$Phenotype$Phenotypic_trait){
            #for many phenotypic columns the reassignation of datos may be slow (big data being copied)
            datos <- fix_Phenotype_NA(datos,phenotypic_trait = Phenotypic_trait,method = config$Phenotype$NA_substitution,level = config$Phenotype$NA_substitution_value)
        }
    }
    datos
}


fix_Phenotype_NA <- function(datos,phenotypic_trait,method = "NA2level",level = "Unknown"){
    # The function change the NA for another value (unless set as none) - the replacement value can be chosen with NA_substitution_value
    method <- match.arg(method,choices = c("NA2level","None"),several.ok=F) 
    if(method != "None"){
        if(method == "NA2level"){
            id.na <- which(is.na(datos$Phenotype[,phenotypic_trait]))
            levels(datos$Phenotype[,phenotypic_trait]) <- c(levels(datos$Phenotype[,phenotypic_trait]),level)
            datos$Phenotype[id.na,phenotypic_trait] <- level
        }
    }
    datos
}

    # Genome DB : 
    # The Genome database may have undesired columns that so far can be remove by positive filtering
    # It removes duplicated genomes or average|max|min|...

process_GenomeDB <- function(datos,config){

    # Remove 'genes' (as in KO, COG or Pfam entries) that are not really genes but other info in the database
    datos$Genome<-datos$Genome[,Gene_names(names(datos$Genome),config$Genome$Pattern_genomic_table_columns)]
    # Here is a potentiall source of problems - I should not assume that the GenomeID is gonna be always a "/" - it maybe another symbol
    names(datos$Genome)[names(datos$Genome) == "/"] <- "GenomeID"
    #It shouldn't be any genome twice - and I haven't found it. Yet, I'll check for it.
    if(anyDuplicated(datos$Genome$GenomeID)){ datos$Genome[!duplicated(datos$Genome$GenomeID),]}
    datos
}

fix_multiple_Genomes <-function(datos,config){
    ### Solve the issue of multiple genomes per individual
    ## NOTE: The split is done base on TaxID - so this function needs to be used only after the Big-join
    # helper functions
    pickone <-function(x,na.rm = F){
        #This function just take one of the repeated genomes - Aggregate is an efficient way to apply columnwise and group wise a function
        # but in this case it has an issue - it does not permit for a constant passing through the groups and columns - so how do I pick a random
        #genome? - get the first- (It does allow for passing a parameter, using ...)
        na.omit(x)[1]
    }
    # Identity do not perform well with aggregate and I can't see why

    if(config$Genome$Fix_multiple_genomes == "identity"){return(datos)}

    operation <- config$Genome$Fix_multiple_genomes
    operation <- match.arg(arg = operation,choices = c("median","mean","var","sum","identity","pickone"),several.ok=F)
    operation <- Filter(is.function,ifelse(sapply(c("median","mean","var","sum","identity","pickone"),function(x,y){x == y},y = operation),
                                           c(median,mean,var,sum,identity,pickone),NA))[[1]]
        
    PK_metadata <- unique(datos[,!Gene_names(names(datos),config$Genome$Pattern_gene_columns)])
    PK_genome <- datos[,Gene_names(names(datos),config$Genome$Pattern_gene_columns)]
    row.names(unique(PK_metadata[,c("TaxID",config$Phenotype$Phenotypic_trait)]))->r.n
    aggregate(PK_genome,by=list(TaxID = datos$TaxID),operation)->out
    merge(PK_metadata[row.names(PK_metadata) %in% r.n,],out,by = "TaxID")
    
}


#####################################################
# Tier 2b: Database processing & joins - functions  #
#####################################################

# Now the time comes, peform sequentially inner joins Names -> TaxID -> GenomeID -> Pheno-Geno Table
# In this case all the column names are hardcode, and they will remain so

Big_Join <- function(datos){
        # 1st Associate the phenotype with its TaxID using the name
        datos$Phenotype <- merge(datos$Phenotype,datos$Taxonomy,all=F,by = "Name")
        datos$Taxonomy <- NULL

        # 2nd Associate the phenotypes with the Genome ID using TaxID and Species_TaxID -just in case
        merge(datos$Phenotype,datos$Metadata,by.x = "TaxID",by.y ="taxid",all=F)->tmp1
        merge(datos$Phenotype,datos$Metadata,by.x = "TaxID",by.y ="species_taxid",all=F)-> tmp2

        rbind(tmp1[intersect(names(tmp1),names(tmp2))],tmp2[intersect(names(tmp1),names(tmp2))])->datos$Phenotype
        rm(tmp1,tmp2)
        datos$Metadata <- NULL

        # 3rd Finally associate the phenotypes with the genomes using the the assembly_accession/genomeID
        merge(y = datos$Phenotype,x =datos$Genome,by.y="assembly_accession",by.x ="GenomeID" )->datos$Phenotype
        datos$Genome <- NULL

        # Now, there's only a larger table so the list is not needed any longer
        datos <- datos$Phenotype

        datos
}

Build_triple_DBs <- function(datos,configuration){
    # Split the data in three dataframes with a single Primary Key as row-names -
    # NOTE: I need to check the original data, in the case of Hiding there's duplicated GenomeID (I did not expect that and should not happen)
    # For the sake of carrying on - I make a UUD index

    
   
    IDX <- UUIDgenerate(n = nrow(datos))

    DB <- list()
    DB$Phenotype <- get_phenotype(datos)
    row.names(DB$Phenotype) <- IDX

    DB$Genome <- get_genome_columns(datos)
    row.names(DB$Genome) <- IDX
    
    DB$Metadata <- datos[configuration$Phenotype$Bacterial_metadata_columns]
    row.names(DB$Metadata) <- IDX
    
    DB
}

Save_triple_DB <- function(datos,folder_output){
    # It is expected - that datos have 3 dataframes - I won't take control of this here
    write.table(x = datos$Phenotype, file = file.path(folder_output,"Phenotype.csv"),sep = ";")
    write.table(x = datos$Genome, file = file.path(folder_output,"Genome.csv"),sep = ";")
    write.table(x = datos$Metadata, file = file.path(folder_output,"Metadata.csv"),sep = ";")
}

###################################################################


######## VARIABLES ########################

args = commandArgs(trailingOnly=TRUE)

if(length(args)  != 3){
    print("Not enough arguments provided, I carry on with presets")

    folder_config_file          <-  "./"
    config_file                 <-  "Hiding_severalpheno_COG_BoW_identity_genomes.dat"
    folder_output               <-  "./"
}else{ 
    print(args);
    folder_config_file <- args[1]
    config_file <- args[2]
    folder_output <-  args[3]
}

####################################################################
##################################### SCRIPT #######################
####################################################################

log("Starting the Script")
log("READING CONFIGURATION FILES ...")
#### Tier 1a : Get configuration files  #####

read_config_file(folder_config_file,config_file)->setup
Documents <- extract_DB_from_config_files(setup)
configuration <-extract_specific_config(setup)

log("CONFIG FILES READ")
log("Preparing some helper functions")

get_genome_columns <- get_genome_columns_function(configuration$Genome$Pattern_gene_columns)
get_metadata_columns <- get_metadata_columns_function(configuration$Genome$Pattern_gene_columns)
get_phenotype <- get_phenotype_function(configuration$Phenotype$Phenotypic_trait)

#### Tier 1b : Read all data

log("QUICK CHECKS")
##Check all exists; These checks are now more like class constructors!

print(check_exist(Documents))
stopifnot(all(check_exist(Documents)))

log("CHECKS PASSED")

log("COLLECTING THE DATA")
## Collect the data
collect_data(Documents)->datos

log("PHENOTYPE, TAXONOMY, METADATA AND GENOME DATABASES READ")


##### Tier 2a:  Process databases ####
log("START OF DATA BASES PROCESSING")

#bypass to make tests!
#configuration$Phenotype$NA_substitution <- "NA2level";configuration$Phenotype$NA_substitution_value <- "U"
datos <- process_PhenotypeDB(datos,configuration)

log("PHENOTYPE DB PROCESSED")

datos <- process_GenomeDB(datos,configuration)

log("GENOME DB PROCESSED")

##### Tier 2b:  create the working database ####
log("GENERATING THE PHENO-GENO TYPE DATABAES")
datos <- Big_Join(datos)

log("DATABASE GENERATED AFTER JOIN")

#Now that I have joined the Taxonomy Indentifiers and the Genomes, I can trim down the number of genomes per bacteria

log("FIXING MULTIPLE GENOMES")

datos <-fix_multiple_Genomes(datos,configuration)

# NOTE: From here onwards, the variable datos does not contain a list but a data frame just with the bacterial name metadata, the column with the phenotype and all columns with genes
log("MULTIPLE GENOMES FIXED")

log(" SAVE THE OUTPUT")    
datos <- Build_triple_DBs(datos,configuration)
Save_triple_DB (datos,folder_output)

log("OUTPUT SAVED")
log("THE END")
