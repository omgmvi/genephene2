#!/usr/bin/env Rscript

# Odin Moron-Garcia 23 July 2021 - continuation on 20-30 August
# 30 Sept. 2021 : Update documentation 

# This script reads a JSON document as a configuration file  with the information for building the (see the script expand_json_model_configuration.R for a description)
# Currently the models are hard-coded in this script but soon to be an option - always thinking on paralelization

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

## NOTES:
# Currently, there is the possibility that a microorganisms have more than one genome rows, wich expand the phenotypes - there's here function to either leave as is or aggregate them or pick the first

# At the moment I uses S3 classes to take advantage of the type system in R, that is, I am using class(object) <- "type" to tag wether a variables is of certain type or other (for example if it is a ssingle model or multiple models)

#TODO: 
    ## Make a model document
    ## Test the DB has the expected structure
    ## the fix_multiple_genomes needs inspection - also needs an argument for the column to use and/or check that it exists
    ## Define the output
    ## Define what to do with the errors.

###################### SETUP ###################

## S3 classes definitions
# Class Documents:
    # list of 4 Document for Phenotype, Taxonomy, Metadata and Genome
# Class Document : list of pairs (folder,file) with the address of a DB in the hard drive

# Class Data_Setup: a list of two objects Phenotype_config and Genome_Config
# Class Phenotype_config : list with 4 elements (Phenotypic_trait, Bacterial_metadata_columns, NA_substitution, NA_substitution_value)    
# Class Genome_config : list of 4 elements (Pattern_gene_columns, Pattern_genomic_table_columns, Fix,multiple_genomes, Genome_type)


# Class model_setup
    #subclases (or 2nd classes) glmnet, Naive_Bayes and XGBoost.
    # Both types are a list containing package,model, train_CV_method <CV|LOOCV>, train_test_split <positive numeric>
## Class multi_model_setup
    # A list of many model_setup objects
# Class multi_train a list of several train objects (train class is from caret package see methods(class = "train")
##### LIBRARIES #####
require(rjson)
#require(RJSONIO)
require(plyr) # use of dlply and ddply to summarize values
require(caret) # For models, in particular partition of data

find.package("pROC",quiet =F)
find.package("ROSE",quiet =F)
find.package("smotefamily",quiet =F)
##### FUNCTIONS ####

    ## Helper functions ##

# Assuming each configuration file is in proper format (JSON and with the corresponding structure)
read_config_file <- function(folder,file){
    fil<- file.path(folder,file)
    rjson::fromJSON(file = fil)
}


# This fuction will search a regexp pattern given by the "genome" substructure to help finding the gene names like K00001

## Search the names of a dataframe for those obbeying a regular expression.
Gene_names <- function(col_names,search_str){
	grepl(col_names,pattern =search_str,ignore.case=T,perl = T)
}

        ## Data input and manipuation functions - They require a configuration file to have been read

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

check_model_config <- function(config){
  # An individual check of every type of model - here class would work perfectly - to assert that all the models for options are available and with correct values - 
    #Check that the root of the list is called "model" and it contains an element called "model" (as well) with values glmnet, nb and xgbTree - then I can do the checkings according to the type of model
    #helper function
    checker <- function(object,existance,available_options){
        if(!exists(x=existance,where=object)){stop(simpleError("existance"))}
        
        if(!is.null(available_options)){
            if(!object[[existance]] %in%  available_options){stop(simpleError("content"))}
        }
    }

    check_model <- function(config){
        
        # Top  'node' is "model"
        checker(config,"model",NULL)

        # It exist a node 'model' - I do not check if valid one - but help the user telling the options
        checker(config$model,"model",c("glmnet","Naive_Bayes","XGBoost"))

        if(config$model$model %in% c("glmnet","Naive_Bayes","XGBoost")){
            checker(config$model,'package',c("caret"))
            checker(config$model,'PCA',c("None"))
            checker(config$model,'train_CV_method',c("CV","LOOCV"))
            checker(config$model,"train_test_split",NULL)
            tryCatch(stopifnot(is.numeric(config$model$train_test_split)),error = function(e){stop(simpleError(paste("parameter train_test_split from  model configuration is not a number but",config$model$train_test_split)))})
            tryCatch(stopifnot(config$model$train_test_split>=0 && config$model$train_test_split<=1),error = function(e){stop(simpleError(paste("Value must be between 0 and 1 and is",config$model$train_test_split)))})

            checker(config$model,"train_n_repeat",NULL)

            tryCatch(stopifnot(is.integer(as.integer(config$model$train_n_repeat))),error = function(e){stop(simpleError(paste("parameter train_test_split from model configuration is not an integer but",config$model$train_n_repeat)))})
            tryCatch(stopifnot(config$model$train_n_repeat>=0 && config$model$train_n_repeat<=100),error = function(e){stop(simpleError(paste("Value must be between 0 and 100 and is",config$model$train_n_repeat)))})

            
            checker(config$model,"train_balance",NULL)
            checker(config$model,"train_balance",c("downSample","upSample","none"))
            #it has pass all the checks -> granted the class


            class(config$model) <- c("model_setup", config$model$model) # Either glmnet or Naive_Bayes thinking of subclass and inheritance - to read about
            config$model
        }
    }


    if(exists("model",where = model_configuration)){
        check_model(config)->config
    }else if(exists("models",where = model_configuration)){
        lapply(config$models,check_model)->config
        # It has passed all the checks, we can give it the multi_model class
        class(config) <- "multi_model_setup"
    }else{
        stop("the model configuration does not have neither a single model nor multiple - need to be a list element named either model or models")
    }

   config 

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

        ## Some dedicate treatment of the databases

## Phenotype DB: 
# Select the colum of interest (microorganism names and the Phenotypic trait) and get rid of everything else.
# Choose what to do with the NA at the phenotypic trait - enconded in the NA_substitution <NA2level|NonNA2level|None> and NA_substitution_value parameters of Phenotype structure

process_PhenotypeDB <- function(datos,config){
            # Remove the Phenotype columns we are not using here
    datos$Phenotype <- subset(datos$Phenotype,select=c(config$Phenotype$Bacterial_metadata_columns,config$Phenotype$Phenotypic_trait))
    fix_Phenotype_NA(datos,phenotypic_trait = config$Phenotype$Phenotypic_trait,method = config$Phenotype$NA_substitution,level = config$Phenotype$NA_substitution_value)
}

# The function change the NA for another value (unless set as none) - the replacement value can be chosen with NA_substitution_value

fix_Phenotype_NA <- function(datos,phenotypic_trait,method = "NA2level",level = "Unknown"){
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
# The Genome database may have undesired columns that so far can be remove by positive filtering (take all columns that fills a pattern rather than remove cols that fullfill a pattern)
# It remove duplicated genomes as well.
# A decision to be taken is about what to do with those microorganisms with more than one reference genomes - so far I have aggregated them by using max() but It deserves a look -
# particularly because is the slowest piece of code after the model itself.

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
    
    # slow code - delete when clear
    # Add and index to "find their row" and split the dataset in two types of columns for metdata only and genes only
    # Then split - apply - aggregate - put back
    #    tmp1<-cbind(tempIndex =c(1:nrow(datos)),datos[,!Gene_names(names(datos),configuration$Genome$Pattern_gene_columns)])    
    #    tmp2<-cbind(tempIndex =c(1:nrow(datos)),datos[, Gene_names(names(datos),configuration$Genome$Pattern_gene_columns)])
    #    tmp2 <- do.call(rbind,lapply(split(x = tmp2,f = tmp1$TaxID),function(x){data.frame(t(c(tempIndex = x[1,"tempIndex"],sapply(x[which(!names(x) %in% "tempIndex")],max,na.rm =T,simplify=T))))}))
    #    tmp2<- merge(tmp1,tmp2,by = "tempIndex",all = F)
    #    tmp2['tempIndex']<- NULL
    #    tmp2

    # the next code seems to be simpler but slower and create difficulties with the phenotypic names
    
    # I need a decision on what to do with the multiple genomes 
    # at the moment I consider using the maximum value per gene (since they are counting the number of repetitions of a single ortholog gene)
    #ddply(datos,"TaxID",function(x){
    #    number_var <- unlist(lapply(x,function(y){is.integer(y) | is.numeric(y)}))
    #    computed<-unlist(lapply(x[,number_var],max,na.rm = T))
    #    cbind(data.frame(x[1,!number_var]),data.frame(t(computed)))
    #    })
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


#####################################################
# Tier 3: Execute the models - functions  #
#####################################################

exec_model<-function(Predictor,Response,model_config){
    # Fit a model according to a model configuration file and some extra parameters
        #At the moment I am  using only functions at the caret package for modelling but some  changes would be needed if we resort to different machine learning wrapping packages
        # Therefore, I leave a "per model" kind of stream so it is easier to expand to other packages or functions
    
                tryCatch(expr = stopifnot(min(nrow(Predictor),ncol(Predictor))!=0),
                        error = function(e){stop("The x columns of the model are empty",call. =F)})

                tryCatch(expr = stopifnot(length(Response)!=0),
                          error = function(e){stop("The y columns of the model is empty",call. =F)})


    stopifnot(any(class(model_config) %in% c("glmnet","Naive_Bayes","XGBoost")))
    
    if(any(class(model_config) == "glmnet")){
        log(" Fitting a Glmnet Caret model")
        tryCatch({
                
                train(x = Predictor,
                      y = Response,
                      method = "glmnet",
                      trControl=trainControl(method = model_config$train_CV_method,number = model_config$train_n_repeat))->model
                log("Glmnet fitted")
                model
        },error = function(e){log("ERROR in Glmnet fit");e},silent = T)->result
    }
    if(any(class(model_config) == "Naive_Bayes")){
        log(" Fitting a Naive Bayes Caret model")
        tryCatch({ 
                train(  x = Predictor,
                        y = Response,
                        method = "nb",
                        trControl=trainControl(method = model_config$train_CV_method,number=model_config$train_n_repeat))->model
                        
                log("Naive Bayes fitted")
                model
        },error = function(e){log("ERROR in Naive Bayes fit");e},silent = T)->result
    } 
    if(any(class(model_config) == "XGBoost")){
        log(" Fitting a XGBoost model")
        tryCatch({ 
                train(  x = Predictor,
                        y = Response,
                        method = "xgbLinear",
                        trControl=trainControl(method = model_config$train_CV_method,number=model_config$train_n_repeat))->model
                        
                log("XGBoost model fitted")
                model
        },error = function(e){log("ERROR in XGBoost fit");e},silent = T)->result
    }
    result
}


send2modelling <- function(model_data,model_config,data_config){

    # Input: 
        #   model_data: data.frame (if an object would be a pheno-geno object) 8 columns for TaxID,GenomeID,Name, Genus, species, Strain, <Phenotype>, organism_name, <Gene 1>, <Gene 2>, ...,<Gene N>
        #   model_config: a list of values for the models - it must contain an element "model" with one of the available models, the rest are associate with the model
        #   data_config:  a list of values for organizing the data - these are object-like structures (for future conversion into classes)

        # Split data in validation - training - using caret

    if(any(class(model_config) == "model_setup")){
            
            exec_model(Predictor = model_data[,Gene_names(names(model_data),data_config$Genome$Pattern_gene_columns)],
                        Response = model_data[,data_config$Phenotype$Phenotypic_trait],
                        model_config = model_config)

    }else if(any(class(model_config) == "multi_model_setup")){
                        
        lapply(model_config,
                function(model_config,Predictor,Response){
                    exec_model(Predictor,Response,model_config)},
                            Predictor = model_data[,Gene_names(names(model_data),data_config$Genome$Pattern_gene_columns)],
                            Response = model_data[,data_config$Phenotype$Phenotypic_trait]) ->model
                            class(model) <- "multiple_train"
                            model
    }
}

balance_data_classes <- function(model_data,data_config,model_config){
    #theoption <- model_config$train_balance
    ifelse(test= class(model_config) == "multi_model_setup",theoption <- model_config[[1]][["train_balance"]],theoption <- model_config[["train_balance"]])
    match.arg(theoption,choices = c("downSample","upSample","none"),several.ok = F)
    switch(theoption,
        downSample = downSample(x = cbind(get_genome_columns(model_data),get_metadata_columns(model_data)),y = get_phenotype(datos),yname=data_config$Phenotype$Phenotypic_trait),
        upSample = upSample(x = cbind(get_genome_columns(model_data),get_metadata_columns(model_data)),y = get_phenotype(datos),yname=data_config$Phenotype$Phenotypic_trait))
    }

generate_train_validation_data <- function(model_data,model_config){

        if(any(class(model_config) == "multi_model_setup")){
            createDataPartition(get_phenotype(model_data),p=model_config[[1]]$train_test_split,list = F)->trainData
        }else{
            createDataPartition(get_phenotype(model_data),p=model_config$train_test_split,list = F)->trainData
        }

            list(validation_dataset = model_data[-trainData,],
            train_dataset =  model_data[trainData,] )
            
}
#####################################################
# Tier 4: interpret and save the models - functions #
#####################################################

evaluate_model <- function(model,validation_data,model_config,data_config){

    individual_model <- function(themodel,theconfig,thedata){
 
        if(class(themodel)[[1]] == "simpleError"){return(themodel)}

        if(any(class(theconfig) %in% c("glmnet","Naive_Bayes","XGBoost"))){
            predict(themodel,thedata)->theprediction
            #confusionMatrix(theprediction,as.factor(thedata[,data_config$Phenotype$Phenotypic_trait]))->CF
            if(nlevels(droplevels(get_phenotype(thedata)))>2){
               AUROC<- pROC::multiclass.roc(response = as.factor(get_phenotype(thedata)),predictor = as.numeric(theprediction)) 
            }else{
                AUROC <- pROC::roc(response = as.factor(get_phenotype(thedata)),predictor = as.numeric(theprediction))
            }

               list(CF = confusionMatrix(data = theprediction,reference = as.factor(get_phenotype(thedata)),mode = "everything")
               ,AUROC = AUROC)
        }
    }

    if(any(class(model) == "multiple_train")){
        structure(mapply(FUN = individual_model,model,model_config,MoreArgs = list(thedata = validation_data),SIMPLIFY = F),class = "multi_model_evaluation")
    
    }else{individual_model(model,model_config,validation_data)}
}

save_models <- function(config_file,folder_output,config_data,config_model,model){

       individual_model_case <- function(config_model,model){         
            stopifnot(any(class(config_model) %in% c("glmnet","Naive_Bayes","XGBoost")))
            
            if(config_model$model == "glmnet"){
                file_output <-gsub(x= config_file,pattern= "\\.dat",replacement = "\\.model.elasticnet.dat")
                log(paste("Saving results in ",file.path(folder_output,file_output)))
               #save(list = c("configuration","Documents","model"),file = file.path(folder_output,file_output),ascii = F)
            }
            if(config_model$model == "Naive_Bayes"){
                file_output <-gsub(x= config_file,pattern= "\\.dat",replacement = "\\.model.NaiveBayes.dat")
                log(paste("Saving results in ",file.path(folder_output,file_output)))
                #save(list = c("configuration","Documents","model_naiveBayes"),file = file.path(folder_output,file_output),ascii = F)
            }
            if(config_model$model == "XGBoost"){
                file_output <-gsub(x= config_file,pattern= "\\.dat",replacement = "\\.model.XGBoost.dat")
                log(paste("Saving results in ",file.path(folder_output,file_output)))
                #save(list = c("configuration","Documents","model_naiveBayes"),file = file.path(folder_output,file_output),ascii = F)
            }
        }

    if(any(class(config_model) %in% "model_setup")){
        individual_model_case(config_model,model)
    }else if(any(class(config_model)== "multi_model_setup")){
        mapply(FUN=individual_model_case,config_model,model)
    }

}

save_stats <- function(thedata,theevaluation,themodel,theconfig){
    #NOTE: The function extract_metrics put together the model set up and the results
        # This only work at the moment because glmnet, Naive_Bayes and XGBoost have the same options at the moment - I meant on making a single data.framefor all the models.
        # When new models are included, it is quite possible that this function need to keep the "stats" as a list and not as a data.frame (see the do.call)
    extract_metrics <- function(evaluation_data,model_config){

        if(class(evaluation_data)[[1]] == "simpleError"){
            thenames<-c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull",
            "AccuracyPValue", "McnemarPValue", "Sensitivity", "Specificity",
            "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1",
            "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy"
            )
            c(unlist(model_config),array(dim = length(thenames),dimnames = list(thenames)))
        }else{        
            data.frame( c(model_config,
                        as.list(evaluation_data$CF$overall),
                        as.list(unlist(evaluation_data$CF$byClass))),
                        AUROC = evaluation_data$AUROC[["auc"]][[1]])
        }
    }
   
    Phenotype_stats <- function(phenotype){
        Density <- round(table(phenotype,useNA ="always")/length(phenotype)*100,digits =2)
        Balance <- round(table(phenotype)/margin.table(table(phenotype))*100,digits = 2)
        Density_agg = sum(Density[!is.na(names(Density))])
        Balance_agg <- (range(Balance) %*% c(-1,1))/Density_agg * 100
        N = length(phenotype)
        list(Nphenotype = N,Density = Density,Density_agg = Density_agg,Balance = Balance,Balance_agg = Balance_agg)
    }
    Phenotype_stats(get_phenotype(thedata)) -> Stats_1
    ## The evaluation
    if(any(class(theevaluation) == "multi_model_evaluation") & any(class(themodel)== "multi_model_setup")){
        mapply(theevaluation,themodel,FUN = extract_metrics,SIMPLIFY = F)->evaluated
        as.data.frame(do.call(rbind,evaluated)) -> evaluated
    }else{extract_metrics(theevaluation,themodel)->evaluated}
    ### HERE : It may be a problem with the many models - if it is save a file per model.
    output <- cbind(evaluated,as.data.frame(as.list(unlist(Stats_1))),as.data.frame(as.list(unlist(theconfig,recursive=T))))
    file_output <-gsub(x= config_file,pattern= "\\.dat",replacement = "\\.results.dat")
    write.table(x = output,file =file.path(folder_output,file_output) ,sep = ";",dec = ",",row.names = F)
    #cat(RJSONIO::toJSON(list(Stats_1,evaluated,themodel,theconfig),pretty = T))
}
### Extras ###

#Little logging
log <- function(txt){
    print(paste("[",Sys.time(),":",txt,"]"))
}
# Some "getters"

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
    function(model_data){
        model_data[[phenotype]]
    }
}

###################################################################


######## VARIABLES ########################

args = commandArgs(trailingOnly=TRUE)
print(args)

if(length(args)  != 5){
    print("Not enough arguments provided, I carry on with presets")

    folder_config_file          <-  "/home/ubuntu/Models/GenePhene2/test.files"
    config_file                 <-  #"GenePhene2_Catalase_activity_KEGG_D2V50_pickone_genome.dat"
                                    #"GenePhene2_Catalase_activity_KEGG_D2V50_multiple_genomes.dat"
                                    #"GenePhene2_Catalase_activity_KEGG_BoW_multiple_genomes.dat"
                                    "GenePhene2_Catalase_activity_KEGG_BoW_pickone_genomes.dat"
    folder_model_config_file    <-  "/home/ubuntu/GenePhene2/Models/config.files"
    model_config_file           <-  "models.json"
                                    #"model.glmnet_elasticnet.json"
                                    #"model.Naive_Bayes.json"
                                    #"model.XGBoost.json"
    folder_output               <-  "/home/ubuntu/Models/GenePhene2/test.results"

#    folder_config_file <- "/home/ubuntu/Models/GenePhene2/data.files"
#    config_file <- "GenePhene2_Output_1,2,4-trihydroxybenzene_D2V_COG.dat"
#    folder_model_config_file <-   "/home/ubuntu/GenePhene2/Models/config.files"
#    model_config_file <- "model.glmnet_elasticnet.json"
#    folder_output <-  "/home/ubuntu/Models/GenePhene2/model.results"
    
    #folder_config_file <- "/home/ubuntu/Models/FAPROTAX/test.files"
    #config_file <- "FAPROTAX_photoautotrophy_KEGG.dat"
    #folder_model_config_file <- "/home/ubuntu/GenePhene2/Models/config.files/"
    #model_config_file <- "model.glmnet_elasticnet.json"
    #folder_output <- "/home/ubuntu/Models/FAPROTAX/test.results/"

#    folder_config_file <- "~/Models/GenePhene2/data.files/"
#    config_file <- "GenePhene2_Input_amino_acid_D2V_pFam.dat"
#    folder_model_config_file <- "../config.files/"
#    model_config_file <- "model.glmnet_elasticnet.json"
#    folder_output <- "~/Models/GenePhene2/model.results/"


}else{ 

    print(args);
    folder_config_file <- args[1]
    config_file <- args[2]
    folder_model_config_file <- args[3]
    model_config_file <- args[4]
    folder_output <-  args[5]
}

##################################### SCRIPT #######################
log("Starting the Script")
log("READING CONFIGURATION FILES ...")
#### Tier 1a : Get configuration files  #####

read_config_file(folder_config_file,config_file)->setup
Documents <- extract_DB_from_config_files(setup)
configuration <-extract_specific_config(setup)

read_config_file(folder_model_config_file,model_config_file)->model_configuration

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
model_configuration <- check_model_config(model_configuration)

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

log("CLASS BALANCING")

datos <- balance_data_classes(datos,configuration,model_configuration)

log("CLASS BALANCING FINISHED")

log("DATA SPLITTING")
#NOTE: It seems that generating the validation data from the artificially balanced dataset may modify its utility (I can't see why, I can't find where I read it)
with(generate_train_validation_data(datos,model_configuration),{validation_dataset <<- get("validation_dataset");datos <<- get("train_dataset")})

log("DATA SPLITTED")
#### Tier 3. Modelling ###

log("START MODELLING")

#Bypass: Few modifications to accelerate and to reduce errors
#datos[c(7,9:20)]->datos
#datos[configuration$Phenotype$Phenotypic_trait]<- as.factor(sample(x = c("yes","no"),size = nrow(datos),replace = T))

send2modelling(model_data = datos,model_config = model_configuration,data_config = configuration)-> model

log("END OF MODELLING STEPS")

# Evaluation of models
log("EVALUATION OF MODELS")

evaluate_model(model,validation_dataset,model_configuration,configuration)->model_evaluation

log("MODELS EVALUATED")
## Recording results
log("SAVING THE MODELS")

save_models(config_file,folder_output,configuration,model_configuration,model)
save_stats(thedata = datos, theevaluation = model_evaluation, themodel = model_configuration, theconfig = configuration) ->evaluated
log("MODELS SAVED")

log("END OF SCRIPT")

stop("END OF SCRIPT")
    # To make the model possible - PCA on the predictors - I think caret can do the same.
        #prcomp(x = datos[,Gene_names(names(datos),pattern_gene_columns)],center = F,scale. = F)->PCA
        #datos<- data.frame(datos[,phenotypic_trait],data.frame(PCA$x))   
        #names(datos)[1] <- phenotypic_trait
        #pattern_gene_columns <- "PC\\d+"



    
