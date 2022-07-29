#!/usr/bin/env Rscript

# Odin Moron-Garcia
# Developed from
# 23 July 2021 - continuation on 20-30 August
# 30 Sept. 2021 : Update documentation 
# 18 July 2022 - Separation in two or more scripts
#       1. reading JSON onto a 3 tables files (18-20 July 20220)
#       2. fitting model (22- July-2022)

# NOTES:
# At the moment I uses S3 classes to take advantage of the type system in R, that is, I am using class(object) <- "type" to tag wether a variables is of certain type or other (for example if it is a ssingle model or multiple models)

#TODO: 
    ## Define the output - almost
    ## Define what to do with the errors.

###################### SETUP ###################

## S3 classes definitions
# Class model_setup
    #subclases (or 2nd classes) glmnet, Naive_Bayes and XGBoost.
    # Both types are a list containing package,model, train_CV_method <CV|LOOCV>, train_test_split <positive numeric>

## Class multi_model_setup
    # A list of many model_setup objects

# Class multi_train a list of several train objects (train class is from caret package see methods(class = "train")

##### LIBRARIES #####
require(rjson)
require(plyr) # use of dlply and ddply to summarize values

require(caret) # For models, in particular partition of data

#find.package("pROC",quiet =F)
#find.package("ROSE",quiet =F)
#find.package("smotefamily",quiet =F)

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



########################
# Tier 1: Read the DBs #
########################


read_DB <- function(folder,file){
    if(file.exists(file.path(folder,file))){
        read.table(file.path(folder,file),sep = ";")
    }else{stop(paste("File not found",file.path(folder,file)))}
}

read_DBs <-function(file,fPheno,fGeno,fMeta =NULL){
    Dbs <- list()
    Dbs$Phenotype <- read_DB(folder,fPheno)
    Dbs$Genome <- read_DB(folder,fGeno)
    if(!is.null(fMeta)){Dbs$Metadata <- read_DB(folder,fMeta)}
    Dbs
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

fix_balance_data_classes <- function(model_config){
    
    if(any(class(model_config) == "multi_model_setup")){

        theoption <- model_config[[1]][["train_balance"]]
        warning("In a multi model experiment the balance is done according to the options set at the first model.\n
                The configuration structure is changed to reflect that")
        cl <- class(model_config)
        lapply(model_config,function(model_struct){model_struct[["train_balance"]] <- theoption;model_struct})->model_config
        class(model_config) <- cl
    }
    model_config
}

balance_data_classes <- function(model_data,data_config,model_config){

    ifelse(test= any(class(model_config) == "multi_model_setup"),
                theoption <- model_config[[1]][["train_balance"]],
                theoption <- model_config[["train_balance"]])

    match.arg(theoption,choices = c("downSample","upSample","None"),several.ok = F)
    switch(theoption,
        downSample  = downSample(x = cbind(get_genome_columns(model_data),get_metadata_columns(model_data)),y = get_phenotype(datos),yname=data_config$Phenotype$Phenotypic_trait),
        upSample    = upSample(x = cbind(get_genome_columns(model_data),get_metadata_columns(model_data)),y = get_phenotype(datos),yname=data_config$Phenotype$Phenotypic_trait),
        None        = model_data)
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

        #Another tinkering step - what tod with the errors at modelling?
        if(any(class(themodel) == "model_setup") & any(class(theevaluation) == "error")){return()}
        if(any(class(theevaluation) == "multi_model_evaluation")){
            cl <- class(theevaluation) 
            lapply(theevaluation,function(x){if(any(class(x)=="error")){}else{x}})->theevaluation
              class(theevaluation) <- cl 
        }
        if(all(unlist(lapply(theevaluation,is.null)))){return()}

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
    }else{
        extract_metrics(theevaluation,themodel)->evaluated
    }
        
        ### HERE : It may be a problem with the many models - if it is save a file per model.
        output <- cbind(evaluated,as.data.frame(as.list(unlist(Stats_1))),as.data.frame(as.list(unlist(theconfig,recursive=T))))


        # Let's build the file name in a way that provide some information of what is inside and also a Date code  to avoid any ambiguity that may happen (repetition of the same model with the same data in two instances, for example)

        if(any(class(themodel) == "model_setup")){
            Model   <- themodel$model
        }else if(any(class(themodel) == "multi_model_setup")){
            Model   <- "multi_model_setup"
        }
        Database        <- theconfig$Database
        Trait           <- theconfig$Phenotype$Phenotypic_trait
        Trait           <- paste(strsplit(Trait,split=" ")[[1]],collapse = "-")

        Genome_type     <- theconfig$Genome$Genome_type
        Ortholog_type   <- theconfig$Genome$Genome_ortholog_type
 
        Date            <- format(Sys.time(),"%b%d-%H:%M:%S")
        file_output <- paste(paste(Database,Model,Trait,Genome_type,Ortholog_type,sep = "_"),Date,"results","dat",sep = ".")
        
        log(paste("Saving evaluation results in ",file.path(folder_output,file_output)))
        #file_output <-gsub(x= config_file,pattern= "\\.dat",replacement = "\\.results.dat")
        write.table(x = output,file =file.path(folder_output,file_output) ,sep = ";",dec = ",",row.names = F)
        #cat(RJSONIO::toJSON(list(Stats_1,evaluated,themodel,theconfig),pretty = T))
}
### Extras ###

#Little logging
log <- function(txt){
    print(paste("[",Sys.time(),":",txt,"]"))
}

###################################################################


######## VARIABLES ########################

args = commandArgs(trailingOnly=TRUE)
log(paste("arguments:",args))


if(length(args)  != 7){
    print("Not enough arguments provided, I carry on with presets")

    folder                      <-  "/home/ubuntu/GenePhene2/Rpackage/Hiding_DB3/"
    Pheno_file                  <-  "Phenotype.csv"
    Geno_file                   <-  "Genome.csv"
    Meta_file                   <-  "Metadata.csv"
    folder_model_config_file    <-  "/home/ubuntu/GenePhene2/Models/config.files"
    model_config_file           <-  "models_downSample.json"
                                    #"model.glmnet_elasticnet.json"
                                    #"model.Naive_Bayes.json"
                                    #"model.XGBoost.json"
    folder_output               <-  "/home/ubuntu/GenePhene2/Rpackage/Hiding_DB3/"



}else{ 

    folder <- args[1]
    Pheno_file <- args[2]
    Geno_file <- args[3]
    Meta_file <- args[4]
    folder_model_config_file <- args[5]
    model_config_file <- args[6]
    folder_output <-  args[7]
}

##################################### SCRIPT #######################
log("Starting the Script")
log("READING DATA FILES")
read_DBs(folder,Pheno_file,Geno_file,Meta_file)->DBs
log("DATA FILES READ")
log("READING MODEL FILE")
read_config_file(folder_model_config_file,model_config_file)->model_configuration
log("MODEL FILE READ")
log("CLASS BALANCING")
#Fix the case for multi model experiments, I can't manage each balance separately (this script is supposed to be for testing the same dataset with several models

model_configuration <- fix_balance_data_classes(model_configuration)
datos <- balance_data_classes(datos,configuration,model_configuration)

log("CLASS BALANCING FINISHED")

log("DATA SPLITTING")
#NOTE: It seems that generating the validation data from the artificially balanced dataset may modify its utility (I can't see why, I can't find where I read it)
with(generate_train_validation_data(datos,model_configuration),{validation_dataset <<- get("validation_dataset");datos <<- get("train_dataset")})
log("DATA SPLITTED")
#tryCatch(stopifnot(nrow(validation_dataset) !=0),error = function(e){log("VALIDATION DATA EMPTY - LEAVING THE PROGRAM");stop()})
log("CHECKPOINT - VALIDATION DATA NOT EMPTY - PASSED")
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
#The models are not really being saved to avoid space problme - actually I just want the parameters
save_models(config_file,folder_output,configuration,model_configuration,model)
save_stats(thedata = datos, theevaluation = model_evaluation, themodel = model_configuration, theconfig = configuration) ->evaluated
log("MODELS SAVED")

log("Tagging the input file")
# I do not like this way doing this but at the moment for the sake of finishing this.
# I need to know which tasks has been finished or not - due to some error and I do it by moving the files
file.rename(file.path(folder_config_file,config_file),to = file.path(folder_config_file,paste(config_file,"complete",sep = ".")))
log("END OF SCRIPT")

#stop("END OF SCRIPT")
    # To make the model possible - PCA on the predictors - I think caret can do the same.
        #prcomp(x = datos[,Gene_names(names(datos),pattern_gene_columns)],center = F,scale. = F)->PCA
        #datos<- data.frame(datos[,phenotypic_trait],data.frame(PCA$x))   
        #names(datos)[1] <- phenotypic_trait
        #pattern_gene_columns <- "PC\\d+"



    
