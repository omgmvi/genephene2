# 24th August - Odin Moron-Garcia

# Important error correction on 11 October


#Goal:
# Simple script to generate a table that 'translate' TaxIDs to Reference genomes 

# Our phenotypic databases have entries with either bacterial names or bact. names and type strains. 
# Using the ncbi taxonomy database we have made a list -per phenotype database - that translate those names to ncbi tax_id (in the taxonomy database files there's a single tax_id per entry)
# On the other hand, Seb's pipelines operate in reference genomes from the Refseq database with the curious fact that Refseq db has two types of Taxon ID - TaxID and Species_TaxId
# Since they seem to be related, it happens that most entries have the same taxid and species taxide (as a single entry) but some has a different taxid as if they were moved entries (from taxa revisions or any other type fo change.

# Therefore, here we collect as many Refseq whose tax_id OR species_tax_id coincides with the corresponding Taxonomy DB taxid (per phenotypic database like MDB, FAPROTAX, MIDAS ...)

#On the week before 11th October I found that not all genomes were being collected (I do not really  need all genomes, but I need to be able to interact with other pre-stablished databases
# like the first genephene - meaning I need to match species phenotypes to their genomes, whatever they have used.



#########################
######### SETUP #########
#########################

### DATABASE SELECTION ###
## Select a Phenotypic Databases to set the variables up
analysis <- ""
while(!analysis %in% c(1:4)){

    cat(prompt="Choose the databases\n\t1)GenePhene2\n\t2)FAPROTAX\n\t3)Hiding\n\t4)MDB\n")
    analysis <- readLines("stdin",n=1)
}

analysis <- c("GenePhene2","FAPROTAX","Hiding","MDB")[as.numeric(analysis)]


# Printing the selection - for testing purposes
print("Database:")
print(analysis)


#analysis <- "GenePhene2"
#analysis <- "FAPROTAX"


### VARIABLES ###

if(analysis == "GenePhene2"){
    Taxonomy_folder <- "/home/ubuntu/Taxonomy/GenePhene2"
    Taxonomy_file <- "SpeciesTaxID.tsv"

    output_folder <-  "/home/ubuntu/Metadatada/GenePhene2"
    output_file_genomes   <-  "genome_metadata.tsv"
    output_file_models   <-  "TaxID2Genome_metadata.tsv"
}

if(analysis == "FAPROTAX"){

    Taxonomy_folder <- "/home/ubuntu/Taxonomy/FAPROTAX"
    Taxonomy_file <- "SpeciesTaxID.tsv"
    
    output_folder <-  "/home/ubuntu/Metadata/FAPROTAX"
    output_file_genomes   <-  "genome_metadata.tsv"
    output_file_models   <-  "TaxID2Genome_metadata.tsv"
}

if(analysis == "Hiding"){

    Taxonomy_folder <- "/home/ubuntu/Taxonomy/Hiding"
    Taxonomy_file <- "SpeciesTaxID.tsv"
    
    output_folder <-  "/home/ubuntu/Metadata/Hiding"
    output_file_genomes   <-  "genome_metadata.tsv"
    output_file_models   <-  "TaxID2Genome_metadata.tsv"
}


if(analysis == "MDB"){

    Taxonomy_folder <- "/home/ubuntu/Taxonomy/MDB"
    Taxonomy_file <- "SpeciesTaxID.tsv"
    
    output_folder <-  "/home/ubuntu/Metadata/MDB"
    output_file_genomes   <-  "genome_metadata.tsv"
    output_file_models   <-  "TaxID2Genome_metadata.tsv"
}

# REFSEQ DATABASE - independent of the phenotype

refseq_folder <- "/home/ubuntu/Metadata/refseq"
refseq_file   <- "assembly_summary_refseq.txt"


### FUNCTIONS #####

#Wrappers 
read_data <- function(folder,file,sep = "\t",header = T,comment.char ="",quote="",...){
    read.table(file.path(folder,file),header = header,sep = sep,comment.char = comment.char,quote = quote,...)
}


read_refseq <- function(folder,file){
    read_data(folder,file,comment.char = "",header = F,nrow = 1,skip = 1,stringsAsFactor = F)->header
    header[1] <- "assembly_accession"

    read_data(folder,file,header =F,skip = 2)-> Refseq
    names(Refseq) <- header
    
    Refseq[,c("taxid","species_taxid","organism_name","ftp_path","assembly_accession")]
}

write_data<- function(x,folder,file){
    write.table(x,file.path(folder,file),sep = "\t",row.names = F, col.names = T,quote = F)
}


##################
##### SCRIPT #####
##################

## Load the data
read_data(Taxonomy_folder,Taxonomy_file)->Taxonomy
read_refseq(refseq_folder,refseq_file)->Refseq

#Inner Join - 

merge(Refseq,Taxonomy,by.x = "taxid",by.y= "TaxID",all.x = F,all.y = F)->result_A
merge(Refseq,Taxonomy,by.x = "species_taxid",by.y= "TaxID",all.x = F,all.y = F)->result_B
rbind(result_A,result_B)->result

## write a table for Seb with the FTP address and another one for the modelling with just taxid and  

write_data(result[c("assembly_accession","taxid","species_taxid","ftp_path")],output_folder,output_file_genomes)
write_data(result[c("assembly_accession","taxid","species_taxid","organism_name")],output_folder,output_file_models)


