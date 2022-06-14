#### Odin M. Moron-Garcia - 10th May 2022 


## Code that read a genome, or containing a genome (and in fact any df with numerical columns)
## and encode colums as 0/1

###################
#### functions ####
###################

log <- function(...) cat(sprintf(paste0(...)), sep='\n', file=stderr())

#########################
#### Arguments parsing ##
#########################

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    print("usage: Rscript remove_duplicates_genomes.R ./Hiding_COG.tsv ./Hiding_COG.reduced.tsv ./Hiding_COG.presence.tsv 'COG\\d+\' GenomeID")
    stop("Bye Bye")
}

#completely careless command line input- no parsing

input_file <- args[1]
output_file <- args[2]
output_file_2 <- args[3]
regexp_columns <- args[4]
PK_column <- args[5] # Primary key column

# Some verbose output
log("input:",input_file)
log("output aggregated file:",output_file)
log("output 0/1 file :",output_file_2)
log("regexp to find columns:",regexp_columns)
log("Primary key column:",PK_column)

##################################
##### PROCESSING THE GENOMES #####
##################################


###### Reading the input ####
read.table(input_file,header = T,sep = "\t",comment.char = "")->input

####### Prepare the data for processing
## Split the data in metadata and genome - the only variable that will hold them together is the PK_column

names(input)-> nm
gene_nm <- grep(x=nm,pattern = regexp_columns,value = T)
meta_nm <- grep(x=nm,pattern = regexp_columns,value = T,invert=T)
Meta <- input[meta_nm]
Gene <- input[gene_nm]

row.names(Gene) <- Meta[[PK_column]]

rm(input,nm,gene_nm,meta_nm)


#########################
## 0/1 rows ##

### make every row a 0/1 for presence absence
apply(Gene,c(1,2),function(x)min(x,1))->Gene

## save it now before aggregating the genomes of a single species
Gene2<-merge(Meta,Gene,by.x = PK_column,by.y = "row.names")
write.table(Gene2,output_file_2,row.names=F,quote = F,sep ="\t")
rm(Gene2)
#########################


######### Group genomes of a single bacterial name #####

## some quicker solution to group genomes
list()->genomes

# notice the <<-; I am not using the apply nor the list as it should be (R and its lists ...)
apply(Meta,1,function(r){genomes[[ r["Name"] ]] <<- as.vector(c( genomes[[ r["Name"] ]],r["assembly_accession"]))})
# |1|: Here the names of the list are the names of the bacteria and later it will be the row.names (I believe?)

### Aggregate the single bacterial names
## The aggregated matrix - assume it is a 0/1 matrix if not, it will take the maximum of every column.
lapply(genomes,function(bac){if(!is.null(dim(Gene[bac,]))){apply(Gene[bac,],2,max)}else{Gene[bac,]}})->Gene
rm(genomes)
do.call(rbind,Gene)->Gene

###Prepare the matrix to write - now I need there's no duplicated metadata for the bacterias
# Wouldn't be enough a inner joing on the assembly_accession?
# Where did I change the row.names in Gene to be the name of the bacteria??? see |1|

unique(Meta[!names(Meta) %in% c("TaxID","assembly_accession","organism_name")])->Meta
Gene<-merge(Meta,Gene,by.x = "Name",by.y = "row.names")
rm(Meta)
write.table(Gene,output_file,row.names=F,quote = F,sep ="\t")
#######################################################
