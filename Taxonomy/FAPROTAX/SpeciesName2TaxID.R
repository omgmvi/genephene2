# Getting the TaxID for the FAPROTAX Species .

#FAPROTAX contain Species names and strain names and I need to obtain the TaxID for them.
# In a first round I ignore the Strain type, to avoid the use of aproximate regular expressions
# There are two files in NCBI taxonomy new dump that provides TaxID code. The shortest and cleaner is typematerial.dmp.
# It contains just TaxID, Name (as Genus + Species) and Strain. So, it is the easiest to search (assuming no mispellings in any database)
# Those entries from FAPROTAX not found in this file are usually spelling errors or outdated names. In those cases, a manually written metadata-like table translates the mispelt to correct ones
# Some FAPROTAX entries are names that were not successfully accepted (and that was on the 1980's list). Those are in the metadata-like file with an NA in the traslation column.


######## SETUP ########

#FILE OPTIONS


PhenoDB_folder <- "/home/ubuntu/Databases/FAPROTAX/FAPROTAX_1.2.4"
PhenoDB_Species_file <- "FAPROTAX.tsv"
name_translation <- "known_mistakes_or_variants.txt"

TaxonomyDB_folder <- "/home/ubuntu/Genomes/Taxonomy"
TaxonomyDB_file <- "typematerial.dmp"
TaxonomyDB_long_file <- "names.dmp"

output_folder <- file.path(TaxonomyDB_folder,"FAPROTAX")
output_file <- "SpeciesTaxID.tsv"

#READING FILES 
#TO DO
   # I should check for exist and so on - not now

    #Phenotypic database 

# Read the whole database and keep only the columns corresponding to Name and Strain

SpDB <- read.table(file.path(PhenoDB_folder,PhenoDB_Species_file),sep = "\t",header = T)

# Faprotax column name is yet not well parsed and contain some strains within.
# The columns to work are the concatenation of Genus plus Species columns (sep by space) and a second column with the Strain name.
SpDB$Name <- trimws(paste(SpDB$Genus,SpDB$Species))
SpDB <- SpDB[,c("Name","Strain")]

    # NCBI Taxonomy database - At this moment just the typematerial.txt file - So, I can do  a left join with the well written species
read_TaxDB_types <- function(){
    TaxDB<- read.table(file.path(TaxonomyDB_folder,TaxonomyDB_file),sep="|",header = F,strip.white = T,comment.char ="",fill = T,quote = "")
    names(TaxDB) <- c("TaxID","Name","Type","Strain","V5")
    TaxDB$V5 <- NULL
    TaxDB
}

read_TaxDB_names <- function(){
    TaxDB<- read.table(file.path(TaxonomyDB_folder,TaxonomyDB_long_file),sep="|",header = F,strip.white = T,comment.char ="",fill = T,quote = "")
    names(TaxDB) <- c("TaxID","Name","UniqueName","NameClass")
    TaxDB$V5 <- NULL
    TaxDB
}
TaxDB <- read_TaxDB_types()
################### SCRIPT ################# 

###### Merging (inner-joining) Species names with Taxonomy

results <- list()

log1 <- nrow(SpDB)

##| Search 1: Merge Species names against typematerial Species names -It would work with correct spellings only |##


merge(SpDB,TaxDB,by = "Name",all = F)->results[[1]]
setdiff(SpDB$Name,unique(results[[1]]$Name))->missingSp

#Wrap up the search
merge(SpDB,data.frame(Name = missingSp))->SpDB
unique(results[[1]][c("Name","TaxID")]) -> results[[1]]

    #It has found 2722 species in the file, worthy to notice that there were a single TaxID for each bacterial name, so the Strain here wouln't have matter
print("log of round 1 - Searching Species DB against the typematerial database")
print(paste("the original database had",log1,"rows"))
print(paste("Missing Species",length(missingSp)))
print(paste("results then have",nrow(results[[1]])))
print(paste("and SpDB remain with",nrow(SpDB),"rows and cols"))
## Search 2 ## Correct the names of few mistakes 

read.table(file.path(TaxonomyDB_folder,'FAPROTAX',name_translation),sep = ";",header = F)->tmp
names(tmp)<-c("Name","CorrectName")
merge(tmp,TaxDB,by.x="CorrectName",by.y="Name",all =F)->tmp

merge(SpDB,tmp,all=F,by="Name")->results[[2]]
setdiff(SpDB$Name,unique(results[[2]]$Name))->missingSp

log1 <- nrow(SpDB)

merge(SpDB,data.frame(Name = missingSp))->SpDB
unique(results[[2]][c("Name","TaxID")]) -> results[[2]]


print("log of round 2 - Searching Species DB against the typematerial database")
print(paste("the original database had",log1,"row"))
print(paste("Missing Species",length(missingSp)))
print(paste("results then have",nrow(results[[2]])))
print(paste("and SpDB remain with",nrow(SpDB),"rows"))

## Search 3 ## The longer file
rm(TaxDB)
TaxDB <- read_TaxDB_names()


merge(SpDB,TaxDB,by = "Name",all = F)->results[[3]]

log1 <-nrow(SpDB)

setdiff(SpDB$Name,unique(results[[3]]$Name))->missingSp

merge(SpDB,data.frame(Name = missingSp))->SpDB

unique(results[[3]][c("Name","TaxID")]) -> results[[3]]

print("log of round 3 - Searching Species DB against the names database")
print(paste("the original database had",log1,"rows"))
print(paste("Missing Species",length(missingSp)))
print(paste("results then have",nrow(results[[3]])))
print(paste("and SpDB remain with",nrow(SpDB),"rows")
)
## Search 4 ## We keep searching on the long files - names.dmp but now we use the name corrector file


read.table(file.path(TaxonomyDB_folder,'FAPROTAX',name_translation),sep = ";",header = F)->tmp
names(tmp)<-c("Name","CorrectName")
merge(tmp,TaxDB,by.x="CorrectName",by.y="Name",all =F)->tmp

log1 <- nrow(SpDB)

merge(SpDB,tmp,all=F,by="Name")->results[[4]]
setdiff(SpDB$Name,unique(results[[4]]$Name))->missingSp


merge(SpDB,data.frame(Name = missingSp))->SpDB
unique(results[[4]][c("Name","TaxID")]) -> results[[4]]

print("log of round 4 - Searching Species DB against the typematerial database")
print(paste("the original database had",log1,"rows"))
print(paste("Missing Species",length(missingSp)))
print(paste("results then have",nrow(results[[4]])))
print(paste("and SpDB remain with",nrow(SpDB),"rows"))

read.table(file.path(TaxonomyDB_folder,'FAPROTAX',name_translation),sep = ";",header = F)->tmp5
names(tmp5) <- c("Name","CorrectName")
merge(data.frame(Name = missingSp),tmp5,by = "Name",all.x = F)

#In here I have to take a decision whether ignore the strain names or not of the Reference Genome. At the moment and for the sake of moving faster I will ignore the Strain of bacterias.
unique(Reduce(f=rbind, x=results))->results

# There are around 230 missing species where most of them are due to mistakes and mispellings and I will do a manual version for them

#### Export the table

write.table(results,file.path(output_folder,output_file),row.names = F, sep = "\t",quote =F )
