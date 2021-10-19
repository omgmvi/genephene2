# Getting the TaxID for the Hiding on plaing sight  Species .

#Hiding contain Species names and strain names and I need to obtain the TaxID for them.
# In a first round I ignore the Strain type, to avoid the use of aproximate regular expressions
# There are two files in NCBI taxonomy new dump that provides TaxID code. The shortest and cleaner is typematerial.dmp.
# It contains just TaxID, Name (as Genus + Species) and Strain. So, it is the easiest to search (assuming no mispellings in any database)
# Those entries from HIDING not found in this file are usually spelling errors or outdated names. In those cases, a manually written metadata-like table translates the mispelt to correct ones


######## SETUP ########

#FILE OPTIONS


PhenoDB_folder <- "/home/ubuntu/Databases/Hiding"
PhenoDB_Species_file <- "ijsem.tsv"
name_translation <- "known_mistakes_or_variants.txt"

TaxonomyDB_folder <- "/home/ubuntu/Genomes/Taxonomy"
TaxonomyDB_file <- "typematerial.dmp"
TaxonomyDB_long_file <- "names.dmp"

output_folder <- file.path(TaxonomyDB_folder,"Hiding")
output_file <- "SpeciesTaxID.tsv"

#READING FILES 
#TO DO
   # I should check for exist and so on - not now

    #Phenotypic database 

# Read the whole database and keep only the columns corresponding to Name and Strain

SpDB <- read.table(file.path(PhenoDB_folder,PhenoDB_Species_file),sep = "\t",header = T,quote = "",comment.char = "")
SpDB <- SpDB[,c("Name","Strain")]
# I need to process the Name column 
    #remove double spaces and trim whitespaces
    #SpDB$Name <- gsub(x = gsub(x = trimws(SpDB$Name),pattern = "\\s{2,}",replacement = " "),pattern = "(\\b\\w*\\b)\\1",replacement = "\\1",perl = T,ignore.case = T)

    # NCBI Taxonomy database 

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

# At this moment just the typematerial.txt file - So, I can do  a left join with the well written species
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
read.table(file.path(output_folder,name_translation),sep = ";",header = F)->tmp
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

read.table(file.path(output_folder,name_translation),sep = ";",header = F)->tmp
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

# Just a check - could be well deleted
read.table(file.path(output_folder,name_translation),sep = ";",header = F)->tmp5
names(tmp5) <- c("Name","CorrectName")
merge(data.frame(Name = missingSp),tmp5,by = "Name",all.x = F)

#In here I have to take a decision whether ignore the strain names or not of the Reference Genome. At the moment and for the sake of moving faster I will ignore the Strain of bacterias.
unique(Reduce(f=rbind, x=results))->results

#### Export the table

write.table(results,file.path(output_folder,output_file),row.names = F, sep = "\t",quote =F )
