# Getting the TaxID for Methanogens databse .
# It is a bit primitive - similar to the one in GenePhene2 with the mistakes or missing names solved "manually" hardcoded in a two-column data.frame (rather than using a file like in FAPROTAX or Hiding
# Yet a single species was only on the names.dmp file, so, I had to read the file again .- yet, I do not do any function or similar

######## SETUP ########

#FILE OPTIONS

PhenoDB_folder <- "/home/ubuntu/Databases/MDB"
PhenoDB_Species_file <- "list_of_species.txt" 
TaxonomyDB_folder <- "/home/ubuntu/Genomes/Taxonomy"
TaxonomyDB_file <- "typematerial.dmp"
TaxonomyDB_long_file <- "names.dmp"

output_folder <- file.path(TaxonomyDB_folder,"MDB")
output_file <- "SpeciesTaxID.tsv"


#Reading files
   # I should check for exist and so on - not now
SpDB <- read.table(file.path(PhenoDB_folder,PhenoDB_Species_file),sep = "\t",header = T)
TaxDB<- read.table(file.path(TaxonomyDB_folder,TaxonomyDB_file),sep="|",header = F,strip.white = T,comment.char ="",fill = T,quote = "")
names(TaxDB) <- c("TaxID","Name","Type","Strain","V5")
TaxDB$V5 <- NULL

## Search 1 ##
merge(SpDB,TaxDB,by = "Name",all = F)->results

# And only three species are missin
setdiff(SpDB$Name,unique(results$Name))->missingSp
# Remove the species found
merge(SpDB,data.frame(Name = missingSp))->SpDB

merge(data.frame(Name = c("Methanohalophilus euhalobius","Methanomicrococcus blatticola","Methanoplanus petrolearius"),Real_Name = c("Methanohalophilus euhalobius","Methanimicroccus blatticola","Methanolacinia petrolearia")),TaxDB,all = F,by.x = "Real_Name",by.y = "Name")->result_tmp

#In here I have to take a decision whether ignore the strain names or not of the Reference Genome. At the moment and for the sake of moving faster I will ignore the Strain of bacterias.
rbind(results[c("Name","TaxID")],result_tmp[c("Name","TaxID")])->results
unique(results)->results

# And only a single species are missing
setdiff(SpDB$Name,unique(results$Name))->missingSp
# Remove the species found
merge(SpDB,data.frame(Name = missingSp))->SpDB


#Search on the names.dmp DB - larger one.

TaxDB<- read.table(file.path(TaxonomyDB_folder,TaxonomyDB_long_file),sep="|",header = F,strip.white = T,comment.char ="",fill = T,quote = "")
names(TaxDB) <- c("TaxID","Name","UniqueName","NameClass")
TaxDB$V5 <- NULL

merge(SpDB,TaxDB,by = "Name",all = F)-> result_tmp

rbind(results[c("Name","TaxID")],result_tmp[c("Name","TaxID")])->results
unique(results)->results

# And no  species are missing
setdiff(SpDB$Name,unique(results$Name))->missingSp
# Remove the species found
merge(SpDB,data.frame(Name = missingSp))->SpDB
#### Export the table

write.table(results,file.path(output_folder,output_file),row.names = F, sep = "\t",quote =F )
