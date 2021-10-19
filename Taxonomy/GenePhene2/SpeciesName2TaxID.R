# Getting the TaxID for the GenePhene2 Species .


######## SETUP ########

#FILE OPTIONS

PhenoDB_folder <- "/home/ubuntu/Databases/GenePhene2"
PhenoDB_Species_file <- "list_of_species.txt" 
TaxonomyDB_folder <- "/home/ubuntu/Genomes/Taxonomy"
TaxonomyDB_file <- "typematerial.dmp"

output_folder <- file.path(TaxonomyDB_folder,"GenePhene2")
output_file <- "SpeciesTaxID.tsv"
#Reading files
   # I should check for exist and so on - not now
SpDB <- read.table(file.path(PhenoDB_folder,PhenoDB_Species_file),sep = "\t",header = T)
TaxDB<- read.table(file.path(TaxonomyDB_folder,TaxonomyDB_file),sep="|",header = F,strip.white = T,comment.char ="",fill = T,quote = "")
names(TaxDB) <- c("TaxID","Name","Type","Strain","V5")
TaxDB$V5 <- NULL


## Search 1 ##
merge(SpDB,TaxDB,by = "Name",all = F)->results
#It is worthless to solve the problem here with the Strains - it is visible that the Strain types are the same
# And only three species are missin
setdiff(SpDB$Name,unique(results$Name))->missingSp
# The issue is a misspelling -> Dechloromonas denitrificans
# A mix of Name and  Strain -> Thauera aromatica AR-1
# Bacteroides merdae is nos Parabacteroides merdae and it would be only found if using the name.dmp file
merge(data.frame(Name = c("Dechloromonas denitrifican","Thauera aromatica AR-1","Bacteroides merdae"),Real_Name = c("Dechloromonas denitrificans","Thauera aromatica","Parabacteroides merdae")),TaxDB,all = F,by.x = "Real_Name",by.y = "Name")->missingSp

#In here I have to take a decision whether ignore the strain names or not of the Reference Genome. At the moment and for the sake of moving faster I will ignore the Strain of bacterias.

rbind(results[c("Name","TaxID")],missingSp[c("Name","TaxID")])->results
unique(results)->results


#### Export the table

write.table(results,file.path(output_folder,output_file),row.names = F, sep = "\t",quote =F )
