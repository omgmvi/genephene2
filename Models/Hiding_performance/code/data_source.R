#### Odin M. Morón-García ####


##### Script to collect the data for the experiment ####


folder_output <- "/home/ubuntu/GenePhene2/Models/Hiding_performance/"
##### Phenotype DB ####

PhenotypesDB <- "/home/ubuntu/GenePhene2/Databases/Hiding/Metabolic.tsv"
Phenotypes <- c("catalase.positive","oxidase.positive","acid.phosphatase")
PK_1<- c("Name")


##### Taxonomy DB ####

TaxonomyDB <- "/home/ubuntu/GenePhene2/Taxonomy/Hiding/SpeciesTaxID.tsv"
PK_2 <- c("TaxID")


##### Metadata DB ####

MetadataDB <- "/home/ubuntu/GenePhene2/Metadata/Hiding/TaxID2Genome_metadata.tsv"
PK_3 <- "taxid"
PK_4 <- "species_taxid"


#### Genome DB ###

GenomeDB_COG <- "/home/ubuntu/GenePhene2/Genomes/Hiding/BagOfWords_Genome/HIDING_BoW_COG.tsv"
GenomeDB_KEGG <- "/home/ubuntu/GenePhene2/Genomes/Hiding/BagOfWords_Genome/HIDING_BoW_KEGG.tsv"
GenomeDB_pFam <- "/home/ubuntu/GenePhene2/Genomes/Hiding/BagOfWords_Genome/HIDING_BoW_pFam.tsv"
Geno_COG <- "Hiding_COG.tsv"
Geno_KEGG <- "Hiding_KEGG.tsv"
Geno_pFam <- "Hiding_pFam.tsv"

PK_5 <- "assembly_accession"
PK_6 <- "GenomeID"

############## JOIN DB #####
# Pheno - Taxon

read.table(PhenotypesDB,sep = "\t",quote ="",comment.char = "" ,header = T) -> PhenotypesDB_f
PhenotypesDB_f <- PhenotypesDB_f[,c(PK_1,Phenotypes)]

read.table(TaxonomyDB,sep = "\t",quote ="",comment.char = "" ,header = T) -> TaxonomyDB_f

merge(PhenotypesDB_f,TaxonomyDB_f,by = PK_1)->PhTaxDB

rm(PhenotypesDB_f,PhenotypesDB,TaxonomyDB,TaxonomyDB_f)

## PhenoTaxon - Metadata

read.table(MetadataDB,sep = "\t",quote ="",comment.char = "" ,header = T) -> MetadataDB_f
MetadataDB_f <- unique(MetadataDB_f)

# Care with the two columns TAXID and species_taxid

merge(PhTaxDB,MetadataDB_f,by.x = PK_2,by.y=PK_3)->Meta1
merge(PhTaxDB,MetadataDB_f,by.x = PK_2,by.y=PK_4)->Meta2
# Maybe just the second is enough since the TaxID that changes for strains of similar taxon remain with the species_id similar

# Care with the names after merging- The merge makes dissapear one of the names (merging by.x=TAXID  by.y = species_taxid make dissapear TAXID column)
valid_names <- intersect(names(Meta1),names(Meta2))
Meta1 <- Meta1[valid_names]
Meta2 <- Meta2[valid_names]

# just row binding the two df.
Meta <- rbind(Meta1,Meta2)

# Now, there's more than once the same genomes due to repetition in Meta1 and Meta2 / an antijoing may have worked to cut off the repeated ones. This is enough anyway
Meta[!duplicated(Meta$assembly_accession),] -> Meta

#clean up variables
rm(Meta1,Meta2,valid_names,MetadataDB_f,MetadataDB,PhTaxDB,PK_1,PK_2,PK_3,PK_4)


## PhenoTaxMetadata - Genome(s) ##

read.table(GenomeDB_COG,sep = "\t",quote ="",comment.char = "" ,header = T) -> GenomeDB_f
merge(Meta,GenomeDB_f,by.x= PK_5,by.y=PK_6)->Genome_f
write.table(x = Genome_f, file.path(folder_output,Geno_COG),sep = "\t",row.names = F)

read.table(GenomeDB_KEGG,sep = "\t",quote ="",comment.char = "" ,header = T) -> GenomeDB_f
merge(Meta,GenomeDB_f,by.x= PK_5,by.y=PK_6)->Genome_f
write.table(x =Genome_f,file.path(folder_output,Geno_KEGG),sep = "\t",row.names = F)

read.table(GenomeDB_pFam,sep = "\t",quote ="",comment.char = "" ,header = T) -> GenomeDB_f
merge(Meta,GenomeDB_f,by.x= PK_5,by.y=PK_6)->Genome_f
write.table(Genome_f,file.path(folder_output,Geno_pFam),sep = "\t",row.names = F)
