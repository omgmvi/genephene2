# Odin Moron-Garcia June 2021
# Code to clean up the Hiding in Plain Sight database following the recomendations from the README.txt plus some extra considerations

#setwd("/home/ubuntu/Databases/Hiding")

####FUNCTIONS #####

#Some columns contain a range of temperature or pH as 5-7 that need to be converted in a single number or two colums. Here I choose the average 
# modified from the authors:

range2average <- function(col){
	col<-as.character(col)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
	col<-sapply(col, simplify=T, function(x){mean(as.numeric(unlist(strsplit(x, split="-", fixed=T))))})
	col
}



##### Read the data into a data.frame
# Authors had used lots of differentn strings as NA data that have to be sanitized - note also the check.names=F that brings column names with dots instead of spaces
ijsem<-read.delim("HidingPlainSightDB.txt", sep="\t", header=T, check.names=F, fill=T,
                  na.strings=c("NA", "", "Not indicated", " Not indicated","not indicated", "Not Indicated", "n/a", "N/A", "Na", "Not given", "not given","Not given for yeasts", "not indicated, available in the online version", "Not indicated for yeasts", "Not Stated", "Not described for yeasts", "Not determined", "Not determined for yeasts"))


### This error is platform dependant??? I keep it here just in case commented with 4#
#First column name is wrong due to a tab.
#### names(ijsem)<- names(ijsem)[-1]
#### ijsem$'NA' <- NULL

#FROM README

#simplify column names
colnames(ijsem)<-c("Habitat", "Year", "DOI", "rRNA16S", "GC", "Oxygen","Length", "Width", "Motility", "Spore", "MetabAssays", "Genus", "Species", "Strain", "pH_optimum", "pH_range", "Temp_optimum", "Temp_range", "Salt_optimum", "Salt_range", "Pigment", "Shape", "Aggregation", "FirstPage", "CultureCollection", "CarbonSubstrate", "Genome", "Gram", "Subhabitat", "Biolog")

#clean Habitat column
ijsem$Habitat <- as.factor(ijsem$Habitat)
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater (river, lake, pond)"]<-"freshwater"
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater sediment (river, lake, pond"]<-"freshwater sediment"

#clean Oxygen column
ijsem$Oxygen <- as.factor(ijsem$Oxygen)
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="aerobic"]<-"obligate aerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="anerobic"]<-"obligate anerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="microerophile"]<-"microaerophile"

#clean pH_optimum column
ijsem$pH_optimum<-range2average(ijsem$pH_optimum)

#remove pH values <0 and >10
ijsem$pH_optimum[ijsem$pH_optimum<0 | ijsem$pH_optimum>10]<-NA

#clean Temp_optimum column
ijsem$Temp_optimum<-range2average(ijsem$Temp_optimum)

#clean Salt_optimum column
ijsem$Salt_optimum<-range2average(ijsem$Salt_optimum)


#there are some formatting issues that should be solved
as.character(ijsem$rRNA16S)->ijsem$rRNA16S
as.character(ijsem$Subhabitat)->ijsem$Subhabitat

ijsem$GC <- as.character(ijsem$GC)
which(ijsem$GC=="DQ987877")->error
ijsem[error,"GC"] <- ijsem[error,"rRNA16S"]
ijsem[error,"rRNA16S"] <- "DQ987877"
rm(error) 

which(ijsem$GC == "BAOL01000001")->error
ijsem[error,"rRNA16S"]<- paste(ijsem[error,"rRNA16S"],ijsem[error,"GC"])
ijsem[error,"GC"] <- NA

which(ijsem$GC == "GU323338")->error
ijsem[error,"rRNA16S"]<- paste(ijsem[error,"rRNA16S"],ijsem[error,"GC"])
ijsem[error,"GC"] <- NA

which(ijsem$GC == "marine organism")->error
ijsem[error,"Subhabitat"]<- ijsem[error,"GC"]
ijsem[error,"GC"] <- NA
rm(error) 


ijsem$GC <- as.factor(ijsem$GC)
gsub(x=levels(ijsem$GC),pattern = "(\xe4\xf3\xf1)",replacement="-")->levels(ijsem$GC)
gsub(x=levels(ijsem$GC),pattern = "(\x8c\xb1.*)",replacement="")->levels(ijsem$GC)

# I take advantage of the range2average function to clean up the other variable that where not cleaned by the authors
ijsem$GC <- range2average(ijsem$GC)
ijsem$Length <- range2average(ijsem$Length)
ijsem$Width <- range2average(ijsem$Width)
#The next few variable may have no sense to get the average but to get them split in minimum and maximum - TO DO
ijsem$pH_range <- range2average(ijsem$pH_range)
ijsem$Temp_range <- range2average(ijsem$Temp_range)
ijsem$Salt_range <- range2average(ijsem$Salt_range)

# To conform with other databases in this study, a column Name being the aggregation of Genus and Species is required
ijsem$Name <- gsub(x = gsub(x = paste(trimws(ijsem$Genus),trimws(ijsem$Species),sep = " "),pattern ="(^\\b\\S+\\b) \\1 (\\b\\S+\\b)$",replacement = "\\1 \\2",perl = T,ignore.case = F),pattern = "(\\w+) (\\w+)",replacement = "\\1 \\L\\2",perl=T,ignore.case = F)
# There´s an specific problem with some compounds 2,3- ... as in 2,3 hydroxy something. Luckily, it only happen with the 2 prefix and is easyly clean up by substituting it 
# The issue will be soon covered, we want to split the string in CarbonSubstrate that is of the kind "Glucose,Fructose" in a vector and add every compound as a row (since it would break
# database normalization by having many NULL columns in some rows")

gsub(x=ijsem$CarbonSubstrate,pattern="(.*2),(.*)",replacement="\\1;\\2")->ijsem$CarbonSubstrate

# To make the vector and then the "long" format data.frame (also Entity-Atribute-Value format) I need to use the same technique than to make anti-joins in R without packages
# that is to make a column with indexes (may well be the row names), split the column entry with the carbon sources in vectors per 

# A function to make a colum of consecutive numeric ID's as primary key based on unique values of several colums
makeIDs <- function(df,cols){
	KeyPrim<-unique(subset(df,select=cols))
	ID <- 1:nrow(KeyPrim)
	cbind(ID,KeyPrim)
}

# A function to get a column whose values are a string that is itself a collection of strings separated by comma (hard coded at the moment). It does the split and transpose, giving
# a new row for each value in the string. That is, if a row with id 11 has "Glucose, Fructose, aminoacids" it will return 3 rows with id 11 and a column phenotrait with values Glucose then Gructose and then aminacids. 
# The function return the same data frame used to find the phenotypes and the IDs but now after a inner join on the phenotypes
split_phenolist<-function(df,pheno,ID){#browser()
	as.character(df[[pheno]])->df[[pheno]]
	transform(do.call(rbind,mapply(function(x,y){(cbind(cbind(phenoTrait=x),ID=y))},strsplit(subset(df,select=pheno,drop=T),split=","),subset(df,select=ID,drop=T))),ID=ID)->pheno_db
	na.omit(pheno_db)->pheno_db
	#print(paste("dim :",dim(pheno_db)))
	gsub(pattern="^ | $",x=pheno_db$phenoTrait,replacement="")->pheno_db$phenoTrait
	unique(pheno_db)->pheno_db
	#print(paste("dim :",dim(pheno_db)))
	merge(pheno_db,unique(df[,!names(df) %in% pheno]),all.y =F,all.x = F)->pheno_db
	#print(paste("dim :",dim(pheno_db)))
}

#### Application of the phenotypic list splitting to the metabolic assays column
# Using as Primary Key colums the Genus Species and Strain
KeyCols <- c("Genus","Species","Strain","Name")
phenotype_Col <- "MetabAssays"

# Generate ID´s
SpeciesDimensionTable<-makeIDs(ijsem,KeyCols)

##### Metabolic assays table #####

MetabAssays<-merge(subset(ijsem,select=c(phenotype_Col,KeyCols)),SpeciesDimensionTable)
# Get phenotypes by splitting the strings
MetabolicTraits<-split_phenolist(MetabAssays,phenotype_Col,KeyCols)

# And in order to make it usable I put it back to wide format
MetabolicTraits$default_val <- "yes"
MetabolicTraits<-(reshape2::dcast(data=MetabolicTraits,formula= Genus + Species + Strain + Name + ID ~ phenoTrait,value.var= "default_val"))
MetabolicTraits$default_val <- NULL
MetabolicTraits$ID <- NULL

##### Carbon Substrate table ####
# Almost the same code for the column Carbon Substrate (won't refactor at the moment)

phenotype_Col <- "CarbonSubstrate"

Carbon<-merge(subset(ijsem,select=c(phenotype_Col,KeyCols)),SpeciesDimensionTable)

CarbonSubstrates<-split_phenolist(Carbon,phenotype_Col,KeyCols)
CarbonSubstrates$default_var = "yes"
CarbonSubstrates<-(reshape2::dcast(data=CarbonSubstrates,formula= Genus + Species + Strain + Name + ID ~ phenoTrait,value.var="default_var" ))
CarbonSubstrates$default_var <- NULL  
CarbonSubstrates$ID <- NULL  

############### OUTPUT #########

# Export all files out to tsv formats
write.table(ijsem, file = "ijsem.tsv",row.names = F,sep = "\t",quote = F)#Almost the same than the original one but with some mistakes corrected
write.table(CarbonSubstrates, file = "Substrates.tsv",row.names = F,sep = "\t",quote = F)#wide format phenotypes related with substrates used
write.table(MetabolicTraits, file = "Metabolic.tsv",row.names = F,sep = "\t",quote = F)# wide format phenotypes related with metabolic reactions
