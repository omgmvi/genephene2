# Code to make the MDB.csv to compy with my own business rules.
# All tables are tabulation separated values
# All table has unquoted strings
# All strains must be in a single field, even if that breaks Codd's Normal Form (A single line records is used in bioinfo)
# Strain names will be duplicated as 1 field called Name and 2 fields Genus and Species
# If the original database have any record with multiple data in, usually corresponding to several substrates, it will be split in several columns


##### Helper variables to split the table by name and type of data.
### Note: These variables will change throughout the script to adapt to new variables being generated

METADATA <- c("Name", "Main publication", "Other publication", "Type strain", "DSM strain number",
	      "ATTC strain number", "Other collections", "Description", "Link","Max. cell GC content")

NUMERIC <- c("Growth doubling time [h]","Growth rate", "Min. growth temp.", "Max. growth temp.", "Min. optimal growth temp.",
	     "Max. optimal growth temp.", "Min. cell GC content",
	     "Min. cell length", "Max. cell length", "Min. cell width", "Max. cell width",
	     "Min. growth NaCl", "Max. growth NaCl", "Min. optimal growth NaCl",
	     "Max. optimal growth NaCl", "Min. growth pH", "Max. growth pH",
	     "Min. optimal growth pH", "Max. optimal growth pH")

METABOLIC <- c( "Cell shape", "Gram reaction","Motility",
	       "1-butanol", "2-butanol", "Isobutanol", "2-propanol", "Acetate",
	       "Additional growth requierments", "Butanol", "Carbon Monoxide",
	       "Cyclopentanol", "Dimethylamine", "Dimethyl sulfide", "Ethanol",
	       "Formate collections", "H2+CO2", "H2+methanol", "Methanol", "Methylamine",
	       "Other substrate", "Propanol", "Propionate", "Trimethylamine","Min. growth requierments")

ENVIRONMENTS <- c("Intestinal tracks", "Other environment", "Reactor environment","Soil environment", "Volcanic environment", "Water environment")

#setdiff(names(datos),c(ENVIRONMENTS,METABOLIC,NUMERIC,METADATA))

# Strings that will be substituted by NA, "yes" or "no". Special care has to be taken with the 0 as no for the numerical phenotypes since 0 is a legitimate number there

NA_STRINGS <- c("not determined","unknown (cell lysis)",
		"not applicable","unavailable","nd",
		"variable","None","none",'-1',"indifferent","no data")
YES_STRINGS <- c("yes", "yes ","positive","positiv","feeds","feedsa")
NO_STRINGS  <- c("0","no","do not feed on","negative")


##### Helper functions #######

exchange_strings <- function(df,str,replacement){
	# Search by column for those records valued as the string (whole string, not a substring)
	# and once found, replace the whole string with the new one.
	apply(df,2,function(x){x[x %in% str] <- replacement;x})->result
	data.frame(result,check.names = F)
}
save_file <- function(df,filename){
	#wrapper to write files without quotes and as tsv
	write.table(df,filename,sep = "\t",row.names = F,quote = F)
}
check_non_numeric_value <- function(df){
	# Search cell by cell if a factor is a number or has been coerced to NA.
	# In that way, I get the non numerical values back. Notice that just unique values are returned
	apply(df,c(1,2),function(x){
		      if(is.na(as.numeric(as.character(x)))){x}})->tmp;
	unique(unlist(tmp))->tmp
	tmp
}

exchange_strings_pattern <- function(str_array,subst_pattern,replacement,ignore.case =T){
	#This function look into a string array (note that it is different than exchange string that worked in a data.frame searching colum by colum
	# It replace any substring given by the pattern with a replacement. Note that it does not use Perl regexp with the idea of using simple regex
	gsub(x = str_array,pattern = subst_pattern,replacement = replacement,ignore.case=ignore.case)->replaced_str
	replaced_str
}


########################## Fixing the databases #########
# The code fix the issues with the databases (mostly erratas) and split fields with more than one value into many columns, keeping the number of rows constant
# Therefore I have to update the names in the helper arrays and keep an ID to avoid row movement problems (mostly the ID will be used during the field string splitting.
# I am separating the data into 4 databases according to the type of phenotypes it records.

# Metadata contains Strain names, and descriptions. Notice that a variable Max. cell GC content is also included due to be a text variable although it was expected to be a numerical. I well could have put it in the metabolic phenotypes. Anyway that variable is useless for any purpose.
# Environmental variable contains phenotypes regarding whether a bacteria lives or not in certain env. Won't be use neither in this project.

# numerical_phenotypes contain phenotypes concerning temperatures and Ph. and therefore numbers. The effort here is to keep in case in the future we use these values 
# but initially we will do classification. It needed a task of cleaning NA's and non-numerical values

# metabolic_phenotypes contain substrate utilization by the methanogenics. It require lots of processing for messy uncontrolled vocabulary plus the aditional and minimum 
# requirements being code as comma separated strings.

#### Read original data ####

read.csv("MDB.csv",check.names = F,)->datos
ID <- 1:NROW(datos)


### DB 1: METADATA & Max. cell GC content

metadata <- datos[,METADATA]

# Substitution of string for NA
metadata <- exchange_strings(metadata,NA_STRINGS,NA)
metadata <- setNames(metadata,nm = METADATA)
## Max. cell GC content has more values to be substituted by NA or yes and I do with manual care
exchange_strings(metadata["Max. cell GC content"],"0",NA)->metadata["Max. cell GC content"]
exchange_strings(metadata["Max. cell GC content"],"present","yes")->metadata["Max. cell GC content"]
exchange_strings(metadata["Max. cell GC content"],"absent","no")->metadata["Max. cell GC content"]


#I need to put all Strain type strings together in a single variable. I do it by "folding" a paste function
Reduce(f = function(x,y){paste(x,y,sep = ",")},
       x = metadata[,c("Type strain","DSM strain number","ATTC strain number","Other collections")],
       right=F,accumulate=F)->metadata$Strain

METADATA <- c(METADATA,"Strain")#Update METADATA array

# Get the species Name colum splitted as Genus + Species; Horrible onliner
cbind(setNames(
	       data.frame(do.call(rbind,strsplit(as.character(metadata$Name),split = " "))),
	       nm = c("Genus","Species")),metadata)->metadata

METADATA <- c(METADATA,"Genus","Species")#Update METADATA array

## Delete unneeded variables
metadata[,c("Type strain","DSM strain number","ATTC strain number","Other collections")]<- NULL

METADATA <- METADATA[-which(METADATA %in% c("Type strain","DSM strain number","ATTC strain number","Other collections"))]#UPDATE METADATA

#Sort manually the order of columns
metadata[,c("Genus", "Species", "Strain", "Name", "Main publication", "Other publication",
	    "Description", "Link", "Max. cell GC content")] -> metadata

metadata <- cbind(ID,metadata)# Link the ID for later uses

#save_file(metadata,"metadata.tsv")
#rm(metadata)


####################### DB 2: Numerical data related with Temperature and pH

numerical_phenotypes <- datos[,NUMERIC]

# Change the missing data to NA
numerical_phenotypes <- exchange_strings(numerical_phenotypes,"no data",NA)

#Make digits numbers again rather than characters, NA become NA values and non-digits would be coerced to NA's (latter should not happen now) 
#non_numerical <- check_non_numeric_value(numerical_phenotypes)
numerical_phenotypes <-apply(numerical_phenotypes,c(1,2),function(x){as.numeric(as.character(x))})
# Recover the data frame structure
numerical_phenotypes <- setNames(data.frame(numerical_phenotypes),nm=NUMERIC)

# for later merge/joins plug the ID column - not used at the moment
numerical_phenotypes <- cbind(ID,numerical_phenotypes)
#save_file(numerical_phenotypes,"numeric.tsv")
rm(numerical_phenotypes)

####### DB 3: Metabolic substrates for methanogenesis- including growth requisites

metabolic_phenotypes <- datos[,METABOLIC]

### Change symbols to NA, "yes", "no"
metabolic_phenotypes <- exchange_strings(metabolic_phenotypes,NA_STRINGS,NA)
metabolic_phenotypes <- exchange_strings(metabolic_phenotypes,YES_STRINGS,"yes")
metabolic_phenotypes <- exchange_strings(metabolic_phenotypes,NO_STRINGS,"no")

# Ensure  names (I think that is unneeded at the moment, due to a change in the function exchange_strings that does not check.names anymore
metabolic_phenotypes <- setNames(metabolic_phenotypes,nm = METABOLIC)

# Add now the ID to avoid problems with the cleaning.
metabolic_phenotypes <- cbind(ID,metabolic_phenotypes)

#non_numerical <- check_non_numeric_value(metabolic_phenotypes)


		## Fixing the columns with more than 1 value per field. ##

## First Column to manually split  - growth requirements

# Clean up poor curation of compound names

# A map-reduce would be good here, but I can´t see clear how to do it and I need to run faster; a for loop over a list of pairs would work
tmp <- levels(metabolic_phenotypes$'Additional growth requierments')

tmp <- exchange_strings_pattern(tmp," mono-, di-, a trimethylamine","monomethylamine,dimethylamine,trimethylamine")
tmp <- exchange_strings_pattern(tmp, "yea?st extract|easyt extract","yeast extract")
tmp <- exchange_strings_pattern(tmp,"\\(?acety?ate\\)?","acetate")
tmp <- exchange_strings_pattern(tmp,"Vitamins","vitamins")
tmp <- exchange_strings_pattern(tmp,"Selenium","selenium")
tmp <- exchange_strings_pattern(tmp,"rument","rumen")
tmp <- exchange_strings_pattern(tmp,"peptones","peptone")
tmp <- exchange_strings_pattern(tmp,"and",",")

# Apply the changes in the levels to the df
tmp -> levels(metabolic_phenotypes$'Additional growth requierments')

  ##-- Preparation of data.frame with the new colums --##

#NOTE: tmp was for levels - now is for the data.frame 

#Get the column and force it to be a character rather than factor - and set it as a temporal variable
tmp <- as.character(metabolic_phenotypes$'Additional growth requierments')
# Albeit temporal, keep the ID - and therefore tmp is a df with two colums ID and the phenotype
tmp <- cbind(ID,tmp)
tmp <- data.frame(tmp)

#Horrible one-liner to split the commas, line by line, into transposed vectors and add the ID (that is vicious in me but works)
tmp <- do.call(rbind,#-- The result of the next few lines is a list of dataframes with ID and pheno columns - just put them together by rows.

	       apply(tmp,1,function(x){#Apply to every row - tmp is now a df with column ID and colum Additional -
			     # Look in reverse: split the pheno list with commas, then unlist it to make a vector, column bind with its corresponding ID
			     # all this can be structure as a data.frame with two columns - and finally plug the ID and pheno to columnames
			     setNames(
				      data.frame(cbind(x[1],(unlist(strsplit(x[2],","))))),
				      nm = c("ID","pheno"))
	       }))

#Another idiomatic way sanitize whitespaces.
#It uses/abuses the mechanism for factors in R
levels(tmp$pheno)<-trimws(levels(tmp$pheno))#Trimming whitespaces on the new pheno column
levels(tmp$pheno)<-paste("Additional",levels(tmp$pheno))# Keep separated the Additional from the usual substrates

#Default value to show in the field for any of these
cbind(tmp,FLAG = "yes")->tmp

## the data frame is now ready to become a wide table with the same number of rows as the original one
reshape2::dcast(data =tmp, ID ~ pheno,value.var="FLAG")->tmp

# A rogue new variable has better sense with new nameing
names(tmp)[names(tmp) == "Additional no"] <- "Additional no supplement"

tmp$'NA'<-NULL# Those that had NA as value become a new variable
##---The data frame is now ready to be added to its parental data source

merge(metabolic_phenotypes,tmp,all.x = T)->metabolic_phenotypes
metabolic_phenotypes[,"Additional growth requierments"]<- NULL

#Update array with names 
METABOLIC <- c(METABOLIC,names(tmp)[!names(tmp) %in% "ID"])
METABOLIC[!(METABOLIC %in% "Additional growth requierments")]->METABOLIC

### Second manually splitted column - wont' go into details. Look the previous section for code explanation

#Clean up poor vocabulary 

tmp <- levels(metabolic_phenotypes[,"Min. growth requierments"])

# Commented line: have a look to the unique substrates before making the patterns
#sort(unique(trimws(unlist(strsplit(paste(levels(metabolic_phenotypes[,"Min. growth requierments"]),collapse = ","),",")))))[1:10]
tmp <- exchange_strings_pattern(tmp,"yest|yeast","yeast")
tmp <- exchange_strings_pattern(tmp,"(vitamin)* B12","vitamin B12")
tmp <- exchange_strings_pattern(tmp,"tungstate","tungsten")
tmp <- exchange_strings_pattern(tmp,"tr(i|Y)pticase( peptone){0,0}([,/]{1,})","tripticase peptone\\3")
tmp <- exchange_strings_pattern(tmp,"^tripticase$","tripticase peptone")
tmp <- exchange_strings_pattern(tmp,"rumen f(l*)uid","rumen fluid")
tmp <- exchange_strings_pattern(tmp,"or vitamins for growth",",vitamins")
tmp <- exchange_strings_pattern(tmp,"peptone(s*)","peptone")
tmp <- exchange_strings_pattern(tmp,"vitamin(s?)","vitamin")
tmp <- exchange_strings_pattern(tmp,"Co-M|coenzym M|coenzyme M","coenzyme M")
tmp <- exchange_strings_pattern(tmp,"casamino acids","casamino acids")
tmp <- exchange_strings_pattern(tmp,"and|or",",")
# apply changes to the metabolites
tmp -> levels(metabolic_phenotypes[,"Min. growth requierments"])

# Data frame preparatio ID + pheno
tmp <- as.character(metabolic_phenotypes[,"Min. growth requierments"])
cbind(ID,tmp)->tmp
tmp <- data.frame(tmp)

#one liner from multi-valued text to data.frame ID + pheno
tmp <- do.call(rbind,apply(tmp,1,function(x){setNames(data.frame(cbind(x[1],(unlist(strsplit(x[2],","))))),nm = c("ID","pheno"))}))

levels(tmp$pheno)<-trimws(levels(tmp$pheno))
levels(tmp$pheno)<-paste("Minimum",levels(tmp$pheno))

#Default value to show in the field for any of these
cbind(tmp,FLAG = "yes")->tmp
reshape2::dcast(data =tmp, ID ~ pheno,value.var="FLAG")->tmp
tmp$'NA'<- NULL
tmp$'Minimum '<- NULL
# Data frame is ready to be merged with the parental phenotypes:
merge(metabolic_phenotypes,tmp,by="ID")-> metabolic_phenotypes

#Deleted undesired columns
metabolic_phenotypes[,"Min. growth requierments"]<- NULL
metabolic_phenotypes[,"Minimum "]<- NULL

#Update array with names 
METABOLIC <- c(METABOLIC,names(tmp)[!names(tmp) %in% "ID"])
METABOLIC[!(METABOLIC %in% "Min. growth requierments")]->METABOLIC


# The third column that has to be splitted, other substrates, is more tricky
#since it has also a yes/no structure in the line(not supporting: in all line except 1).
# It has a compound separated by commas, the 2,3 butanediol

# 1st make a vector of yes/no/NA according to the not supporting (also written with mistakes

tmp2 <- metabolic_phenotypes[,"Other substrate"]
ifelse(grepl(x= levels(tmp2),"not suppo(r?)ting.*"),"no","yes")->levels(tmp2)
tmp2 <- cbind(ID = ID,FLAG = as.vector(tmp2))# Leave it aside while cooking the rest of it

# tmp2 has the yes/no flag since most of the column is not supporting except 131.


## Now, I have to clean up the errors in this one (that does not seem to be many)
tmp <- levels(metabolic_phenotypes[,"Other substrate"])

tmp <- exchange_strings_pattern(tmp,"cyclohexsanol","cyclohexanol")
tmp <- exchange_strings_pattern(tmp,"2,3-butanodiol","2;3-butanodiol")
tmp <- exchange_strings_pattern(tmp,"not suppo(r?)ting\\s*:","")
# Apply the changes to the original df
tmp -> levels(metabolic_phenotypes[,"Other substrate"])

# Now, lets prepare the data frame ID + pheno
tmp <- as.character(metabolic_phenotypes[,"Other substrate"])
cbind(ID,tmp)->tmp
tmp <- data.frame(tmp)
tmp <- do.call(rbind,apply(tmp,1,function(x){setNames(data.frame(cbind(x[1],(unlist(strsplit(x[2],","))))),nm = c("ID","pheno"))}))
# Clean up the new phenotypes
levels(tmp$pheno)<-trimws(levels(tmp$pheno))
#Distinguish these substrates from the other ones (it would be good to check if they are already on the list
levels(tmp$pheno)<-paste("Other",levels(tmp$pheno))

reshape2::dcast(merge(tmp,data.frame(tmp2),all = T),formula = ID ~ pheno,value.var ="FLAG")->tmp
tmp$'NA' <- NULL
tmp$'Other ' <- NULL
merge(metabolic_phenotypes,tmp,by="ID")->metabolic_phenotypes

metabolic_phenotypes[,"Other substrate"]<- NULL
metabolic_phenotypes[,"NA.x"]<- NULL
metabolic_phenotypes[,"NA.y"]<- NULL
metabolic_phenotypes[,"NA"]<- NULL
metabolic_phenotypes$'Other '<-NULL

#Update array with names 
METABOLIC <- c(METABOLIC,names(tmp)[!names(tmp) %in% "ID"])
METABOLIC[!(METABOLIC %in% "Other substrate")]->METABOLIC

# It seems that most of the phenotypes are now empty!!!!
#> str(data.frame(do.call(cbind,lapply(names(metabolic_phenotypes),function(x,df){if(is.character(df[[x]])){as.factor(df[[x]])}else{df[[x]]}},df = metabolic_phenotypes)),check.names=F))

#save_file(metabolic_phenotypes,"metabolic.tsv")
#rm(metabolic_phenotypes)


########## Environment related phenotypes #######
# Colums as "live in such environment: Intestinal.tracks" and FLAGGED as yes/no
# Irrelevant for genephene project
environmental_phenotypes <- datos[,ENVIRONMENTS]
environmental_phenotypes <- exchange_strings(environmental_phenotypes,c("no data","not applicable"),NA)
environmental_phenotypes <- exchange_strings(environmental_phenotypes,c("yes","1"),"yes")
#non_numerical <- check_non_numeric_value(environmental_phenotypes)
#save_file(environmental_phenotypes,"environmental.tsv")
rm(environmental_phenotypes)
## WRAP UP - calculate what I wanted 

# Finally, I can make the database I wanted with just the Substrates (There were not metabolic process in this DB)
merge(metadata[,c("ID","Name","Genus","Species","Strain")],metabolic_phenotypes,by="ID")->Substrates_DB
Substrates_DB$ID <- NULL
save_file(Substrates_DB,"Substrates.tsv")
# This work should be done in bash for the sake of clarity, but it is easy to do it here
cat(METABOLIC,sep ="\n",file = "list_of_phenotypes.txt")
write.table(metadata[,c("Name","Strain")],"list_of_species.txt",quote = T,row.names=F,sep ="\t")
#summary of phenotype coverage
apply(metabolic_phenotypes[METABOLIC],2,function(x){paste("Coverage: ",round(sum(!is.na(x))/length(x)*100,digits =2),"%")})->tmp
write.table(x = data.frame(tmp),"summary_phenotypes_coverage.txt",quote = F,row.names=T,sep = "\t")
apply(metabolic_phenotypes[METABOLIC],2,function(x){paste("Balance :",round(sum(x=="yes",na.rm=T)/sum(!is.na(x),na.rm=T)*100,digits =2),"%")})->tmp
write.table(x = data.frame(tmp),"summary_phenotypes_balance.txt",quote = F,row.names=T,sep ="\t")
