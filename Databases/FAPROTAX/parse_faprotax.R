### Preamble ####
#FAPROTAX database parsing based on regexp per line

## TODO: correct species not included

############################################### PARSING THE DATABASE ONTO A LIST STRUCTURE #######################
# Variables and counters 
DB <- list()

line <- 0
pheno_counter <- 0
counter <- 0


### Set up of files ###
file("FAPROTAX.txt",open="rt")->FHdata


# Init reading the file. as a do-while
seek(FHdata,where = 0)
readLines(FHdata,n=1)->line

# carry on reading the file
while(length(line)>0){
	counter <- counter + 1
	#print(line)
	#print(counter)
        ## Four types of lines: empty, those starting by #, Those starting by other symbol than # or *, and those starting by *
        #corresponcdence # -> comments, * -> species | [^#*]->phenotype. But the latter has an issue, there are lines to extend or contract the phenotypes
        # The logic seems to be that they have use set operators to simplify the database (simplify as in make it short
        # Then add_group is a 'contains' operation. substract_group seems to be set difference and instersect_group 
        # gladly it seems that there are not substract or intersect
        ## The lines starting with a group (followed by #----) define a major phenotype group plus a set of metadata entries that are phenotypes with functional dependencies of the major group in the databases sense. That means, the same taxa has both phenotypes, likely to be because they are redundant.

            ## Comments and empty lines are not needed
	#if(grepl(x = line,pattern="^$",ignore.case=T,perl=T)){print("empty line")}
	#if(grepl(x = line,pattern="^\\#",ignore.case=T,perl=T)){#print("comment found")}

                    #set operation lines
	if(grepl(x = line,pattern="^[^\\#\\*]",ignore.case=T,perl=T)){
		if(grepl(x = line,pattern="^add_|^subtract_|^instersect_",ignore.case=T,perl=T)){
			#print("A group here")
			DB[[pheno_counter]]$groups <- c(DB[[pheno_counter]]$groups,line)
		
		}else if(!grepl(x =line,pattern="\\*")){            #Phenotype group line- there are errors where the species do not have an * to start

			#print("phenotype found")
			pheno_counter <- pheno_counter + 1

			#declare the new variables
			DB[[pheno_counter]] <- list()
                        DB[[pheno_counter]]$species <- c()
			DB[[pheno_counter]]$groups <- c()

                        # assign the current phenotypic line  
			DB[[pheno_counter]]$phenotype <- line
		}
	}else if(grepl(x = line,pattern="^\\*",ignore.case=T,perl=T)){
		#print("species found")
		DB[[pheno_counter]]$species <- c(DB[[pheno_counter]]$species,line)
	}
	
        #if(counter >= 200){break}
        
        ## #Read a new line###
	readLines(FHdata,n=1)->line
}

#seek(FHdata,where = 0)

close(FHdata)
rm(list = ls()[!ls() %in% "DB"])

############################## PROCESSING OF THE DATA STRUCTURE #############################
### parse the strings 


# manual mistakes in the DB - Species that are not spelt right
## code to help correction during the curation  grep(x =lapply(processed_DB,function(ll){ll$species$unlikelySpecies})[48][[1]],pattern ="(ceae)|(ales)|^(\\w+)$",value =T,inv=T)[1:40]
mistaken  <- c("Pseudomonas Oleovorans","Methylomonas Clara","Bacteroidetes Gracilimonas","Thauera Azoarcus")
mispelt   <- c("treptococcus difficile","higella sonnei")
strange   <- c("Desulfosporosinus OT","Desulfitobacterium GBFH","Bacillus SXB","Pantoea IMH","Clostridium OhilAs","Clostridium 9.1","Thiomonas 3As",
    "Thiomonas WZW","SUP05","Arctic96BD-19","Thiomonas SS","Geobacteraceae Dfrl","SAR86 clade","SAR11 clade","Salmonella Typhimurium SL 1344","Marinobacter CAB","Rhizobium BTAil","Gordonia MTCC 4818") 
#There's at least 1 extra mispelt "Cor ynebacterium freneyi" that I solve at the very end of the script to do not make it more complex.


processed_DB <- list()
elem_n <- 0

for (element in DB){

    elem_n <- elem_n +1 
     processed_DB[[elem_n]]<-list()


    #process the phenotype name - that is actually a group definition according to FAPROTAX file structure and their readme file
    if(exists("phenotype",element)){
        processed_DB[[elem_n]]$phenotype <- list()

        # The list phenotypes tend to be <group_name>\t<list of redundant groups in a format class1:a1,a2;class2:b1,b2;... >    
        # or <group_name>\t# <bibliography>

        strsplit(element$phenotype,"\t+")->tokenized
        unlist(tokenized)->tokenized
        #Note: phenotypic names may contain spaces that are not in later add_groups operations - therefore trimwhite space
        processed_DB[[elem_n]]$phenotype$main <- trimws(tokenized[[1]])

        #parse the line with extra information, although I won't use it.
        if(length(tokenized)>1){
            tokenized[2]->tokenized
            #separate the groups <class>:type1,type2
            trimws(unlist(strsplit(tokenized,";")))->tokenized #Must be now an array class:t1,t2...
            lapply(strsplit(tokenized,":"),function(tok){paste(tok[1],trimws(unlist(strsplit(tok[-1],","))),sep=":")})->tokenized
            unlist(tokenized)->tokenized    

            processed_DB[[elem_n]]$phenotype$subcategories <- tokenized
        }
    }
    
    if(exists("species",element)){

        processed_DB[[elem_n]]$species <- list()
            # Each species line starts with the species name without spaces but with * before and after every word *Escherichia*coli*
            # then several tabs and a comment, usually with literature \t\t\t# blah blah
            # cut the  line i two, by splitting with tabs
        unlist(lapply(element$species,function(x){strsplit(x,"\t+",perl=T,useBytes=T)[[1]][1]->y
        #Remove the *, and sanitize
        trimws(gsub(x = y,pattern="\\*+",replacement=" "))}))->possibleSpecies
        # The DB is full of name pairs that does not correspond to species but to a kind of phylogenetic graph e.g. Gammaproteobacteria Methylococcales
        # Although the DB contain errors - see below - Assuming that Genus name starts in capital letters and the specific name in non capital make a nice distinction 
        grep(x = possibleSpecies,pattern= "[A-Z]+\\D+\\s+[a-z]+",value = T,perl=T,ignore.case = F,invert=F)->likelySpecies
        grep(x = possibleSpecies,pattern= "[A-Z]+\\D+\\s+[a-z]+",value = T,perl=T,ignore.case = F,invert=T)->unlikelySpecies
        #The last two sentences -repeated except the invert flag - generates two vectors that I will retain
        
        # Yet, I found many bacterias that has a Genus and then a strain type-like string. The next regexp find those that the specific name start with words and
        # carry on with either words or numbers. Then, remove those from the unlikelySpecies vector.

        likelySpecies<-c(likelySpecies,grep(x = unlikelySpecies,pattern = "^[A-Z]+\\w+\\s+[\\w]+[-\\d+]",ignore.case=F,perl=T,value = T))
        unlikelySpecies <- setdiff(unlikelySpecies,likelySpecies)
        
        # On the other hand there were strings in likelySpecies with the Taxa1 taxa2 format, in particular the mistake is having genus names with no capital.
        # Those has to be moved from likelySpecies to unlikelySpecies
        unlikelySpecies <- c(unlikelySpecies,grep(x = likelySpecies,pattern = "((ceae)|(ales))\\s+(([a-z]\\w+)|(\\D+-\\D+))",ignore.case=F,perl=T,value = T))
        likelySpecies <- setdiff(likelySpecies,unlikelySpecies)

        # Yet, there are some other mistakes like strain-types without numbers or all numbers - Species name with specific epithet in capital letters ...
        # missing letters ...
        likelySpecies<-c(likelySpecies,intersect(unlikelySpecies,c(mistaken,strange)))
        unlikelySpecies <- setdiff(unlikelySpecies,likelySpecies)#probably inefficient - but it does not matter at the moment.
        
        #Two species have a missing letter at the beggining - So I take advantage of being the same letter (if another bacteria happen to be mispelt, I need to change it)
        likelySpecies<-c(likelySpecies,if(length(Sps<-intersect(unlikelySpecies,mispelt))>0){paste('S',Sps,sep="")})
        unlikelySpecies <- setdiff(unlikelySpecies,likelySpecies)#probably inefficient - but it does not matter at the moment.
        
        
        likelySpecies ->processed_DB[[elem_n]]$species$likelySpecies
        unlikelySpecies ->processed_DB[[elem_n]]$species$unlikelySpecies
        
        rm(unlikelySpecies,likelySpecies,possibleSpecies,Sps)
    }

    if(exists("groups",element)){
	#parse the operations
		#Sanitize
	gsub(x = element$groups,pattern="\t*#.*",replacement = "")->groups
	read.table(text = groups,sep = ":",header =F,col.names = c("operation","out_group"))->groups

        processed_DB[[elem_n]]$group_set_op <- groups 
        rm(groups)
    }


    rm(tokenized)
}

rm(elem_n)
rm(element)
rm(mistaken,mispelt,strange)
######################################### GROUPS OPERATIONS ####################################################################
### Force the group operation now that the species are curated and all strings are supposedly parsed.

## Some ad-hoc helper functions

process_group_operations <-function(DB,element){

    search_element <-function(text){
        #by not testing if it exists or not, I am catching errors that sould not be happening, like whitespaces
        DB[[which(sapply(DB,function(x,y){x$phenotype$main == y},y = text,simplify=T,USE.NAMES=F))]]
    }

    for(row in 1:nrow(element[["group_set_op"]])){
        if (element[["group_set_op"]][row,"operation"] == "add_group"){
            included_element <- search_element(element[["group_set_op"]][row,"out_group"])
            element$species$likelySpecies <- c(element$species$likelySpecies,included_element$species$likelySpecies)
        }else if (element[["group_set_op"]][row,"operation"] == "subtract_group"){
            included_element <- search_element(element[["group_set_op"]][row,"out_group"])
            element$species$likelySpecies <- setdiff(element$species$likelySpecies,included_element$species$likelySpecies)
        }else if (element[["group_set_op"]][row,"operation"] == "intersect_group"){
            included_element <- search_element(element[["group_set_op"]][row,"out_group"])
            element$species$likelySpecies <- interesect(element$species$likelySpecies,included_element$species$likelySpecies)
        }#There should not exist neither intersect nor any other set operation
    }
    #return
    element
}


for(element_n in 1:length(processed_DB)){
    element <- processed_DB[[element_n]]

    if(exists("group_set_op",element)){
          processed_DB[[element_n]]<- process_group_operations(processed_DB,element)
    }
}

rm(element,element_n)
################################## MAKE SPECIES TABLE! ##########################

final_DB <- list()
for(element in processed_DB){
    if(length(element$species$likelySpecies)>0){
        final_DB <- c(final_DB,list(as.data.frame(cbind(element$phenotype$main,element$species$likelySpecies))))
        }
}
do.call(rbind,final_DB)->final_DB
unique(final_DB)->final_DB
names(final_DB) <- c("Trait","Name")
reshape2::dcast(cbind(final_DB,value = "yes"),Name~Trait)->DB
################################ REPARE THE NAMES ######################

tmp<-do.call(rbind,lapply(strsplit(sub(x=as.vector(DB[["Name"]]),pattern =" ",replacement=";"),";"),function(x){if(length(x)<2){x <- c("",x)};setNames(c(x),c("Genus","Species"))}))
tmp <- data.frame(tmp)
#tmp <- cbind(tmp,Strain = gsub(x = tmp$Species,pattern = "(?!(w*\\s*)[-\\w]*\\d+\\w*)^.*",replacement = "*",perl=T,ignore.case=T,useBytes=T))
#tmp$Species <- gsub(x = tmp$Species,pattern = "(?:(\\w*\\s*)([-\\w]*\\d+\\w*))",replacement = "*",perl=T,ignore.case=T,useBytes=T)

tmp <-cbind(tmp,Strain = gsub(x = tmp$Species,pattern = "([-\\w]*\\d+\\w*)|(\\b[A-Z]+\\b)|(\\w+-\\w+)|(\\w+)",replacement = "\\1\\2\\3",perl=T,ignore.case=F,useBytes=T))
tmp$Species <- cbind(apply(tmp,1,function(x){ifelse(x["Strain"]!="",gsub(x = x["Species"],pattern = x["Strain"],replacement="",useBytes=T,ignore.case=F,fixed = T),x["Species"])}))

DB <- cbind(tmp,DB)

################ A last minute manual species corrections 

levels(DB$Name)[grep(x = levels(DB$Name),pattern="Cor ynebacterium freneyi")] <- "Corynebacterium freneyi"

################################ Finally Save it!!! ####################
write.table(x = DB,file = "FAPROTAX.tsv",sep = "\t",quote = F,row.names = F)
