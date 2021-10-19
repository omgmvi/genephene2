data <- read.table("ijsem.tsv",header = T,sep = "\t",comment.char = "")[c("Genus","Species","Strain")]
data <- (as.data.frame(apply(data,2,function(x){as.character(x)}),stringsAsFactors = F))
data$Name <- paste(data[,1],data[,2])
data[c("Genus","Species")]<-NULL
write.table(data[c("Name","Strain")],"list_of_species.txt",row.names = F,sep = "\t",quote = T)
