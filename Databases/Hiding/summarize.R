

read.table(file = "Metabolic.tsv",header = T,sep = "\t",comment.char="")-> MetabolicTraits
summary_metabolic <- apply(MetabolicTraits[,-c(1:3)],2,function(x){paste("Coverage: ",round(sum(x == "yes",na.rm=T)/NROW(x)*100,digits = 2))})
write.table(summary_metabolic, file = "summary_metabolic.tsv",row.names = T,sep = "\t",quote = F,col.names = F) #The percentage of bacterias phenotyped for a given metabolism


read.table(file = "Substrates.tsv",header = T,sep = "\t",comment.char="",check.names = F)-> CarbonSubstrates
summary_carbonsubstrates <- apply(CarbonSubstrates[,-c(1:3)],2,function(x){paste("Coverage: ",round(sum(x == "yes",na.rm=T)/NROW(x)*100,digits = 2))})
write.table(summary_carbonsubstrates, file = "summary_carbonsource.tsv",row.names = T,sep = "\t",quot = F,col.names = F)#same for carbon substrates
