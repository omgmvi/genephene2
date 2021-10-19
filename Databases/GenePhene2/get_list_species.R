read.table("./Metabolic_features.tsv",header = T,sep = "\t") ->Mf
Mf[,c("Name","Strain")]->Mf

read.table("./Metabolic_traits.tsv",header = T, sep="\t")->Mt
Mt[,c("Name","Strain")]->Mt

merge(Mf,Mt,all=T)->All

write.table(All,"list_of_species.txt",sep = "\t",row.names = F)
