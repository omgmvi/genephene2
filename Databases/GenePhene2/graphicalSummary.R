

read.table("Metabolic_traits.tsv",sep = "\t",header =T)->datos
datos[-which(names(datos)%in% c("Genus","Species","Strain","Name"))]->datos


sort(apply(datos,2,function(x){100*(sum(!is.na(x))/length(x))}),decreasing=T)->coverage
sort(apply(datos,2,function(x){sum(x == "yes",na.rm = T)/sum(!is.na(x))}),decreasing=T)->balance

print("Metabolic traits")
txtplot::txtplot(y =coverage,x = 1:length(coverage))
txtplot::txtplot(y =balance,x = 1:length(balance))

invisible(readline(prompt="Press [enter] to continue"))

read.table("Metabolic_features.tsv",sep = "\t",header =T)->datos
datos[-which(names(datos)%in% c("Genus","Species","Strain","Name"))]->datos


sort(apply(datos,2,function(x){100*(sum(!is.na(x))/length(x))}),decreasing=T)->coverage
sort(apply(datos,2,function(x){sum(x == "yes",na.rm = T)/sum(!is.na(x))}),decreasing=T)->balance

print("Metabolic features")
txtplot::txtplot(y =coverage,x = 1:length(coverage))
txtplot::txtplot(y =balance,x = 1:length(balance))
