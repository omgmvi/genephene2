# Count the level of sparsity
read.csv('Metabolic_features.csv',check.names = F)->datos
cat(paste(names(datos)[-1],apply(datos[-1],1,function(x){round(length(na.omit(x))/length(x)*100,digits = 2)}),sep = " Coverage:=",collapse = "\n"),'\n')
read.csv('Metabolic_traits.csv',check.names = F)->datos
cat(paste(names(datos)[-1],apply(datos[-1],1,function(x){round(length(na.omit(x))/length(x)*100,digits = 2)}),sep = " Coverage:=",collapse = "\n"),'\n')
