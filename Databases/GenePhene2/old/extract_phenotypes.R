read.csv('Metabolic_features.csv',check.names = F)->datos
cat(paste(names(datos)[-1],collapse = "\n"),'\n')
read.csv('Metabolic_traits.csv',check.names = F)->datos
cat(paste(names(datos)[-1],collapse = "\n"),'\n')
