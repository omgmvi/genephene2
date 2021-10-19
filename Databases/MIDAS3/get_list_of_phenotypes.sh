head MiDAS_Metadata.csv -n 1|tr ";" "\n" |tr -d  '\"' > list_of_phenotypes.txt
