#cat <(cat Metabolic_features.csv |head -n 1|tr ',' '\n')|tail -n +2 > list_of_phenotypes.txt
#cat <(cat Metabolic_traits.csv |head -n 1|tr ',' '\n') |tail -n +2>> list_of_phenotypes.txt
Rscript extract_phenotypes.R > list_of_phenotypes.txt
