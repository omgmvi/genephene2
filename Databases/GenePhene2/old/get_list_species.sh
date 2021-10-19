cat Metabolic_features.csv  |cut -d , -f 1|tail -n +2 > list_of_species.txt
cat Metabolic_traits.csv  |cut -d , -f 1|tail -n +2 >> list_of_species.txt
sort list_of_species.txt|uniq > list_of_species1
mv list_of_species1 list_of_species.txt
