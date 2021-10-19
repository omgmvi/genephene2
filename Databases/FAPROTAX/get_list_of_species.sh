#!/bin/bash

#cat FAPROTAX.tsv |cut -d $'\t' -f 1|tail -n +2|awk '{ print "\""$0"\""}' > list_of_species.txt


echo -e '"Name"\t"Strain"' > list_of_species.txt

#paste <(cat FAPROTAX.tsv |cut -d $'\t' -f 1,2|tail -n +2|awk '{ print "\""$0"\""}'|tr '\t' ' ') <(cat FAPROTAX.tsv |cut -d $'\t' -f 3|tail -n+2| awk '{ print "\""$0"\""}') >>list_of_species.txt

paste <(cat FAPROTAX.tsv |cut -d $'\t' -f 1,2|tail -n +2|awk '{ print "\""$0"\""}'|tr '\t' ' '|sed -r -e "s/ *\"$/\"/g") <(cat FAPROTAX.tsv |cut -d $'\t' -f 3|tail -n+2| awk '{ print "\""$0"\""}') >> list_of_species.txt

# This destroy the other parts of the script, specifically the SpeciesName2TaxID.R
#paste <(cat FAPROTAX.tsv |cut -d $'\t' -f 4|tail -n +2|awk '{ print "\""$0"\""}') <(cat FAPROTAX.tsv |cut -d $'\t' -f 3|tail -n +2| awk '{ print "\""$0"\""}') >> list_of_species.txt

