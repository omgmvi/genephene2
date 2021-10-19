#!/bin/bash
#head Substrates.csv -n 1|cut -d ',' -f 1-4 --complement|tr "," "\n"
#head Metabolic.csv -n 1|cut -d ',' -f 1-4 --complement|tr "," "\n"
{ (head ./Metabolic.tsv -n 1|cut -f 1-4 --complement|tr "\t" "\n"|tr -d '\"' ) && (head ./Substrates.tsv -n 1|cut  -f 1-4 --complement|tr '\t' '\n'|tr -d '\"') }>list_of_phenotypes.txt
#|nl

