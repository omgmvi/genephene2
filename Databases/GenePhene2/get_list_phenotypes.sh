#!/bin/bash
head ./Metabolic_features.tsv -n 1 |cut -d $'\t' -f1-4 --complement|tr '\t' '\n' > list_of_phenotypes.txt
head ./Metabolic_traits.tsv -n 1|cut -d $'\t' -f1-4 --complement|tr '\t' '\n' >> list_of_phenotypes.txt
