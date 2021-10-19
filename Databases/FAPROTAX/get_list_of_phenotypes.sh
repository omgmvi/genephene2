#!/bin/bash
head FAPROTAX.tsv -n 1|tr '\t' '\n' |tail -n +2 > list_of_phenotypes.txt
