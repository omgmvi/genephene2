#!/bin/bash
cat FAPROTAX.tsv |column -s $'\t' -t |less -#2 -S -N
