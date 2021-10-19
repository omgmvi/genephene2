#!/bin/bash

python3 see_csv.py > db_modified.csv

column -s ";" -t < db_modified.csv |less -#2 -N -S
rm db_modified.csv
