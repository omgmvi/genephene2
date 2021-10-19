#!/bin/bash
cat MiDAS_Metadata.csv |column -s ";" -t |less -#2 -N -S
