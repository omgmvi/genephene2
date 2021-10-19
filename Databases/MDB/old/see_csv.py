#!/usr/bin/env python3

import csv
with open('MDB.csv', newline='') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for row in linereader:
         print('; '.join(row))

