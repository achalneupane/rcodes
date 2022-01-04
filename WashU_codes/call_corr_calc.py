#!/usr/bin/env python
# coding: utf-8

import sys
from csv import reader

with open(sys.argv[1], 'w') as f:
    csv_reader = reader(sys.stdin,delimiter='\t')
    prevLine = ""
    for row in csv_reader:
        row_edit = row[9:]
        aa = sum([i == j for i, j in zip(row_edit, prevLine)])
        corr = (aa*100)/len(row_edit)
        prevLine = row[9:]
        #print(row[2],corr)
        f.write(row[2]+'\t'+str(corr)+'\n')
        

