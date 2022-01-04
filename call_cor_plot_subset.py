#!/usr/bin/env python
# coding: utf-8
import sys
from csv import reader

with open(sys.argv[1], 'w') as f:
    csv_reader = reader(sys.stdin,delimiter='\t')
    prevLine = ("0",0,0,0,"A","G",0,0,0,0,0,0)   
    for row in csv_reader:
        if (row [6]== prevLine[6] and row [10]== prevLine[10]):
            f.write('\t'.join([str(n) for n in row])+'\n')
            prevLine = row
        else:
            prevLine = row
            next

