#! /usr/bin/env python

import os, sys, math
from optparse import OptionParser
from collections import defaultdict
import itertools

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", dest="f", help="Input .sam file.")
parser.add_option("-o", "--output_file", dest="o", help="Output file name.")
parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=5000, help="Bin size (5000)")
parser.add_option("-r", "--chromosome", dest="chromosome", help = "Chromosome for junction")
parser.add_option("-a", "--leftstart", dest="leftstart", type = "int", help="Left start position for junction")
parser.add_option("-b", "--leftend", dest="leftend", type = "int", help="Left end position for junction")
parser.add_option("-y", "--rightstart", dest="rightstart", type = "int", help="Right start position for junction")
parser.add_option("-z", "--rightend", dest="rightend", type = "int", help="Right end position for junction")
parser.add_option("-m", "--mapqual", dest="mapqual", type = "int", default=10, help="Minimum mappign quality (10)")
parser.add_option("-i", "--insert", dest="insert", type = "int", default=2000, help="maximum estimated insert between forward and reverse (2000)")
parser.add_option("-c", "--mincov", dest="mincov", type = "int", default=0, help="Minimum coverage for a bin to be in the output, across all libraries")


(opt, args) = parser.parse_args()

#commands setting
fname = opt.f
mindiff = opt.insert
minmap = opt.mapqual
binsize = opt.binsize
chrom = opt.chromosome
leftstart = opt.leftstart
leftend = opt.leftend
rightstart = opt.rightstart
rightend = opt.rightend
(opt, args) = parser.parse_args()

chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))
o = open(opt.o,'w')



matches = []
f = open(fname)
while 1:
    l = f.readline()
    if l == '':
        break
    #edit data that is not reads sequences
    if l[0] == '@':
        h = l.split()
        name = h[1].replace("SN:","")
        nlen = int(h[2].replace("LN:",""))
        chromlist[name] = nlen
        continue
    x = l.split('\t')
    tup = []
    tup.append([x[0], x[2], int(x[3]), int(x[4]), x[9]])
    templines = []
    while 1:
        temp = f.tell()
        next = f.readline()
        xn = next.split('\t')
        if xn[0] == x[0]:
            tup.append([xn[0], xn[2], int(xn[3]), int(xn[4]), xn[9]])
            templines.append(xn)
        else:
            f.seek(temp)
            break
    combos = []   

    for each in itertools.combinations(tup,2):
        if "Chr" not in each[0][1] or "Chr" not in each[1][1]:
            continue
        if each[0][1] == each[1][1] and abs(each[0][2] - each[1][2]) < mindiff:
            continue
        if each[0][3] < minmap or each[1][3] < minmap:
            continue
        if each[0][1] == each[1][1] == chrom:
            if (each[0][2] > leftstart and each[0][2] < leftend) and (each[1][2] > rightstart and each[1][2] < rightend):
                 combos.append(each)
            elif (each[1][2] > leftstart and each[1][2] < leftend) and (each[0][2] > rightstart and each[0][2] < rightend):
                combos.append(each)
            else: 
                 continue

   #   flag = 0
   #   if len(combos) > 1:
   #      flag = len(combos)
    for item in combos:
        matches.append(list(item))
      
      
f.close()

#put selected reads name and sequences into two tables
list_name = []
list_seq = []      
for item in matches:
    list_name.append(item[0][0])
    list_seq.append(item[0][4])
   
#write tables into fasta form
for i in range(len(list_seq)):
    o.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")
  
#don't forget to close the write file
o.close()





   