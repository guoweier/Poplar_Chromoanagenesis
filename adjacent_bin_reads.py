#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 12:53:02 2018

@author: weier
"""

import os, sys, math
from optparse import OptionParser
from collections import defaultdict

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-o", "--output_file", dest="o", help="Output file name.")
parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=5000, help="Bin size (5000)")
parser.add_option("-m", "--mapqual", dest="mapqual", type = "int", default=10, help="Minimum mappign quality (10)")
parser.add_option("-c", "--mincov", dest="mincov", type = "int", default=0, help="Minimum coverage for a bin to be in the output, across all libraries")

(opt, args) = parser.parse_args()

minmap = opt.mapqual
binsize = opt.binsize
(opt, args) = parser.parse_args()

chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))
o = open(opt.o,'w')

#read in list of sam files, must be in current directory and end in "_aln.sam"
li = os.listdir(os.getcwd())
fileset = filter(lambda x: x.endswith('.sam'), li)
fileset.sort()

for fname in fileset:
   print fname
   matches = []
   f = open(fname)
   while 1:
       l = f.readline()
       if l == '':
           break
       if l[0] == '@' and l[1] != "P":
           h = l.split()
           name = h[1].replace("SN:","")
           nlen = int(h[2].replace("LN:",""))
           chromlist[name] = nlen
           continue
       if l[0] == "@" and l[1] == "P":
           continue
       x = l.split('\t')
       tup = []
       tup.append([x[2], int(x[3]), int(x[4]), x[9]])
       if x[2][0] != "C":
           continue
       else:
           if int(x[4]) < minmap:
               continue
           else:
               if (int(x[3]) // binsize) != (int(x[3]) + len(x[9])) // binsize:
                   name = x[2] + "_" + str((int(x[3]) // binsize + 1) * binsize)
                   bins[name][fname] += 1
               else:
                   continue
   f.close()
    
vbins = bins.keys()
vbins.sort(key = lambda x: int(x.split('_')[1]))
vbins.sort(key = lambda x: x.split('_')[0])

head = ['Ref', 'Bin_Position']
for fn in fileset:
   head.append(fn.split('_aln')[0])

o.write('\t'.join(head)+'\n')

for item in vbins:
   tup = item.split('_')
   line = [tup[0], int(tup[1])]
   calls = []
   for file in fileset:
      temp = bins[item][file]
      line.append(temp)
      calls.append(temp)
   if sum(calls) >= opt.mincov:
      oline = map(lambda x: str(x), line)
      o.write('\t'.join(oline)+'\n')


o.close()
    
             
     
     
     
     
     
     