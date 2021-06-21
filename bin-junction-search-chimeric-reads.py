#! /usr/bin/env python

import os, sys, math
from optparse import OptionParser
from collections import defaultdict
import itertools

# This script searches for sequencing reads that two ends mapped far apart on the genome. 
# Output: A list of paired-regions where the chimeric reads mapped, then followed with the number of chimeric reads in each sample. 
# Output format: 
# Ref1 Bin-Start1 Ref2 Bin-Start2 Sample1 Sample2 Sample3 ... Samplen


usage = ""
parser = OptionParser(usage=usage)
#parser.add_option("-f", "--file", dest="f", help="Input .sam file.")
parser.add_option("-o", "--output_file", dest="o", help="Output file name.")
parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=5000, help="Bin size (5000)")
parser.add_option("-m", "--mapqual", dest="mapqual", type = "int", default=10, help="Minimum mappign quality (10)")
parser.add_option("-i", "--insert", dest="insert", type = "int", default=2000, help="maximum estimated insert between forward and reverse (2000)")
parser.add_option("-c", "--mincov", dest="mincov", type = "int", default=0, help="Minimum coverage for a bin to be in the output, across all libraries")


(opt, args) = parser.parse_args()


mindiff = opt.insert
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
      tup.append([x[2], int(x[3]), int(x[4])])
      templines = []
      while 1:
         temp = f.tell()
         next = f.readline()
         xn = next.split('\t')
         if xn[0] == x[0]:
            tup.append([xn[2], int(xn[3]), int(xn[4])])
            templines.append(xn)
         else:
            f.seek(temp)
            break
      combos = []      
      for each in itertools.combinations(tup,2):
         if "Chr" not in each[0][0] or "Chr" not in each[1][0]:
            continue
         if each[0][0] == each[1][0] and abs(each[0][1] - each[1][1]) < mindiff:
            continue
         if each[0][2] < minmap or each[1][2] < minmap:
            continue
         combos.append(each)
   #   flag = 0
   #   if len(combos) > 1:
   #      flag = len(combos)
      for item in combos:
         matches.append(list(item))
   
   f.close()
   
   for item in matches:
      item.sort(key = lambda x: x[0])
      chr1 = item[0][0]
      pos1 = item[0][1]
      key1 = str(pos1 / binsize)
      chr2 = item[1][0]
      pos2 = item[1][1]
      key2 = str(pos2 / binsize)
      bigkey = '-'.join([chr1,key1,chr2,key2])
      bigkey2 = '-'.join([chr2,key2,chr1,key1])
      bins[bigkey][fname] +=1
      bins[bigkey2][fname] +=1   



vbins = bins.keys()
vbins.sort(key = lambda x: int(x.split('-')[3]))
vbins.sort(key = lambda x: x.split('-')[2])
vbins.sort(key = lambda x: int(x.split('-')[1]))
vbins.sort(key = lambda x: x.split('-')[0])

head = ['Ref1', 'Ref1-BinStart', 'Ref2', 'Ref2-BinStart']
for fn in fileset:
   head.append(fn.split('_aln')[0])

o.write('\t'.join(head)+'\n')

for item in vbins:
   tup = item.split('-')
   #NOTE THIS
#   if tup[0] == tup[2] and tup[1] == tup[3]:
#      continue
   line = [tup[0], int(tup[1])*binsize, tup[2], int(tup[3])*binsize]
   calls = []
   for file in fileset:
      temp = bins[item][file]
      line.append(temp)
      calls.append(temp)
   if sum(calls) >= opt.mincov:
      oline = map(lambda x: str(x), line)
      o.write('\t'.join(oline)+'\n')


o.close()




   
