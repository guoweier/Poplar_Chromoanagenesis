#! /usr/bin/env python

import os, sys, math
from optparse import OptionParser
from collections import defaultdict
import itertools

#Introduction:
#This is the script generating seeds for PRICE assembly. 
#After selecting potential junctions through dosage plots and chimeric reads, PRICE assembly is used to expanding the contigs on the junction regions before PCR. 
#For PRICE assembly, seeds are reads that exactly on the junctions; surrounding reads are reads that near the junctions but not exactly on the junctions.
#For seeds, these reads are exactly cross the junctions, so two ends of these reads are far apart. So they are actually the chimeric reads counted before. 

#Input file: 
#Potential junctions in chimeric reads format with title line. i.e. Ref1 Ref1-BinStart Ref2 Ref2-BinStart sample1 sample2 ... samplen 
#Each element is separated by tab. The file is txt. 

#Output file:
#fasta file
#1st line: potential junction; list all seeds reads within fasta file format (name + sequence) 


usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f", "--input_file", dest="f", help="Txt file for selecting chimeric reads.")
#parser.add_option("-n", "--name", dest="name", help="File with list of sample names")
parser.add_option("-b", "--bin_size", dest="b", type = "int", default=500, help="Bin size used for chimeric reads search (500)")
parser.add_option("-m", "--mapqual", dest="mapqual", type = "int", default=10, help="Minimum mappign quality (10)")
parser.add_option("-i", "--insert", dest="insert", type = "int", default=2000, help="maximum estimated insert between forward and reverse (2000)")
parser.add_option("-c", "--mincov", dest="mincov", type = "int", default=0, help="Minimum coverage for a bin to be in the output, across all libraries")


(opt, args) = parser.parse_args()


txtfile = opt.f
binsize = opt.b
mindiff = opt.insert
minmap = opt.mapqual
(opt, args) = parser.parse_args()

chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))

#set output dictories and tables
list_name = {}
list_seq = {}
sampname_chi = []

#open txt file
txt = open(txtfile,'r')

#read the title
head = txt.readline()
head = head.split("\n")[0]
hd = head.split("\t")
#record sample names
for i in range(4,len(hd)):
    sampname_chi.append(hd[i])

#edit each line
for line in txt:
    print(line)
    line = line.split("\n")[0]
    ln = line.split("\t")
    
    #adjust each column contents
    chr1 = ln[0]
    pos1 = int(ln[1])
    chr2 = ln[2]
    pos2 = int(ln[3])
    
    #select the largest chimeric reads sample
    samchi = 0
    for a in range(len(sampname_chi)):
        if int(ln[a+4]) > samchi:
            samchi = int(ln[a+4])
            chifile = sampname_chi[a]
        else:
            continue
        
    matches = []

    #read in list of sam files, must be in current directory and end in "_aln.sam"
    li = os.listdir(os.getcwd())
    fileset = filter(lambda x: x.endswith('.sam'), li)
    for samfile in fileset:
        if chifile in samfile:
            fname = samfile

    #open target file
    f = open(fname)
    while 1:
        l = f.readline()
        #edit data that is not reads sequence
        if l == '':
            break
        if l[0] == '@' and l[1] == "S":
            h = l.split()
            name = h[1].replace("SN:","")
            nlen = int(h[2].replace("LN:",""))
            chromlist[name] = nlen
            continue
        if l[0] == "@":
            continue
        #edit data that is reads sequence
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
            if chr1 == chr2:
                if each[0][1] == each[1][1] == chr1:
                    if (each[0][2] > pos1 and each[0][2] < (pos1 + binsize)) and ((each[1][2] > pos2) and (each[1][2]) < (pos2 + binsize)):
                        combos.append(each)
                    if (each[1][2] > pos1 and each[1][2] < (pos1 + binsize)) and ((each[0][2] > pos2) and (each[0][2]) < (pos2 + binsize)):
                        combos.append(each)
            elif chr1 != chr2:
                if each[0][1] == chr1 and each[1][1] == chr2:
                    if (each[0][2] > pos1 and each[0][2] < (pos1 + binsize)) and ((each[1][2] > pos2) and (each[1][2]) < (pos2 + binsize)):
                        combos.append(each)
                if each[1][1] == chr1 and each[0][1] == chr2:
                    if (each[1][2] > pos1 and each[1][2] < (pos1 + binsize)) and ((each[0][2] > pos2) and (each[0][2]) < (pos2 + binsize)):
                        combos.append(each)
            else: 
                continue
            
        #put selected reads in table
        for item in combos:
            matches.append(list(item))

    #close chifile and prepare to work for the next line in txtfile
    f.close()
        
    #put selected reads name and sequences into two dictories
    for item in matches:
        if (chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)) not in list_seq:
            list_name[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] = [item[0][0]]
            list_name[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] = [item[1][0]]
            list_seq[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] = [item[0][4]]
            list_seq[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] = [item[1][4]]
        elif (chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)) in list_seq:
            list_name[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] += [item[0][0]]
            list_name[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] += [item[1][0]]
            list_seq[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] += [item[0][4]]
            list_seq[chr1+'+'+str(pos1)+'_'+chr2+'+'+str(pos2)] += [item[1][4]]

  
    
#write tables into fasta form
for i in list_name:
    o = open("seeds" + i + ".fasta", "w")
    for a in range(len(list_name[i])):
        o.write(">" + list_name[i][a] + "\n" + list_seq[i][a] + "\n")
    o.close()
  



   
