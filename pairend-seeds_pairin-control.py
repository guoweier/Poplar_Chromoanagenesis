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
#Potential junctions in chimeric reads format with no title line. i.e. Ref1 Ref1-BinStart Ref2 Ref2-BinStart CK sample1 sample2 ... samplen 
#Each element is separated by tab. The file is txt. 

#Output file:
#fasta file
#1st line: potential junction; list all seeds reads within fasta file format (name + sequence) 


usage = ""
parser = OptionParser(usage=usage)
#parser.add_option("-f", "--file", dest="f", help="Input .sam file.")
parser.add_option("-f", "--input_file", dest="f", help="Text file for selecting chimeric reads.")
#parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=5000, help="Bin size (5000)")
parser.add_option("-m", "--mapqual", dest="mapqual", type = "int", default=10, help="Minimum mappign quality (10)")
parser.add_option("-i", "--insert", dest="insert", type = "int", default=2000, help="maximum estimated insert between forward and reverse (2000)")
parser.add_option("-c", "--mincov", dest="mincov", type = "int", default=0, help="Minimum coverage for a bin to be in the output, across all libraries")


(opt, args) = parser.parse_args()


txtfile = opt.f
mindiff = opt.insert
minmap = opt.mapqual
(opt, args) = parser.parse_args()

chromlist = defaultdict()
bins = defaultdict(lambda: defaultdict(int))

#set output dictories
list_name = {}
list_seq = {}

#open txt file
txt = open(txtfile,'r')

#edit each line
for ln in txt:
    print(ln)
    ls = ln.split('\t') #split each line
    
    #adjust each column contents
    l0 = ls[0]
    l1 = int(ls[1])
    l2 = ls[2]
    l3 = int(ls[3])
    l4 = int(ls[4])
    l5 = int(ls[5])
    l6 = int(ls[6])
    l7 = int(ls[7])
    
    #select the largest chimeric reads sample
    sample_chi = [l4,l5,l6,l7]
    fname_chi = ['rePOP25_72_aln.sam','rePOP26_54_aln.sam','rePOP27_88_aln.sam','rePOP28_86_aln.sam']
    samchi = 0
    for a in range(len(sample_chi)):
        if sample_chi[a] > samchi:
            samchi = sample_chi[a]
            chifile = fname_chi[a]
        else:
            continue
        
    matches = []
    
    f = open(chifile)
    while 1:
        l = f.readline()
        #edit data that is not reads sequence
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
            if l0 == l2:
                if each[0][1] == each[1][1] == l0:
                    if (each[0][2] > l1 and each[0][2] < (l1 + 500)) and ((each[1][2] > l3) and (each[1][2]) < (l3 + 500)):
                        combos.append(each)
                    if (each[1][2] > l1 and each[1][2] < (l1 + 500)) and ((each[0][2] > l3) and (each[0][2]) < (l3 + 500)):
                        combos.append(each)
            elif l0 != l2:
                if each[0][1] == l0 and each[1][1] == l2:
                    if (each[0][2] > l1 and each[0][2] < (l1 + 500)) and ((each[1][2] > l3) and (each[1][2]) < (l3 + 500)):
                        combos.append(each)
                if each[1][1] == l0 and each[0][1] == l2:
                    if (each[1][2] > l1 and each[1][2] < (l1 + 500)) and ((each[0][2] > l3) and (each[0][2]) < (l3 + 500)):
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
        if (l0+'+'+str(l1)+'_'+l2+'+'+str(l3)) not in list_seq:
            list_name[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] = [item[0][0]]
            list_name[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] = [item[1][0]]
            list_seq[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] = [item[0][4]]
            list_seq[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] = [item[1][4]]
        elif (l0+'+'+str(l1)+'_'+l2+'+'+str(l3)) in list_seq:
            list_name[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] += [item[0][0]]
            list_name[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] += [item[1][0]]
            list_seq[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] += [item[0][4]]
            list_seq[l0+'+'+str(l1)+'_'+l2+'+'+str(l3)] += [item[1][4]]

  
    
#write tables into fasta form
for i in list_name:
    o = open("seeds" + i + ".fasta", "w")
    for a in range(len(list_name[i])):
        o.write(">" + list_name[i][a] + "\n" + list_seq[i][a] + "\n")
    o.close()
  



   