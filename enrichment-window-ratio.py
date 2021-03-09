# Weier Guo Python

from __future__ import division
import os, sys, math
from optparse import OptionParser
from scipy import stats
import random2 as random 
import statistics

# Intro:
# this script is used for measuring enrichment ratio.
# for each selected junction, two jointed ends are treated separately, with each one being set as the center of a window.
# gene density = bp of gene in window / total size of window * 100

# Usage: gene-density-window-ratio.py -f gene-file-from-gff.txt -j junction-file.txt -g genome-size-file.txt -s window-size -r random-select-time -d randome-dataset-number -o output.txt

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f", "--input_file", dest="f", help="input annotation file simplified from gff")
parser.add_option("-j", "--junction", dest="j", help="input junction file")
parser.add_option("-F", "--fakejunction", dest="fakejun", help="input adjacent junction file")
parser.add_option("-n", "--samplename", dest="name", help="sample name")
parser.add_option("-s", "--size", dest="s", type="int", default=10000, help="window size (default=10000)")
parser.add_option("-r", "--randomtime", dest="r", type="int", default=10000, help="the times of randomly select breakpoints on the genome (default=10000)")
parser.add_option("-d", "--dataset", dest="d", type="int", default=1000, help="the times of running times for running the random sets (default=1000)")
parser.add_option("-o", "--output", dest="o", help="output file")

(opt, args) = parser.parse_args()

#open files and set parameters
f = open(opt.f)
j = open(opt.j)
fjun = open(opt.fakejun)
samplename = opt.name
s = opt.s
rtime = opt.r
dataset = opt.d
o = open(opt.o, "w")
enrichs = {}
enrichdens = {}
sampdens = []
Chrom = []
Genome = {}
Randden = {}
oRandden = []

#put annotation file into dictionary
fhead = f.readline()
for line in f:
	line = line.split("\n")[0]
	ln = line.split("\t")
	chrm = ln[0]
	start = int(ln[1])
	end = int(ln[2])
	if chrm not in enrichs:
		enrichs[chrm] = []
	enrichs[chrm].append([start,end])

#function for enrichment ratio in selected window
def enrich_den_window(breakchr, breakpos):
	left = breakpos-(s/2)
	right = breakpos+(s/2)
	jun = []
	enrichbp = 0
	#check enrichs bp with annotation file
	for item in enrichs:
		if breakchr	== item:
			for enrich in enrichs[item]:
				if enrich[0] < left and enrich[1] < left:
					continue
				elif enrich[0] < left and enrich[1] > right:
					jun.append([left,right])
				elif enrich[0] < left and enrich[1] >= left and enrich[1] <= right:
					jun.append([left,enrich[1]])
				elif enrich[0] >= left and enrich[0] <= right and enrich[1] >= left and enrich[1] <= right:
					jun.append([enrich[0],enrich[1]])
				elif enrich[0] >= left and enrich[0] <= right and enrich[1] > right:
					jun.append([enrich[0],right])

				else:
					continue
	#calculate enrichment ratio
	for g in jun:
		enrichbp += g[1]-g[0]
	den = float(enrichbp/s)
	if den > 1:
		den = 1.0
	return den

#apply function on junction list
jhead = j.readline()
for line in j:
	line = line.split("\n")[0]
	ln = line.split("\t")
	chrom = ln[0]
	pos = int(ln[1])
	junname = chrom+"-"+str(pos)
	junden = enrich_den_window(chrom, pos)
	enrichdens[junname] = junden
	sampdens.append(junden)

#set genome for random selection
fjun_head = fjun.readline()
for line in fjun:
	line = line.split("\n")[0]
	ln = line.split("\t")
	if ln[0] not in Genome:
		Genome[ln[0]] = []
	Genome[ln[0]] += [int(ln[1])]
	Chrom.append(ln[0])

#randomly select breaks for set times, calculate for enrichment ratio
for i in range(dataset):
	print(i)
	Randden[str(i)] = []
	for ji in range(rtime):
		rchrom = random.choice(Chrom)
		rpos = random.choice(Genome[rchrom])
		randomden = enrich_den_window(rchrom, rpos)
		Randden[str(i)]+=[randomden]

	#run one-sample T-test
	sampmean = statistics.mean(sampdens)
	randmean = statistics.mean(Randden[str(i)])
	testres = stats.ttest_1samp(sampdens, randmean)
	pvalue = testres[1]
	oRandden.append([str(i),str(randmean),str(pvalue)])
	

#output
#breakpoints enrichment ratio
ohead = "Chrom\tPos\tEnrichmentRatio"
o.write(ohead+"\n")
for item in enrichdens:
	ch = item.split("-")[0]
	position = item.split("-")[1]
	text = ch+"\t"+position+"\t"+str(enrichdens[item])
	o.write(text+"\n")

#random dataset summary
op = open(samplename+"-"+"Stats.txt", "w")
op.write("Dataset\tRandom-mean\tP-value\n")
for item in oRandden:
	op.write(item[0]+"\t"+item[1]+"\t"+item[2]+"\n")

f.close()
j.close()
fjun.close()
o.close()
op.close()




