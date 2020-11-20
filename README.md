# Poplar Chromoanagenesis
A workflow for locating the genomic position of chromosomal rearrangement junctions by using Illumina genome sequencing. 

### Introduction
Genomic structural variation is widely reported in different species. Extreme scenario such as chromoanagenesis, which is the fragmentation and reorganization of one chromosome in early mitosis, has been detected under various situation. It has shown to correlate with human cancer, but also have been found in plants. <br>

Identifying the novel DNA junctions would be really essential to study genome restructuring, because it provides the foundamental information for studying: 
1. How many chromosomes involved;
2. How does it affect genomic features, such as genes, TEs, etc. 
3. Discover other features such as chromosome inversion, translocation, etc. 

By using the Illumina short-read sequencing data, this workflow is able to locate the novel DNA junctions during genome restructuring. <br>
It is able to find out: 
1. Where are the two ends originally mapped on the reference genome. 
2. The surrounding 100-200bp sequence of novel DNA junction. 
3. The overlaps or unknown insertions between two joined fragment ends. 

### Working Process
There are 4 parts for this workflow. 

#### Part 1. Record dosage Variation boundaries.
1. Input file is the outcome of Bin-By-Sam.py from Comai Lab (https://github.com/Comai-Lab/bin-by-sam). Check the link for more inforamtion. 
2. Make scatter plots with the input data. The plot is going to show relative read coverage of consecutive set-size bins along the whole genome. 
3. Manually select boundary bins of copy number variation.

#### Part 2. Extract cross-junction reads. 
1. To find out whether there are novel junctions, the program look for regions where have sequencing reads that mapped at least 2,000 bp apart on the genome. 
2. The outcome of this script will generate a list of genomic bin pairs, which represent two far apart genomic locations that have reads mapped simultaneously on them. It follows by the number of cross-junction reads in each examine sample.   
2. To avoid false positive, we need to remove the noises by setting two thresholds. 
3. First threshold: The novel junctions are unique to specific sample, so other peer samples and controls should have 0 cross-junction reads. 
4. Second threshold: The sample containing cross-junction reads should have enough coverage to demonstrate it may be a real junction. Please look at intro-adjacent-bin-threshold.md (https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/intro-adjacent-bin-threshold.md) for detailed information.
6. To sum up, only when one sample have cross-junction reads over the baseline for this sample, meanwhile other peer samples and controls have 0 cross-junction reads, this genomic region pairs will become a potential novel junction location. 
7. Since every sample has its own pseudo-junction coverage, the program did threshold selection one by one. Finally, every examine sample will have its unique list of potential novel junction location. 

#### Part 3. 
