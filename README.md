# Poplar Chromoanagenesis
A workflow for locating the genomic position of chromosomal rearrangement junctions by using Illumina genome sequencing. 

## Introduction
Genomic structural variation is widely reported in different species. Extreme scenario such as chromoanagenesis, which is the fragmentation and reorganization of one chromosome in early mitosis, has been detected under various situation. It has shown to correlate with human cancer, but also have been found in plants. 
Identifying the novel DNA junctions would be really essential to study genome restructuring, because it provides the foundamental information for studying:
a. How many chromosomes involved;
b. How does it affect genomic features, such as genes, TEs, etc. 
c. Discover other features such as chromosome inversion, translocation. 

By using the Illumina short-read sequencing data, this workflow is able to locate the novel DNA junctions during genome restructuring. 
It is able to find out: 
a. Where are the two ends originally mapped on the reference genome. 
b. The surrounding 100-200bp sequence of novel DNA junction. 
c. The overlaps or unknown insertions between two joined fragment ends. 
