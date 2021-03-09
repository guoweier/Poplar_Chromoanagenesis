# Enrichment Ratio
### Introduction
This is a [pipeline](https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/enrichment-ratio/enrichment-window-ratio.py) for calculating the genomic feature density around each of the novel breakpoint, then run a one-sample t-test to check the significance between real breakpoints and pseudo-junction loci. 
### Working Process
#### Calculate density
1. Set a window on the genome with the breakpoint at the center. Here we run two different window sizes: 10kb and 100kb. 
2. Based on the Populus annotation file (downloaded from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#)), we count the base pair number of corresponded genomic features. Here we test two different features: gene and repeated elements. 
3. The feature density is calculated as: Density = bp_of_feature / window_size
4. Do this analysis for all the novel breakpoints. 
#### Set population group
1. We use the previously constructed pseudo-junction pool, and randomly selected 1,000 of these pseudo breaks for genomic feature density calculation. 
2. For each sample, this type of random pseudo-break datasets were established 1,000 times for every examined genomic feature. 
3. For every pseudo-breakpoint, we calculate the density. 
#### Statistical Analysis
1. Enrichment ratios were calculated by taking the means of genomic feature density at real breakpoints, divided by the means of the corresponded features density at random pseudo breakpoints datasets.
2. Run [one-sample t-test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_1samp.html) with sample group and population group. 
### Output File
The script has two output .txt files. 
1. The file with feature density for each novel breakpoint.
2. The file with 1,000 random datasets, followed with the mean of 1,000 pseudo-breakpoint density and p-value after t-test. 
