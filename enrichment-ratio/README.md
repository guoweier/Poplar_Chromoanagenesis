# Enrichment Ratio
### Introduction
This is a pipeline for calculating the genomic feature density around each of the novel breakpoint, then run a one-sample t-test to check the significance between real breakpoints and pseudo-junction loci. 
### Working Process
#### Set sample group
1. Set a window on the genome with the breakpoint at the center. Here we run two different window sizes: 10kb and 100kb. 
2. Based on the Populus annotation file (downloaded from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#)), we count the base pair number of corresponded genomic features. Here we test two different features: gene and repeated elements. 
3. The feature density is calculated as: Density = bp_of_feature / window_size
4. Do this analysis for all the novel breakpoints. 
#### Set population group
1. We use the previously constructed pseudo-junction pool, and randomly selected 1,000 of these pseudo breaks for genomic feature density calculation. 
