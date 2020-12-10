# Poplar Chromoanagenesis
A pipeline for locating the genomic position of chromosomal rearrangement junctions by using Illumina genome sequencing. 

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
There are 4 parts for this pipeline. 

#### Part 1. Record dosage Variation boundaries.
1. Input file is the outcome of [Bin-By-Sam.py](https://github.com/Comai-Lab/bin-by-sam) from Comai Lab. Check the link for more inforamtion. 
2. Make scatter plots with the input data. The plot is going to show relative read coverage of consecutive set-size bins along the whole genome. 
3. Manually select boundary bins of copy number variation.

#### Part 2. Extract cross-junction reads. 
1. To find out whether there are novel junctions, the program look for regions where have sequencing reads that mapped at least 2,000 bp apart on the genome. 
2. The outcome of this script will generate a list of genomic bin pairs, which represent two far apart genomic locations that have reads mapped simultaneously on them. It follows by the number of cross-junction reads in each examine sample.   
2. To avoid false positive, we need to remove the noises by setting two thresholds. 
3. First threshold: The novel junctions are unique to specific sample, so other peer samples and controls should have 0 cross-junction reads. 
4. Second threshold: The sample containing cross-junction reads should have enough coverage to demonstrate it may be a real junction. Please look at [intro-adjacent-bin-threshold.md](https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/intro-adjacent-bin-threshold.md) for detailed information.
6. To sum up, only when one sample have cross-junction reads over the baseline for this sample, meanwhile other peer samples and controls have 0 cross-junction reads, this genomic region pairs will become a potential novel junction location. 
7. Since every sample has its own pseudo-junction coverage, the program did threshold selection one by one. Finally, every examine sample will have its unique list of potential novel junction location. 

#### Part 3. Compare dosage variation boundaries with cross-junction reads regions. (optional)
1. This step is to typically compare regions containing cross-junction reads with dosage boundaries. If only interested in the question like 'where do CNVs insert in the genome', this step is helpful. 

#### Part 4. PRICE assembly and BLAST
1. PRICE assembly can assemble cross-junction reads into one contig. Look at [PRICE](http://derisilab.ucsf.edu/software/price/) for more information. 
2. For running PRICE assembly, every potential junction needs a seed file, which is a fasta file containing its cross-junction reads. A custom [python script](https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/pairend-seeds-general.py) was designed for extracting cross-junction reads. Briefly, it takes a list of potential junction regions, and searched in .sam files to extract typical reads. <br>
   - The input list example:<br>
    ```
    Ref1	Bin1-Start	Ref2	Bin2-Start	POP25_72	POP26_09	POP26_54	POP27_32	POP27_77	POP27_88	POP28_09	POP28_86	POP30_88	POP31_79	POP33_31	female	male
    Chr01	1274500	Chr16	3732500	11	0	0	0	0	0	0	0	0	0	0	0	0
    ```
   - The input list starts with the two crossed regions, then followed by cross-junction reads number for each sample of this specific region. <br>
3. A custom [bash script](https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/Price-assembly-loop.sh) was designed to run PRICE assembly for multiple tasks at once. <br>
   - The fundamental PRICE command example:<br>
    ```
    PriceTI -spfp ${seed}.fasta 50 150 95 -icf ${seed}.fasta 1 1 1 -nc 10 -mol 20 -tol 10 -mpi 80 -o ${seed}.fasta
    ```
   - -spfp a b c d: unpaired reads are split into paired ends, with the scores of double-use nuceotides halved. (a)input file, (b)the length of the 'reads' that will be taken from each side of the input reads, (c)amplicon insert size (including read), (d)required % identity for match (25-100 allowed)
   - -icf a b c: input seed file. (a)initial contig file, (b)number of addition steps, (c)number of cycles per step, (d)const by which to multiply quality scores
   - -nc a: (a)num. of cycles
   - -mol a: (a)minimum overlap length for mini-assembly
   - -tol a: (a)threshold seq num for scaling overlap for contig-edge assemblies
   - -mpi a: (a)minimum % identity for contig-edge assembly
   -o a: (a)output file name (.fasta)
   - The outcome of PRICE script: For every potential junction, there is a fasta file containing its integrated contigs. <br>
4. After PRICE assembly, we blast contigs onto the reference genome. A custom [bash script](https://github.com/guoweier/Poplar_Chromoanagenesis/blob/master/BlastN-multiquery.sh) was designed for running BLAST nucleotides for multiple tasks at once. <br>
   - The fundamental blastN command example:<br>
      ```
      blastn -db /cato2pool/weier-poplar/raw/ref/Ptrichocarpa_v3.0_210.fa -query ${fasta}.cycle10.fasta -out ${fasta}.out
      ```
   - -db: database
   - -query: input contig files
   - -out: output file
   - The criteria for selecting valid junctions:<br>
      a. The contig mapped onto predicted regions on the reference genome
      b. The contig did not map to other genomic positions multiple times (>=4)

### Results
This pipeline predicts the genomic positions of novel DNA junctions by using Illumina short-read sequencing. It finally gives out potential junction sequence, the two genomic mapping positions, as well as other information such as junction orientation. <br>
With these results, we can predict how genome being rearranged during early development. Structural variation such as transloction, inversion, or even extreme scenario such as chromoanagenesis can be characterized with this approach. <br>

