# Weier Guo bash script

# Introduction: 
# A bash script for blast nucleotide with custom database.
# Input:
# Query: sequence waiting for blast. It must be in fasta format. 
# Remeber all query files should be in one directory, and the directory only contain seed files but no others. 
# Database: custom database built by own reference genome. Here, I am using Populus trichocarpa v3.0 (I am trying to find a way that can input commands into bash script). 
# Output:
# Every fasta file will have one output file. It includes similar format of blast nucleotide results as on NCBI.

# More details: https://www.ncbi.nlm.nih.gov/books/NBK279680/

for file in *
do
	Allfiles+=($file)
done

for fasta in ${Allfiles[@]}
do 
	fastafile=$(echo $fasta | cut -d'.' -f 1)
	Myfastafiles+=($fastafile)
done

for fasta in ${Myfastafiles[@]}
do
	blastn -db /cato2pool/weier-poplar/raw/ref/Ptrichocarpa_v3.0_210.fa -query ${fasta}.cycle10.fasta -out ${fasta}.out
done
