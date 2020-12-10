# Weier Guo bash script

#Introduction
#A bash script for PRICE assembly.
#It takes all the seed files in the current folder, runs PRICE, creates contigs. 

#Input
#Seed files: 
#A fasta file containing cross-junction reads.
#Remeber all seed files should be in one directory, and the directory only contain seed files. 

# read all files in the current folder
for file in *
do
	Allfiles+=($file)
done

# record name
for sam in ${Allfiles[@]}
do 
	samfile=$(echo $sam | cut -d'.' -f 1)
	Mysamfiles+=($samfile)
done

# run PRICE
for seed in ${Mysamfiles[@]}
do
	PriceTI -spfp ${seed}.fasta 50 150 95 -icf ${seed}.fasta 1 1 1 -nc 10 -mol 20 -tol 10 -mpi 80 -o ${seed}.fasta
done

