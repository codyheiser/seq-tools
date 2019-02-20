#!/bin/bash

# Author: C Heiser
# Date: 18Feb2018

# fastq preprocessing, QC

# Usage: path_to_script [RAW_DIR]
#	RAW_DIR = the relative path from cwd to directory containing .fastq.gz files for analysis
#
#	outputs are generated in cwd regardless of location of script or RAW_DIR

FILES=`ls $1/*.{fastq,fastq.gz,fq,fq.gz} 2>/dev/null` # get names of all FASTQ files from RAW_DIR 

wd=$(pwd)

# make directories for outputs
mkdir fastQC

printf "\nBegin FASTQ Preprocessing\n"

# loop through all .fastq.gz files
for f in $FILES;
do
	cd $wd # ensure you start in current working directory
	
	file=`basename $f` # get name of FASTQ file instead of path

	# perform fastQC on file $f in directory $wd/$1/ and write to ./1_fastQC/
	printf "\nConducting FastQC analysis on $f\n"
	fastqc -o . --noextrac -t 1 $wd/$1/$file --outdir=fastQC

done

printf "\nDone!\n"

