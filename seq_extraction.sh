#!/bin/bash

# Author: C Heiser
# Date: 20Feb2018

# extract sequence between start and stop sedquences in a FASTQ file 

# Usage: path_to_script [FASTQ] [START_SEQ] [STOP_SEQ]
#	FASTQ = the relative path from cwd to .fastq.gz file for analysis
#	START_SEQ = string to start grabbing sequence from (excluded from output)
#	STOP_SEQ = string to stop grabbing sequence after (excluded from output)
#
#	outputs are generated in cwd regardless of location of script or FASTQ

printf "\nBegin Sequence Extraction from $1 between $2 and $3\n"

file=`basename $1` # get name of FASTQ file instead of path

if [[ $1 =~ \.gz$ ]] # if file is gzipped, use gunzip before awk command
then
	time gunzip -cd $1 | awk '{ if (match($0,/'$2'(.+?)'$3'/,m)) print m[1] }' > ${file%.fastq.gz}.txt
else # if file is not gzipped, run awk normally
	time awk '{ if (match($0,/'$2'(.+?)'$3'/,m)) print m[1] }' $1 > ${file%.fastq.gz}.txt
fi

printf "\nDone!\n"
