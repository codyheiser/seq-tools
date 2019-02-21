#!/bin/bash

# Author: C Heiser
# Date: 20Feb2018

# get reverse complement of sequences in file 

# Usage: path_to_script [FASTQ]
#	FASTQ = the relative path to sequence file for analysis

grep '^[ATCG]' $1 | rev | tr ATCG TAGC
