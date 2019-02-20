#!/bin/bash

# Author: C Heiser
# Date: 04Oct2018

# fastq preprocessing, QC, alignment, and sorting for ChIP-seq data

# Usage: path_to_script [RAW_DIR] [genome]
#	RAW_DIR = the relative path from cwd to directory containing .fastq.gz files for analysis
#	genome = reference genome to use (e.g. mm10, hg19)
#
#	outputs are generated in cwd regardless of location of script or RAW_DIR

FILES=`ls $1/*.{fastq,fastq.gz,fq,fq.gz} 2>/dev/null` # get names of all FASTQ files from RAW_DIR 

# define resources 
if [ "$2" = "hg19" ] # if given genome is human, set genome and chromosome lengths accordingly
then
	genome="/path_to/hg19/genome.fa" # TODO: update path to resources
	chromlens="/path_to/hg19.chr1-22xyM.txt" # TODO: update path to resources
elif [ "$2" = "mm10" ] # if given genome is mouse, set genome and chromosome lengths accordingly
then
	genome="/path_to/mm10/genome.fa" # TODO: update path to resources
	chromlens="/path_to/mm10.chr1-19xyM.txt" # TODO: update path to resources
else # if valid genome name is not given, warn user
	printf "\nPlease provide a valid genome for alignment (e.g. mm10, hg19)\n"
	break
fi

bedtools="/path_to/bedtools" # TODO: update path to resources
ucsc_bg2bw="/path_to/UCSC_tools/bedGraphToBigWig" # TODO: update path to resources
bamcoverage="/path_to/bamCoverage" # TODO: update path to resources
wd=$(pwd)

# make directories for outputs
mkdir 1_fastQC
mkdir 2_bam
mkdir 3_unique_bam
mkdir 4_fingerprintplot
mkdir 5_deepTOOLS
mkdir 6_ChAsE

printf "\nBegin FASTQ Preprocessing\n"

# loop through all .fastq.gz files
for f in $FILES;
do
	cd $wd # ensure you start in current working directory
	
	file=`basename $f` # get name of FASTQ file instead of path

	# perform fastQC on file $f in directory $wd/$1/ and write to ./1_fastQC/
	printf "\nConducting FastQC analysis on $f\n"
	fastqc -o . --noextrac -t 1 $wd/$1/$file --outdir=1_fastQC
	
	# perform Burrows-Wheeler alignment on file $f in directory $wd/$1/ and save as ./2_bam/f.bam
	printf "\nAligning file $f to $2 genome\n"
	time bwa mem -t 8 $genome $wd/$1/$file | samtools view -bS > $wd/2_bam/${file%.fastq.gz}.bam -
	
	# filter for only uniquely mapped reads in each .bam file
	printf "\nFiltering for uniquely mapped reads in .bam file ${f%.fastq.gz}.bam\n"
	time samtools view -F 4 -q 1 -hb $wd/2_bam/${file%.fastq.gz}.bam > $wd/3_unique_bam/${file%.fastq.gz}_unique.bam
	
	# sort reads in unique .bam file
	printf "\nSorting uniquely mapped reads in .bam file ${file%.fastq.gz}_unique.bam\n"
	cd $wd/3_unique_bam # change to unique bam directory
	time samtools sort -@8 ${file%.fastq.gz}_unique.bam -o ${file%.fastq.gz}_unique_sort.bam
	
	# create index (.bam.bai) file for the sorted, unique .bam file
	printf "\nCreating index file for ${file%.fastq.gz}_unique_sort.bam\n"
	time samtools index ${file%.fastq.gz}_unique_sort.bam
	
	# create fingerprint plot of genome coverage and save in ./4_fingerprintplot/
	printf "\nCreating fingerprint plot for ${file%.fastq.gz}\n"
	time plotFingerprint -b ${file%.fastq.gz}_unique_sort.bam -plot $wd/4_fingerprintplot/${file%.fastq.gz}_unique_sort.pdf
	
	# create BigWig files using bamcoverage (path defined at top of script)
	printf "\nCreating bigWig for deepTOOLS analysis\n"
	time $bamcoverage --numberOfProcessors max --normalizeUsingRPKM --binSize 20 --smoothLength 100 --verbose --bam ${file%.fastq.gz}_unique_sort.bam -o $wd/5_deepTOOLS/${file%.fastq.gz}.b20s100dt.bw
	
	# create BigWig using ucsc_bg2bw and bedtools (paths defined at top of script)
	printf "\nCreating bigWig for ChAsE visualization\n"
	time $bedtools genomecov -split -bg -g $chromlens -ibam ${file%.fastq.gz}_unique_sort.bam > temp1.bam.bedgraph # generate bedgraph
	sort -k1,1 -k2,2n temp1.bam.bedgraph > temp2.sort.bam.bedgraph # sort bedgraph
	$ucsc_bg2bw temp2.sort.bam.bedgraph $chromlens $wd/6_ChAsE/${file%.fastq.gz}.bw # generate .bw file
	rm temp* # remove unneccesary bedgraph intermediate files

done

printf "\nDone!\n"

