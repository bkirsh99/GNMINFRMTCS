#!/bin/bash
#
# AUTHOR: Bianca Kirsh <bkirsh.3006@gmail.com>
# DATE: January 26, 2021
# FILE: fastqToFilteredBam.sh
#
# SYNOPSIS: A pipeline created for Assignment 2 (BMEG 400E).
#
# DESCRIPTION: This script pre-processes paired-end raw sequence data provided in FASTQ format to produce
# analysis-ready BAM files containing only uniquely mapped reads. This includes quality control using
# FastQC, alignment to a reference genome using Bowtie 2, and finally conversion of the SAM file into a
# sorted and filtered BAM file using Samtools and Sambamba.
#
# REQUIREMENTS: This scripts requies a CWL account and a connection to the course server (137.82.55.186)
# while logged into the UBC Virtual Private Network (VPN). Furthermore, it must be called from inside
# an active conda environment with the following packages installed:
#       - fastqc (https://anaconda.org/bioconda/fastqc)
#       - bowtie2 (https://anaconda.org/bioconda/bowtie2)
#       - samtools (https://anaconda.org/bioconda/samtools)
#       - sambamba (https://anaconda.org/bioconda/sambamba)
#
# USAGE: To run this script, simply enter the following command prompt in a unix shell:
#
#       ./fastqToFilteredBam.sh <sampleID> <fastqFile1> <fastqFile2>
#
#       Alternatively, it can be passed as a command-line argument to the program available
# for download at https://raw.githubusercontent.com/BMEGGenoInfo/Assignments/main/Assignment_2/runTheseJobsSerially.sh.
#
#       ./runTheseJobsSerially.sh ./fastqToFilteredBam.sh <taskfile>
#
#       In this case, the <taskfile> consists of a tab-delimited file with three columns:
# sampleID, fastqFile1, and fastqFile2. A sample first line is illustrated below.
#
#       iPSC_input  input_iPSC_SRA66_subset_1.fastq.gz  input_iPSC_SRA66_subset_2.fastq.gz
#
# OUTPUTS: This script creates two new directories within your current working directory: "MyLogDirectory" and
# "MyOutputDirectory." While the first holds a set of empty files that represent the internal state of the pipeline,
# the latter contains the actual output of each step. It consists of a two-level directory whose master file directory
# is indexed by the sampleID. The user file directories include the following files:
#       - <fastqFile1>.html, <fastqFile2>.htm, <fastqFile1>.zip, and <fastqFile2>.zip: Outputs of FastQC step
#       - <sampleID>.sam: Output of Bowtie 2 alignment step
#       - <sampleID>.bam: Output of Samtools conversion step
#       - <sampleID>.sorted.bam and <sampleID>.sorted.bai: Outputs of Sambamba sorting step
#       - <sampleID>.filtered.bam and <sampleID>.filtered.bai: Outputs of Sambamba filtering step
#
# SPECIAL NOTES: For the purposes of this assignment, the workflow is designed to operate using pre-downloaded
# fastq files (/usr/local/share/data/assignment_2) and a pre-indexed human genome (/usr/local/share/indexes/hg38_bowtie2_index).
#

set -Eeuo pipefail # Abort the script on errors and unbound variables.

# Variables to hold the arguments given to this script
sample=$1
fq1=$2
fq2=$3
# Variables to indicate input and output file path locations
logDir=MyLogDirectory # this is where all the files to keep track of progress will go.
filePath=/usr/local/share/data/assignment_2 # this is the path to the pre-downloaded fastq files.
genome=/usr/local/share/indexes/hg38_bowtie2_index # this is the path to the pre-indexed genome file.
outDir=MyOutputDirectory/$sample # this is where the output files from the workflow will go.

mkdir -p MyLogDirectory # make the directory where log files will go, if it doesn't exist already
mkdir -p $outDir # make the directory where output files will go, if it doesn't exist already

echo Running pipeline for $sample

if [ ! -e $logDir/$sample.fastqc.done ] #run this code only if $logDir/$sample.fastqc.done is missing
then
        echo Performing fastqc of sample $sample with the following fastqs:
        ls $filePath/$fq1 $filePath/$fq2
        # Run a FastQC analysis on each of the input files.
        for file in $fq1 $fq2; do
                fastqc $filePath/$file -o $outDir
        done
        touch $logDir/$sample.fastqc.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.fastqc.done was not missing
        echo Already performed fastqc of $sample
fi

if [ ! -e $logDir/$sample.sam.done ] #run this code only if $logDir/$sample.sam.done is missing
then
        echo Performing alignment of sample $sample with the following indexed genome file: $genome
        bowtie2 -x $genome -1 $filePath/$fq1 -2 $filePath/$fq2 -S $outDir/$sample.sam
        touch $logDir/$sample.sam.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.sam.done was not missing
        echo Already performed alignment of $sample
fi

if [ ! -e $logDir/$sample.bam.done ] #run this code only if $logDir/$sample.bam.done is missing
then
        echo Performing SAM to BAM conversion of sample $sample
        samtools view -h -b -o $outDir/$sample.bam $outDir/$sample.sam
        touch $logDir/$sample.bam.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.bam.done was not missing
        echo Already performed SAM to BAM conversion of $sample
fi

if [ ! -e $logDir/$sample.sorted.done ] #run this code only if $logDir/$sample.sorted.done is missing
then
        echo Sorting the BAM file $sample.bam
        sambamba sort $outDir/$sample.bam
        touch $logDir/$sample.sorted.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.sorted.done was not missing
        echo Already sorted the BAM file for $sample
fi

if [ ! -e $logDir/$sample.filtered.done ] #run this code only if $logDir/$sample.filtered.done is missing
then
        echo Filtering BAM file $sample.bam to remove ambiguously mapped reads.
        sambamba view -F "[XS] == null and not unmapped and not duplicate" -h -f bam -o $outDir/$sample.filtered.bam $outDir/$sample.sorted.bam
        touch $logDir/$sample.filtered.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.filtered.done was not missing
        echo Already filtered the BAM file for $sample
fi

echo Pre-processing complete for sample $sample. The output files are located in: $outDir
