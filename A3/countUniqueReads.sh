#!/bin/bash
#
# AUTHOR: Bianca Kirsh <bkirsh.3006@gmail.com>
# DATE: February 1, 2021
# FILE: countUniqueReads.sh
#
# SYNOPSIS: A mini-pipeline created for Assignment 3 (BMEG 400E).
#
# DESCRIPTION: This script uses Sambamba to get the number of uniquely mapped reads from SAM and/or
# BAM files that have already been trimmed and aligned against the hg38 genome build.
#
# REQUIREMENTS: This scripts requies a CWL account and a connection to the course server (137.82.55.186)
# while logged into the UBC Virtual Private Network (VPN). Furthermore, it must be called from inside
# an active conda environment with the following packages installed:
#       - sambamba (https://anaconda.org/bioconda/sambamba)
#
# USAGE: To run this script, simply enter the following command prompt in a unix shell:
#
#       ./fastqToFilteredBam.sh <path> <fileName> <readLength>
#
#       Alternatively, it can be passed as a command-line argument to the program available
# for download at https://raw.githubusercontent.com/BMEGGenoInfo/Assignments/main/Assignment_2/runTheseJobsSerially.sh.
#
#       ./runTheseJobsSerially.sh ./countUniqueReads.sh <taskfile>
#
#       In this case, the <taskfile> consists of a tab-delimited file with three columns:
# path, fileName, and readLength. A sample first line is illustrated below.
#
#       /usr/local/share/data/assignment_3     H3K27me3_iPSC_SRA60_subset_1_LEN_50_mapped.bam  50
#
# OUTPUTS: This script outputs the read length, as well as the the number of uniquely mapping reads.
#
# SPECIAL NOTES: For the purposes of this assignment, the workflow is designed to operate using pre-downloaded
# bam files (/usr/local/share/data/assignment_3).

set -Eeuo pipefail # Abort the script on errors and unbound variables.

# Variables to hold the arguments given to this script
path=$1 #path to file
file=$2 #file name
len=$3 #read length

echo Running pipeline for $file
echo -e "\t" Read length: $len

if [[ $file =~ \.sam$ ]]
then
        unique=`sambamba view -F "[XS] == null and not unmapped and not duplicate" -q -S -c $path/$file`
fi

if [[ $file =~ \.bam$ ]]
then
        unique=`sambamba view -F "[XS] == null and not unmapped and not duplicate" -q -c $path/$file`
fi

echo -e "\t" Uniquely mapping reads: $unique "\n"
