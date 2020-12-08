#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Count Cameron et al. (2013) reads mapped by HISAT2 to genes using 
# HTSeq-count.
#-------------------------------------------------------------------------------
# Inputs:
# HISAT2 alignments in .bam format, which have been sorted. Apis mellifera (Amel)
# GTF file.
#
# Outputs:
# Text files with read counts per gene.
#-------------------------------------------------------------------------------

# LOAD MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load python/anaconda/2019.10/3.7  # where HTSeq is installed.

# STEP 1: MAKE OUTPUT DIRECTORIES FOR RESULTS ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

mkdir ../02_outputs/30_cameron_htseq_count

# STEP 2: COUNT READS PER GENE IN BAM FILES ----

# Change working directory.

cd ../02_outputs/20_cameron_hisat2_aligned_reads

# Perform counts.
# Args:
# -f: format of the input data.

for NAME in *hisat2.bam
do
	OUTFILE=${NAME:0:11}
htseq-count -f bam ${NAME} ../../00_data/02_gtf/Apis_mellifera.Amel_4.5.48.gtf.gz > ../30_cameron_htseq_count/${OUTFILE}_HTSeq.txt
done
