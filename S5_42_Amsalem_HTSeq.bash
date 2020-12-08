#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 3, 4 and 5
# Tasks: Count Amsalem et al. (2015) reads mapped by HISAT2 to genes using 
# HTSeq-count.
#-------------------------------------------------------------------------------
# Inputs:
# HISAT2 alignments in .bam format, which have been sorted. Bter GTF file.
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

mkdir ../02_outputs/32_amsalem_htseq_count

# STEP 2: COUNT READS PER GENE IN BAM FILES ----

# Change working directory.

cd ../02_outputs/24_amsalem_hisat2_aligned_reads

# Remove results for M_line207.
# Inspection of the FastQC analysis and HISAT2 alignment summary shows that this set
# of data should be excluded from further analysis as the mRNA-seq reads have are
# different in several respects compared to the other mRNA-seq libraries, and only
# a small percentage align to the Bombus terrestris genome.

rm -f M_line207_hisat2.bam

# Perform counts on remaining .bam files.
# Args:
# -f: format of the input data.

for NAME in *hisat2.bam
do
	OUTFILE=${NAME:0:9}
htseq-count -f bam ${NAME} ../../00_data/02_gtf/Bombus_terrestris.Bter_1.0.47.gtf.gz > ../32_amsalem_htseq_count/${OUTFILE}_HTSeq.txt
done
