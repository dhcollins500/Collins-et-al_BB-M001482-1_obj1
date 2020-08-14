#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Align RNA-Seq reads to Bombus terrestris genome using HISAT2.
#-------------------------------------------------------------------------------
# Inputs:
# HISAT2 index of the Bombus terrestris genome, raw RNA-Seq read files.
#
# Outputs:
# .txt files summarizing the stats of each alignment command. Sorted .bam files
# of the reads. MultiQC report of HISAT2 alignment statistics.
#-------------------------------------------------------------------------------

# LOAD MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add HISAT2/2.1.0
module add python/anaconda/2019.10/3.7  # where MultiQC is installed.
module add samtools/1.10

# STEP 1: MAKE OUTPUT DIRECTORIES FOR RESULTS ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

mkdir ../02_outputs/20_hisat2_aligned_reads
mkdir ../02_outputs/21_hisat2_alignment_summaries

# STEP 2:  ALIGN RNA-SEQ TO BOMBUS TERRESTRIS GENOME ----

# Change working directory.

cd ../00_data/00_fastq

# HISAT2 arguments: 
# -x: the basename of the index of the reference genome.
# --summary-file fileName and --new-summary: Print alignment summary in machine-friendly style to fileName.
# -U: file containing unpaired reads  

for FILE in *.fastq.gz
do
OUTFILE=${FILE:0:5}
hisat2 -x ../01_fasta/Bter_HISAT2_index --summary-file ../../02_outputs/21_hisat2_alignment_summaries/${OUTFILE}_hisat2.txt --new-summary -U ${FILE} | samtools sort -o ../../02_outputs/20_hisat2_aligned_reads/${OUTFILE}_hisat2.bam
done

# STEP 2: GENERATE MULTIQC REPORT OF ALIGNMENT SUMMARIES ----

# Change working directory.

cd ../../02_outputs/21_hisat2_alignment_summaries

# Generate report.
# Args:
# -n: rename report
# -o: save report in directory 

multiqc . -n hisat2_multiqc.html -o ../../02_multiqc_reports
