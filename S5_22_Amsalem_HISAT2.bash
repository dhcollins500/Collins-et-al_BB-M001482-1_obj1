#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 3, 4, 5
# Tasks: Align mRNA-seq reads from Amsalem et al. (2015) to Bombus terrestris
# (Bter) genome.

#-------------------------------------------------------------------------------
# Inputs:
# mRNA-seq files from Amsalem et al. (2015). HISAT2 index for Bter genome.
#
# Outputs:
# Folders with results of alignments.
# Multiqc summary of alignments.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULE ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load HISAT2/2.1.0
module load python/anaconda/2019.10/3.7  # where MultiQC is installed.
module load samtools/1.10

# STEP 1: MAKE OUTPUT DIRECTORIES FOR RESULTS ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

mkdir ../02_outputs/24_amsalem_hisat2_aligned_reads
mkdir ../02_outputs/25_amsalem_hisat2_alignment_summaries

# STEP 2: ALIGNMENT WITH HISAT2 ----

# Change directory

cd ../00_data/00_fastq

# HISAT2 arguments: 
# -x: the basename of the index of the reference genome.
# --summary-file fileName and --new-summary: Print alignment summary in machine-friendly style to fileName.
# -U: file containing unpaired reads  

for FILE in amsalem*.fastq 
do
OUTPUT=${FILE:13:9}
hisat2 -x ../01_fasta/Bter_HISAT2_index -p 6 --summary-file ../../02_outputs/25_amsalem_hisat2_alignment_summaries/${OUTPUT}_hisat2.txt --new-summary -U ${FILE} | samtools sort -o ../../02_outputs/24_amsalem_hisat2_aligned_reads/${OUTPUT}_hisat2.bam
done

# STEP 3: INDEX ALL .BAM FILES ----

# Change working directory.

cd ../../02_outputs/24_amsalem_hisat2_aligned_reads

# Looping through all bam files in the directory and generating an index for them.

for BAMFILE in *hisat2.bam
do
samtools index ${BAMFILE}
done

# STEP 4: GENERATE MULTIQC REPORT OF ALIGNMENT SUMMARIES ----

# Change working directory.

cd ../25_amsalem_hisat2_alignment_summaries

# Generate report.
# Args:
# -n: rename report
# -o: save report in directory 

multiqc . -n amsalem_hisat2_genome_multiqc_report.html -o ../10_multiqc_reports
