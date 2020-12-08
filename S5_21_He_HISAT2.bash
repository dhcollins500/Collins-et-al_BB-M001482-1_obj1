#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1, 2, 4
# Tasks: Align mRNA-seq reads from He et al. (2017) to Apis mellifera (Amel) 
# genome.
#-------------------------------------------------------------------------------
# Inputs:
# mRNA-seq files from He et al. (2017). HISAT2 index for Amel genome.
#
# Outputs:
# Folders with results of alignments.
# Multiqc summary of alignments.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULE ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load HISAT2/2.1.0
module load python/anaconda/2019.10/3.7  # where Multiqc is installed.
module load samtools/1.10

# STEP 1: MAKE OUTPUT DIRECTORIES FOR RESULTS ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

mkdir ../02_outputs/22_he_hisat2_aligned_reads
mkdir ../02_outputs/23_he_hisat2_alignment_summaries

# STEP 2: ALIGNMENT WITH HISAT2 ----

# Change directory

cd ../00_data/00_fastq

# HISAT2 arguments: 
# -x: the basename of the index of the reference genome.
# --summary-file fileName and --new-summary: Print alignment summary in machine-friendly style to fileName.
# --dta: Tailor alignment report for Downstream Transcriptome Assembly, in this case with StringTie.
# -1: name of file containing read 1.
# -2: name of file containing read 2.

for FILE in he*R1.fastq 
do
OUTPUT=${FILE:8:9}
FILEHANDLE=${FILE%?????????}  # Remove last 9 characters of the string.
hisat2 -x ../01_fasta/Amel_genome_HISAT2_index --summary-file ../../02_outputs/23_he_hisat2_alignment_summaries/${OUTPUT}_hisat2.txt --new-summary --dta -1 ${FILE} -2 ${FILEHANDLE}_R2.fastq | samtools sort -o ../../02_outputs/22_he_hisat2_aligned_reads/${OUTPUT}_hisat2.bam
done

# STEP 3: INDEX ALL .BAM FILES ----

# Change working directory.

cd ../../02_outputs/22_he_hisat2_aligned_reads

# Looping through all bam files in the directory and generating an index for them.

for BAMFILE in *hisat2.bam
do
samtools index ${BAMFILE}
done

# STEP 4: GENERATE MULTIQC REPORT OF ALIGNMENT SUMMARIES ----

# Change working directory.

cd ../23_he_hisat2_alignment_summaries

# Generate report.
# Args:
# -n: rename report
# -o: save report in directory 

multiqc . -n he_hisat2_genome_multiqc_report.html -o ../10_multiqc_reports
