#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Pseudoalign mRNA-seq reads to Bombus terrestris transcriptome using
# Kallisto.
#-------------------------------------------------------------------------------
# Inputs:
# Kallisto index of Bombus terrestris transcriptome, raw mRNA-seq read files.
#
# Outputs:
# Abundances of transcripts in each mRNA-seq file and summaries of the
# pseudoalignments. MultiQC report of Kallisto pseudoalignment statistics.
# Text file with study design for tximport.
#-------------------------------------------------------------------------------

# LOAD MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add kallisto/0.46.1
module add python/anaconda/2019.10/3.7  # where MultiQC is installed.

# STEP 1: MAKE OUTPUT DIRECTORIES FOR RESULTS ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

mkdir ../02_outputs/02_multiqc_reports
mkdir ../02_outputs/10_kallisto_pseudoalignments
mkdir ../02_outputs/10_kallisto_pseudoalignments/abundances
mkdir ../02_outputs/10_kallisto_pseudoalignments/summaries 

# STEP 2: PSEUDOALIGNMENT WITH KALLISTO ----

# Change directory

cd ../00_data/00_fastq

# Psuedo align reads to transcriptome
# Args:
# -i: fileName for index to be used for quanitification.
# -o: directory to write output to.
# --single: quantify single-end reads.
# -l: estimated average fragment length.
# -s: estimated standard deviation of fragment length. 

for FILE in *fastq.gz
do
OUTPUT=${FILE:0:5}
kallisto quant -i ../01_fasta/Bter_v1_cDNA_index.idx -o ../../02_outputs/10_kallisto_pseudoalignments/abundances/${OUTPUT} --single -l 200 -s 20 ${FILE} 2> ../../02_outputs/10_kallisto_pseudoalignments/summaries/${OUTPUT}.txt 
done

# STEP 3: GENERATE MULTIQC REPORT OF PSEUDOALIGNMENT SUMMARIES ----

# Change working directory.

cd ../../02_outputs/10_kallisto_pseudoalignments/summaries

# Generate report.
# Args:
# -n: rename report
# -o: save report in directory 

multiqc . -n kallisto_multiqc.html -o ../../02_multiqc_reports

# STEP 4: GENERATE STUDY DESIGN FILE ----

# Initiate the file.
# Args:
# -e: enable interpretation of backslash escapes.

echo -e "sample\ttime_point\tcaste" > caste_study_design.txt

# Loop over file names to fill in the details.

for NAME in ?????.txt
do
	SAMPLE=${NAME:0:5}
	TIMEPOINT=${NAME:3:1}
	CASTE=${NAME:4:1}
echo -e "${SAMPLE}\t${TIMEPOINT}\t${CASTE}" >> caste_study_design.txt
done
