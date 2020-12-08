#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Generate and combine FastQC reports on mRNA-seq files from each of Cameron et 
# al. (2013), He et al. (2017) and Amsalem et al. (2015).
#-------------------------------------------------------------------------------
# Inputs:
# mRNA-seq files from Cameron et al. (2013), He et al. (2017) and 
# Amsalem et al. (2015).
#
# Outputs:
# MultiQC report combining the FastQC reports from the mRNA-seq reads for
# Cameron et al. (2013), He et al. (2017) and Amsalem et al. (2015).
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULE ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load fastqc/0.11.9
module load python/anaconda/2019.10/3.7  # Where MultiQC module is installed.

# STEP 1: RUNNING FASTQC ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Change directory to the location of the data.

cd ../00_data/00_fastq

# Run fastQC program on all mRNA-seq files.

fastqc *fastq

# STEP 2: MAKING DIRECTORIES FOR RESULTS ----

# Making output directories.

mkdir ../../02_outputs/01_cameron_fastqc
mkdir ../../02_outputs/02_he_fastqc
mkdir ../../02_outputs/03_amsalem_fastqc

# STEP 3: PARSING FASTQC RESULTS ----

# Moving Results ----

# Move all fastqc output files (.html and .zip) from 
# Cameron et al. (2013) to the appropriate output folder.

mv cameron*.zip ../../02_outputs/01_cameron_fastqc

mv cameron*.html ../../02_outputs/01_cameron_fastqc

# Move all fastqc output files (.html and .zip) from 
# He et al. (2017) to the appropriate output folder.

mv he*.zip ../../02_outputs/02_he_fastqc

mv he*.html ../../02_outputs/02_he_fastqc

# Move all fastqc output files (.html and .zip) from 
# Amsalem et al. (2015) to the appropriate output folder.

mv amsalem*.zip ../../02_outputs/03_amsalem_fastqc

mv amsalem*.html ../../02_outputs/03_amsalem_fastqc

# Combining FastQC Data into MultiQC Reports ----

# Make output directory for MultiQC reports.

mkdir ../../02_outputs/10_multiqc_reports

# Change working directory to the fastqc data for Cameron et al. (2013).

cd ../../02_outputs/01_cameron_fastqc

# Run multiqc python module to collate all fastqc reports in folder together into one report.
# -n fileName: call report fileName rather than default of multiqc_report.html.
# -o file/path: create report in subdirectory file/path.

multiqc . -n cameron_mRNA_seq_reads_fastqc_multiqc_report.html -o ../10_multiqc_reports

# Change working directory to the fastqc data for He et al. (2017).

cd ../02_he_fastqc

# Run multiqc python module to collate all fastqc reports in folder together into one report.

multiqc . -n he_mRNA_seq_reads_fastqc_multiqc_report.html -o ../10_multiqc_reports

# Change working directory to the fastqc data for Amsalem et al. (2015).

cd ../03_amsalem_fastqc

# Run multiqc python module to collate all fastqc reports in folder together into one report.

multiqc . -n amsalem_mRNA_seq_reads_fastqc_multiqc_report.html -o ../10_multiqc_reports
