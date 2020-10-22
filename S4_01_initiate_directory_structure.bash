#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Make a series of directories/subdirectories to conduct the analysis in.
# Move/copy files to directories.
#-------------------------------------------------------------------------------
# Inputs:
# The contents of the GitHub repository 
# https://github.com/dhcollins500/Collins-et-al_BB-M001482-1_obj1
#
# Outputs:
# A series of named directories/subdirectories containing the scripts from
# the GitHub repository. Gene lists and RNA-Seq data in the appropriate folder.
#-------------------------------------------------------------------------------

# STEP 1: SET UP THE DIRECTORY STRUCTURE BY MAKING DIRECTORIES ----

# This script assumes that the location that the GitHub repository was 
# downloaded to is the working directory.

# Make the initial directory for the analysis containing the analysis name.

mkdir S4_RNA_seq_comparison

# Change into the new directory to make further directories.

cd S4_RNA_seq_comparison

# Make directories to store the data, scripts and outputs.

mkdir 00_data
mkdir 01_scripts
mkdir 02_outputs

# Change into the data directory to make further directories.

cd 00_data

# Make directories to store different types of data.

mkdir 00_fastq
mkdir 01_fasta
mkdir 02_gtf
mkdir 03_sums
mkdir 10_transcripts2genes
mkdir 20_table_s10_csvs

# STEP 2: MOVE SCRIPTS TO CORRECT SUBDIRECTORY ----

# Change directory to the GitHub repository.

cd ../../

# Move the scripts to the 01_scripts subdirectory.

mv S4_* S4_RNA_seq_comparison/01_scripts

# STEP 3: MOVE RNA-SEQ RAW READS TO CORRECT SUBDIRECTORY ----

# "PATH/TO/RNA-SEQ/READS" needs to be CHANGED to the directory
# on your machine that contains the RNA-Seq reads from GEO GSE90751.

mv "PATH/TO/RNA-SEQ/READS"/*.fastq.gz S4_RNA_seq_comparison/00_data/00_fastq

# Remove RNA-Seq files that are not used in this analysis (09_MQ3, 11_MW1 and 16_LQ3).

rm -f S4_RNA_seq_comparison/00_data/00_fastq/09_MQ3.fastq.gz

rm -f S4_RNA_seq_comparison/00_data/00_fastq/11_MW1.fastq.gz

rm -f S4_RNA_seq_comparison/00_data/00_fastq/16_LQ3.fastq.gz

# STEP 4: MOVE CURRENT STUDY GENE LISTS TO CORRECT SUBDIRECTORY ----

# Move the gene lists (Table_S10_xx.csv).

mv Table_S10*.csv S4_RNA_seq_comparison/00_data/20_table_s10_csvs
