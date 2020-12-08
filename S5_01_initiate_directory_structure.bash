#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Make a series of directories/subdirectories to conduct the analysis in.
# Move/copy files to directories.
#-------------------------------------------------------------------------------
# Inputs:
# The contents of the GitHub repository 
# https://github.com/dhcollins500/Collins-et-al_BB-M001482-1_obj1
#
# Outputs:
# A series of named directories/subdirectories containing the scripts from
# the GitHub repository.
#-------------------------------------------------------------------------------

# STEP 1: SET UP THE DIRECTORY STRUCTURE BY MAKING DIRECTORIES ----

# This script assumes that the location that the GitHub repository was 
# downloaded to is the working directory.

# Make Overall Analysis Directory ----

# Make the initial directory for the analysis containing the analysis name.
# Amel = Apis mellifera

mkdir S5_analyses

# Make Subdirectories ----

# Change into the new directory to make further directories.

cd S5_analyses

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
mkdir 10_SRA_Acc_Lists
mkdir 20_table_s7_csv

# STEP 2: MOVE SCRIPTS TO CORRECT SUBDIRECTORY ----

# Change Directory to the GitHub Repository ----

cd ../../

# Move Scripts to 01_scripts Subdirectory ----

mv S5_* S5_analyses/01_scripts

# STEP 3: MOVE SRA ACCESSION LISTS TO CORRECT SUBDIRECTORY ----

# Change working directory.

cd S5_analyses/01_scripts

# Move .txt Files to 00_data/20_SRA_Acc_Lists Subdirectory ----

mv *AccList.txt ../00_data/10_SRA_Acc_Lists

# STEP 4: MOVE BOMBUS TO APIS GENE LIST TO CORRECT SUBDIRECTORY ----

# Move the gene list (Table_S7.csv).

mv ../../Table_S7.csv ../00_data/20_table_s7_csv
