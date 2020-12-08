#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Generate indexes for genome related files from Apis mellifera (Amel) and
# Bombus terrestris (Bter).
#-------------------------------------------------------------------------------
# Inputs:
# Amel genome fasta file and GTF file.
# Bter genome fasta file and GTF file.
#
# Outputs:
# HISAT2 indexes for Amel and Bter fasta file.
# ALFA indexes for Amel and Bter GTF file.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULE ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load HISAT2/2.1.0
module load python/anaconda/2019.10/3.7  # Where ALFA is installed.
module load bedtools/2.29.2  # Needed for ALFA.

# STEP 1: GENERATE HISAT2 INDEX OF AMEL GENOME ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Change working directory.

cd ../00_data/01_fasta

# Indexing command.

hisat2-build Amel_v4_5_genome_DNA.fasta Amel_genome_HISAT2_index

# STEP 2: GENERATE HISAT2 INDEX OF BTER GENOME ----

# Indexing command.

hisat2-build Bter_toplevel_dna.fasta Bter_HISAT2_index

# STEP 3: GENERATE ALFA INDEX OF AMEL GTF FILE ---

# Change working directory.

cd ../02_gtf

# Index the GTF file with alfa.
# Args:
# -a: path to annotation file.
# -g: alfa index basename.

alfa -a Amel_v4_5_48.gtf -g Amel_v4_5_48_gtf_alfa_index

# STEP 4: GENERATE ALFA INDEX OF GTF FILE ---

# Index the GTF file with alfa.

alfa -a Bter.gtf -g Bter_gtf_alfa_index
