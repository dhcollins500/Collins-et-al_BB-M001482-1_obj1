#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Index transcriptome and genome with appropriate software.
#-------------------------------------------------------------------------------
# Inputs:
# Bombus terrestris transcriptome and genome.
#
# Outputs:
# Kallisto index for Bombus terrestris transcriptome and HISAT2 index for Bombus 
# terrestris genome.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add kallisto/0.46.1
module add HISAT2/2.1.0

# STEP 1: UNZIP THE GENOME AND TRANSCRIPTOME FILES ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Change working directory.

cd ../00_data/01_fasta

# Unzip genome fasta file.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz > Bter_toplevel_dna.fasta

# Unzip transcriptome fasta file.

gunzip --stdout Bombus_terrestris.Bter_1.0.cdna.all.fa.gz > Bter_cDNA.fasta

# STEP 2: GENERATE KALLISTO INDEX OF TRANSCRIPTOME ----

# Indexing command.
# Args:
# -i: fileName for index to be created.

kallisto index -i Bter_v1_cDNA_index.idx Bter_cDNA.fasta

# STEP 3: GENERATE HISAT2 INDEX OF GENOME ----

# Indexing command.

hisat2-build Bter_toplevel_dna.fasta Bter_HISAT2_index
