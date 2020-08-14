#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Isolate list of transcripts and corresponding genes from GTF file.
#-------------------------------------------------------------------------------
# Inputs:
# Bombus terrestris GTF file.
#
# Outputs:
# Text document linking transcript and gene IDs.
#-------------------------------------------------------------------------------

# STEP 1: UNZIP THE GTF FILE ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Change working directory to correct data subdirectory.

cd ../00_data/02_gtf

# Unzip GTF file.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout Bombus_terrestris.Bter_1.0.47.gtf.gz > Bter.gtf

# STEP 2: PRINT TRANSCRIPT AND GENE IDS FROM GTF TO A NEW FILE ----

# Using awk to search the GTF.
# If the third column equals transcript then print the twelth column (transcript id) and tenth column (gene id) separated by a tab.

awk '$3 == "transcript" {print $12"\t"$10}' Bter.gtf > transcripts2genes.txt

# STEP 3: EDIT THE FILE ----

# Remove ; from the file using sed.

sed 's/;//g' transcripts2genes.txt > transcripts2genesBter.txt

# STEP 4: MOVE THE FILE TO THE CORRECT SUBDIRECTORY ----

# Move file to appropriate folder

mv transcripts2genesBter.txt ../10_transcripts2genes
