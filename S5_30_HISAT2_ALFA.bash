#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: ALFA QC of HISAT2 alignment of Cameron et al. (2013), He et al. (2017) 
# to the Apis mellifera (Amel) genome, and Amsalem et al. (2015) alignment to the
# Bombus terrestris (Bter) genome.
#-------------------------------------------------------------------------------
# Inputs:
# .BAM alignments of Cameron et al. (2013) and He et al. (2017) to Amel genome,
# and Amsalem et al. (2015) to Bter genome. 
# ALFA index of Amel and Bter GTF files.
#
# Outputs:
# Two PDF files of ALFA results for each of the three data sets.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULE ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module load python/anaconda/2019.10/3.7  # Where ALFA is installed.
module load bedtools/2.29.2  # Needed for ALFA.

# STEP 1: ALFA ANALYSIS OF CAMERON ET AL. (2013) ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Change directory.

cd ../02_outputs/20_cameron_hisat2_aligned_reads

# Generate list.

# Loop to create $LIST variable with list of file and names for input files for alfa.

for FILE in *.bam
do
NAME=${FILE:0:11}  # variable:starting_position:length
LIST="${LIST:+$LIST } $FILE $NAME"
done

# Run Alfa ----

# Python module produces graph showing read distribution across genome features.

# alfa arguments:
# -g file/Path/fileName: path and basename of alfa index
# --bam fileName label: file name and associated label (for the plot legends) for the bam files.
# --pdf file/Path/fileName: path to save the plots to in pdf format.

alfa -g ../../00_data/02_gtf/Amel_v4_5_48_gtf_alfa_index --bam ${LIST} --pdf ../10_multiqc_reports/cameron_ALFA_results 

# STEP 2: ALFA ANALYSIS OF HE ET AL. (2017) ----

# Change directory.

cd ../22_he_hisat2_aligned_reads

# Generate list.

# Clear LIST of previous values

LIST=

# Loop to create $LIST variable with list of file and names for input files for alfa.

for FILE in *.bam
do
NAME=${FILE:0:9}  # variable:starting_position:length
LIST="${LIST:+$LIST } $FILE $NAME"
done

# Run Alfa ----

# See above for arguments.

alfa -g ../../00_data/02_gtf/Amel_v4_5_48_gtf_alfa_index --bam ${LIST} --pdf ../10_multiqc_reports/he_ALFA_results 

# STEP 3: ALFA ANALYSIS OF AMSALEM ET AL. (2015) ----

# Change directory.

cd ../24_amsalem_hisat2_aligned_reads

# Generate list.

# Clear LIST of previous values

LIST=

# Loop to create $LIST variable with list of file and names for input files for alfa.

for FILE in *.bam
do
NAME=${FILE:0:9}  # variable:starting_position:length
LIST="${LIST:+$LIST } $FILE $NAME"
done

# Run Alfa ----

# See above for arguments.

alfa -g ../../00_data/02_gtf/Bter_gtf_alfa_index --bam ${LIST} --pdf ../10_multiqc_reports/amsalem_ALFA_results 
