#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1, 2, 4 
# Tasks: Download the raw mRNA-seq files for Cameron et al. (2013).
#-------------------------------------------------------------------------------
# Inputs:
# None.
#
# Outputs:
# Renamed mRNA-seq files for Cameron et al. (2013) in the 00_data/00_fastq 
# directory.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module load sra/sra-2.10.2

# STEP 1: DOWNLOAD AND RENAME CAMERON ET AL. (2013) MRNA-SEQ DATA ----

# Change to 00_fastq Directory ----

cd ../00_data/00_fastq

# Configure SRA Toolkit ----

# On first use, the SRA toolkit has to be configured interactively.
# Run: vdb-config --interactive then press X (see https://github.com/ncbi/sra-tools/issues/291)

# Download Data ----

# Prefetching the Cameron et al. (2013) run using the provided list of accession numbers.

prefetch --option-file ../10_SRA_Acc_Lists/S5_sraCameronAccList.txt

# Check the integrity of each downloaded file, then convert the .sra files into 
# fastq format using a loop.
# Args:
#–split-files separates the read into left and right ends, and puts the forward 
# and reverse reads in two separate files.

for FILE in SRR1*
do
cd ${FILE}
vdb-validate ${FILE}.sra
fasterq-dump --split-files ${FILE}.sra
cd ../
done

# Note: check output of the above loop manually - all parameters should be "ok" or
# "consistent".

# Loop over folders to move the .fastq file

for FILE in SRR1*
do
cd ${FILE}
mv ${FILE}.sra.fastq ../${FILE}.sra.fastq
cd ../
done

# Rename Data ---

# Rename .fastq files with more informative names.
# Name format = study_timepoint_phenotype_rep.fastq
# where phenotype Q = queen and phenotype W = worker.

mv SRR1028781.sra.fastq cameron_2013_60hr_Q_rep1.fastq
mv SRR1028782.sra.fastq cameron_2013_60hr_Q_rep2.fastq
mv SRR1028783.sra.fastq cameron_2013_60hr_W_rep1.fastq
mv SRR1028784.sra.fastq cameron_2013_60hr_W_rep2.fastq
