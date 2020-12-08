#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 3, 4, 5
# Tasks: Download the raw RNA-Seq files for Amsalem et al. (2015).
#-------------------------------------------------------------------------------
# Inputs:
# List of SRA accession numbers.
#
# Outputs:
# Renamed mRNA-seq files for Amsalem et al. (2015) in the 00_data/00_fastq 
# directory.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module load sra/sra-2.10.2

# STEP 1: DOWNLOAD AND RENAME AMSALEM ET AL. (2015) MRNA-SEQ DATA ----

# Change to 00_fastq Directory ----

cd ../00_data/00_fastq

# Configure SRA Toolkit ----

# On first use, the SRA toolkit has to be configured interactively.
# Run: vdb-config --interactive then press X (see https://github.com/ncbi/sra-tools/issues/291)

# Download Data ----

# Prefetching the Amsalem et al. (2015) run using the provided list of accession numbers.
# Note: data are 100 bp, single end data.

prefetch --option-file ../10_SRA_Acc_Lists/S5_sraAmsalemAccList.txt

# Check the integrity of each downloaded file, then convert the .sra files into 
# fastq format using a loop.
# Args:
#–split-files separates the read into left and right ends, and puts the forward 
# and reverse reads in two separate files.

for FILE in SRR2*
do
cd ${FILE}
vdb-validate ${FILE}.sra
fasterq-dump --split-files ${FILE}.sra
cd ../
done

# Note: check output of the above loop manually - all parameters should be "ok" or
# "consistent".

# Loop over folders to move the .fastq file

for FILE in SRR2*
do
cd ${FILE}
mv ${FILE}.sra.fastq ../${FILE}.sra.fastq
cd ../
done

# Rename Data ---
# Format of name: study_condition_geneticline.fastq
# Where conditions are: F = founding post diapause, D = diapause, M = mated

mv SRR2396657.sra.fastq amsalem_2015_F_line191.fastq
mv SRR2396656.sra.fastq amsalem_2015_F_line051.fastq
mv SRR2396655.sra.fastq amsalem_2015_F_line296.fastq
mv SRR2396654.sra.fastq amsalem_2015_F_line141.fastq
mv SRR2396653.sra.fastq amsalem_2015_F_line207.fastq
mv SRR2396652.sra.fastq amsalem_2015_D_line191.fastq
mv SRR2396651.sra.fastq amsalem_2015_D_line207.fastq
mv SRR2396650.sra.fastq amsalem_2015_D_line296.fastq
mv SRR2396649.sra.fastq amsalem_2015_D_line051.fastq
mv SRR2396648.sra.fastq amsalem_2015_D_line141.fastq
mv SRR2396644.sra.fastq amsalem_2015_M_line426.fastq
mv SRR2396641.sra.fastq amsalem_2015_M_line141.fastq
mv SRR2396640.sra.fastq amsalem_2015_M_line191.fastq
mv SRR2396634.sra.fastq amsalem_2015_M_line207.fastq
mv SRR2396633.sra.fastq amsalem_2015_M_line296.fastq

# STEP 2: REMOVE TEMPORARY SRA DIRECTORIES AND FILES

# Remove Directories and Files ----
# Args:
# -f: never prompt 
# -r: recursively remove directories and their contents

rm -f -r SRR*
