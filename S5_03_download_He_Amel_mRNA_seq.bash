#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1, 2, 4
# Tasks: Download the raw mRNA-seq files for He et al. (2017).
#-------------------------------------------------------------------------------
# Inputs:
# None.
#
# Outputs:
# Renamed mRNA-seq files for He et al. (2017) in the 00_data/00_fastq directory.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module load sra/sra-2.10.2

# STEP 1: DOWNLOAD AND RENAME HE ET AL. (2017) MRNA-SEQ DATA ----

# Change to 00_fastq Directory ----

cd ../00_data/00_fastq

# Configure SRA Toolkit ----

# On first use, the SRA toolkit has to be configured interactively.
# Run: vdb-config --interactive then press X (see https://github.com/ncbi/sra-tools/issues/291)

# Download Data ----

# Prefetching the He et al. (2017) run using the provided list of accession numbers.

prefetch --option-file ../10_SRA_Acc_Lists/S5_sraHeAccList.txt

# Check the integrity of each downloaded file, then convert the .sra files into 
# fastq format using a loop.
# Args:
#–split-files separates the read into left and right ends, and puts the forward 
# and reverse reads in two separate files.

for FILE in SRR3*
do
cd ${FILE}
vdb-validate ${FILE}.sra
fasterq-dump --split-files ${FILE}.sra
cd ../
done

# Note: check output of the above loop manually - all parameters should be "ok" or
# "consistent".

# Loop over folders to move the .fastq file

for FILE in SRR3*
do
cd ${FILE}
mv ${FILE}.sra_1.fastq ../${FILE}.sra_1.fastq
mv ${FILE}.sra_2.fastq ../${FILE}.sra_2.fastq
cd ../
done

# Rename Data ---

# Rename .fastq files with more informative names.
# Name format = study_timepoint_phenotype_rep_readnumber.fastq
# where phenotype Q = queen and phenotype W = worker.

mv SRR3102934.sra_1.fastq he_2017_2d_W_rep2_R1.fastq
mv SRR3102934.sra_2.fastq he_2017_2d_W_rep2_R2.fastq

mv SRR3123272.sra_1.fastq he_2017_2d_W_rep1_R1.fastq
mv SRR3123272.sra_2.fastq he_2017_2d_W_rep1_R2.fastq

mv SRR3123273.sra_1.fastq he_2017_2d_W_rep6_R1.fastq
mv SRR3123273.sra_2.fastq he_2017_2d_W_rep6_R2.fastq

mv SRR3123275.sra_1.fastq he_2017_2d_W_rep5_R1.fastq
mv SRR3123275.sra_2.fastq he_2017_2d_W_rep5_R2.fastq

mv SRR3123276.sra_1.fastq he_2017_2d_W_rep4_R1.fastq
mv SRR3123276.sra_2.fastq he_2017_2d_W_rep4_R2.fastq

mv SRR3123277.sra_1.fastq he_2017_2d_W_rep3_R1.fastq
mv SRR3123277.sra_2.fastq he_2017_2d_W_rep3_R2.fastq

mv SRR3123279.sra_1.fastq he_2017_4d_W_rep1_R1.fastq
mv SRR3123279.sra_2.fastq he_2017_4d_W_rep1_R2.fastq

mv SRR3123281.sra_1.fastq he_2017_4d_W_rep2_R1.fastq
mv SRR3123281.sra_2.fastq he_2017_4d_W_rep2_R2.fastq

mv SRR3123337.sra_1.fastq he_2017_4d_W_rep3_R1.fastq
mv SRR3123337.sra_2.fastq he_2017_4d_W_rep3_R2.fastq

mv SRR3123340.sra_1.fastq he_2017_4d_W_rep4_R1.fastq
mv SRR3123340.sra_2.fastq he_2017_4d_W_rep4_R2.fastq

mv SRR3123341.sra_1.fastq he_2017_4d_W_rep5_R1.fastq
mv SRR3123341.sra_2.fastq he_2017_4d_W_rep5_R2.fastq

mv SRR3123342.sra_1.fastq he_2017_4d_W_rep6_R1.fastq
mv SRR3123342.sra_2.fastq he_2017_4d_W_rep6_R2.fastq

mv SRR3123355.sra_1.fastq he_2017_2d_Q_rep1_R1.fastq
mv SRR3123355.sra_2.fastq he_2017_2d_Q_rep1_R2.fastq

mv SRR3123357.sra_1.fastq he_2017_2d_Q_rep2_R1.fastq
mv SRR3123357.sra_2.fastq he_2017_2d_Q_rep2_R2.fastq

mv SRR3123359.sra_1.fastq he_2017_2d_Q_rep3_R1.fastq
mv SRR3123359.sra_2.fastq he_2017_2d_Q_rep3_R2.fastq

mv SRR3123361.sra_1.fastq he_2017_2d_Q_rep4_R1.fastq
mv SRR3123361.sra_2.fastq he_2017_2d_Q_rep4_R2.fastq

mv SRR3123362.sra_1.fastq he_2017_2d_Q_rep5_R1.fastq
mv SRR3123362.sra_2.fastq he_2017_2d_Q_rep5_R2.fastq

mv SRR3123364.sra_1.fastq he_2017_2d_Q_rep6_R1.fastq
mv SRR3123364.sra_2.fastq he_2017_2d_Q_rep6_R2.fastq

mv SRR3123372.sra_1.fastq he_2017_4d_Q_rep2_R1.fastq
mv SRR3123372.sra_2.fastq he_2017_4d_Q_rep2_R2.fastq

mv SRR3123380.sra_1.fastq he_2017_4d_Q_rep1_R1.fastq
mv SRR3123380.sra_2.fastq he_2017_4d_Q_rep1_R2.fastq

mv SRR3123385.sra_1.fastq he_2017_4d_Q_rep6_R1.fastq
mv SRR3123385.sra_2.fastq he_2017_4d_Q_rep6_R2.fastq

mv SRR3123388.sra_1.fastq he_2017_4d_Q_rep5_R1.fastq 
mv SRR3123388.sra_2.fastq he_2017_4d_Q_rep5_R2.fastq

mv SRR3123389.sra_1.fastq he_2017_4d_Q_rep4_R1.fastq
mv SRR3123389.sra_2.fastq he_2017_4d_Q_rep4_R2.fastq

mv SRR3123390.sra_1.fastq he_2017_4d_Q_rep3_R1.fastq
mv SRR3123390.sra_2.fastq he_2017_4d_Q_rep3_R2.fastq

# STEP 2: REMOVE TEMPORARY SRA DIRECTORIES AND FILES

# Remove Directories and Files ----
# Args:
# -f: never prompt 
# -r: recursively remove directories and their contents

rm -f -r SRR*
