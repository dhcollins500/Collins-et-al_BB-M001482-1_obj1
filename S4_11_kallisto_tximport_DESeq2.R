#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Summarise transcript-level counts to gene level and create DESeq2
# data set.
#-------------------------------------------------------------------------------
# Inputs:
# Abundances from Kallisto pseudoalignment. Table describing the experimental
# design.

# Outputs:
# dds (DESeq2 data set) object called ddsKallisto.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(tximport)  # tximport()
library(DESeq2)  # DESeqDataSetFromTximport(), DESeq() 

# LOADING DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Transcript to Gene File ----

setwd("../00_data/10_transcripts2genes")

t2g <- read.table(file = "transcripts2genesBter.txt",
                    header = FALSE,
                    col.names = c("TXNAME",
                                  "GENEID"))

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Set Internal Variables ----

baseDir <- "ABSOLUTE/PATH/TO/FILES/02_outputs/10_kallisto_pseudoalignments"
# This string is specific to your computer, PLEASE CHANGE accordingly.

abundanceDir <- "abundances"

summaryDir <- "summaries"

# Generate Character Vector with Sample Names ----

sample_id <- 
  read.table(file.path(baseDir, summaryDir, 
                    "caste_study_design.txt"),
             header = TRUE)

# Add Column Containing Full Sample ID, to Match the Folder Names ----

sample_id$full_id <- paste0(dir(file.path(baseDir, abundanceDir)))

# Add Condition Column to sample_id ----

sample_id$condition <- paste0(sample_id$time_point, sample_id$caste) 

# Generate Named Vector with Path to Quantification Files ----

files <- file.path(baseDir, abundanceDir, sample_id$full_id, "abundance.h5")

names(files) <- paste0(sample_id$sample)

# Import Transcript-level Estimates and Summarise to Gene Level
# Producing "Original Counts and Offsets" ----

txi <- tximport(files, type = "kallisto", tx2gene = t2g)

# Initiate DESeq2 Data set (dds) Object ----

ddsKallisto <- DESeqDataSetFromTximport(txi, sample_id, ~condition)

# Filter Rows ----

# Filter rows of dds object to remove any genes with less than 10 counts across 
# all samples.

ddsKallisto <- ddsKallisto[rowSums(counts(ddsKallisto)) >= 10, ]

# Conduct Differential Expression Analysis ----

ddsKallisto <- DESeq(ddsKallisto)

# Remove Objects that are no Longer Required ----

rm(sample_id, t2g, abundanceDir, baseDir, files, summaryDir, txi)
