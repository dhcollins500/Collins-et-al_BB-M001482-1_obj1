#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Load HTSeq gene counts and create DESeq2 data set.
#-------------------------------------------------------------------------------
# Inputs:
# HTSeq-count count files. 

# Outputs:
# dds (DESeq2 data set) object called ddsHTSeq.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(DESeq2)  # DESeqDataSetFromHTSeqCount(), DESeq()

# LOADING DATA ----

# No data loaded.

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Prepare Variables ----

# Specify variable that points to directory where HTSeq-count results files are.

directory <- "PATH/TO/FILES/02_outputs/30_htseq_count"
# This string is specific to your computer, PLEASE CHANGE accordingly.

# Specify sample files to read in.

sampleFiles <- list.files(directory)

# Specify the condition of each file. 

sampleCondition <- substr(sampleFiles, start = 4, stop = 5)

# Make a sample table with condition as a factor.

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

sampleTable$condition <- factor(sampleTable$condition)

# Initiate the DESeqDataSet ----

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

# Filter Rows ----

# Filter rows of dds object to remove any genes with less than 10 counts across 
# all samples.

ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 10, ]

# Conduct Differential Expression Analysis ----

ddsHTSeq <- DESeq(ddsHTSeq)

# Remove Objects that are no Longer Required ----

rm(directory, sampleFiles, sampleCondition, sampleTable)
