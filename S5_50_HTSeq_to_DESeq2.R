#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Load HTSeq gene counts and create DESeq2 data set.
#-------------------------------------------------------------------------------
# Inputs:
# HTSeq-count count files. 

# Outputs:
# dds (DESeq2 data set) objects called ddsCameron, ddsHe and ddsAmsalem.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(DESeq2)  # DESeqDataSetFromHTSeqCount(), DESeq()

# LOADING DATA ----

# No data loaded.

# FUNCTION DEFINITIONS ----
# Note: File paths in the below function are computer specific and need
#       to be changed for your computer.

LoadHTSeqFiles <- function (x) {
  # Loads HTSeq-count files for Cameron et al. (2013), He et al. (2017) or 
  # Amsalem et al. (2015), and generates a dds object filtered to remove lowly 
  # expressed genes.
  #
  # Args:
  #   x: string denoting whether the data loaded should be from Cameron
  #      et al. (2013), He et al. (2017) or Amsalem et al. (2015) ("cameron",
  #      "he" or "amsalem").
  #
  # Returns:
  #   A dds object filtered to remove lowly expressed genes.
  
  # Set variables based on arguments.
  
  if (x == "cameron") {
    resultsFolder <- "30_cameron_htseq_count"
    y <- 6  # used to specify sample name from sample file name.
    z <- 11 # used to specify sample name from sample file name.
    a <- 6  # used to specify condition from sample file name.
  } else if (x == "he") {
    resultsFolder <- "31_he_htseq_count"
    y <- 1  # used to specify sample name from sample file name.
    z <- 9  # used to specify sample name from sample file name.
    a <- 4  # used to specify condition from sample file name.
  } else if (x == "amsalem") {
    resultsFolder <- "32_amsalem_htseq_count"
    y <- 1  # used to specify sample name from sample file name.
    z <- 9  # used to specify sample name from sample file name.
    a <- 1  # used to specify condition from sample file name.
  } else {
    stop('Argument x must equal "cameron", "he" or "amsalem".')
  }
  
  # Specify variable that points to directory where HTSeq-count results files 
  # are.
  
  directory <- paste0("ABSOLUTE/PATH/TO/FILES/S5_analyses/02_outputs/",
                      resultsFolder)
  # The first string is specific to your computer, PLEASE CHANGE accordingly.
  
  # Specify sample files to read in.
  
  sampleFiles <- list.files(directory)
  
  # Specify the names of each sample.
  
  sampleNames <- substr(sampleFiles, start = y, stop = z)
  
  # Specify the condition of each file. 
  
  sampleCondition <- substr(sampleFiles, start = y, stop = a)
  
  # Make a sample table with condition as a factor.
  
  sampleTable <- data.frame(sampleName = sampleNames,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  
  sampleTable$condition <- factor(sampleTable$condition)
  
  # Initiate the DESeqDataSet
  
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                    directory = directory,
                                    design= ~ condition)
  
  # Filter rows
  
  # Filter rows of dds object to remove any genes with less than 10 counts across 
  # all samples.
  
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Conduct differential expression analysis.
  
  dds <- DESeq(dds)
  
  # Return dds.
  
  return(dds)
  
}

# EXECUTED STATEMENTS ----

# Cameron et al. (2013) Data ----

ddsCameron <- LoadHTSeqFiles("cameron")

# He et al. (2017) Data ----

ddsHe <- LoadHTSeqFiles("he")

# Amsalem et al. (2015) Data ----

ddsAmsalem <- LoadHTSeqFiles("amsalem")
