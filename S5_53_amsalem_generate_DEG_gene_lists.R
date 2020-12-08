#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 3, 4 and 5
# Tasks: Produce lists of differentially expressed genes (DEGs) from 
# Amsalem et al. (2015).
#-------------------------------------------------------------------------------
# Inputs:
# ddsAmsalem object.

# Outputs:
# .csv files containing lists of DEGs for Amsalem et al 2015.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# DESeq2 loaded in the script in the "Data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Amsalem et al. (2015) Data ----

source("S5_50_htseq_to_DESeq2.R")
sessionInfo()

# FUNCTION DEFINITIONS ----

MakeGeneList <- function (x, y = 0) {
  # Produces and saves a list of DEGs that meet a specific log-fold change (LFC)
  # for a given time point, caste and study.
  #
  # Args:
  #   x: string denoting which type of queen from Amsalem et al. (2015) should
  #      be compared to diapause queens (mated = "M" or founder post 
  #      diapause = "F"). 
  #   y: numeric denoting the LFC to use when filter DESeq2 results
  #      (default of 0 (i.e. all significant genes)).
  #
  # Returns:
  #   A .csv file of Bombus terrestris gene names for the significant 
  #   genes that meet the LFC threshold set.
  
  # Note: false discovery rate (FDR) of 0.05 used, but this can be changed by
  # changing the alpha level in the results() call below.
  
  # Extract results.
  
  allExtractedResults <- results(dds, contrast = c("condition", "D", x),
                                 alpha = 0.05, lfcThreshold = y)
  
  # Remove NAs from padj column (represent genes excluded from analysis as all 
  # counts were 0 or it contained an extreme count outlier).
  
  allExtractedResults <- allExtractedResults[!is.na(allExtractedResults$padj), ]
  
  # Split results into genes more expressed in D, and more expressed in "x".
  
  dGenes <- 
    allExtractedResults[allExtractedResults[, "log2FoldChange"] > 0, ]

  xGenes <- 
    allExtractedResults[allExtractedResults[, "log2FoldChange"] < 0, ]

  # Index only significant (padj < 0.05) genes.
  
  sigDGenes <- as.data.frame(dGenes[dGenes[, "padj"] < 0.05, ])
  
  sigXGenes <- as.data.frame(xGenes[xGenes[, "padj"] < 0.05, ])
  
  # Make column with gene names from row names.
  
  sigDGenes$Bter_gene_ID <- row.names(sigDGenes)
  
  sigXGenes$Bter_gene_ID <- row.names(sigXGenes)
  
  # Write results to a csv file.
  
  nameOfDFile <- paste0("amsalem_results_D_vs_", x, "_Diapause_up-regulated_LFC", y, ".csv")
  
  write.csv(sigDGenes, file = nameOfDFile, row.names = FALSE)
  
  nameOfXFile <- paste0("amsalem_results_D_vs_", x, "_Diapause_down-regulated_LFC", y, ".csv")
  
  write.csv(sigXGenes, file = nameOfXFile, row.names = FALSE)
  
}

# EXECUTED STATEMENTS ----

# Amsalem et al. (2015) All DEG Lists ----

# Make a results directory.

dir.create("../02_outputs/70_Amsalem_DEG_lists")
# Will produce a warning if directory already exists.

setwd("../02_outputs/70_Amsalem_DEG_lists")

# Set dds.

dds <- ddsAmsalem

# Index DESeq2 results into individual lists.

# LFC = 0 (all differentially expressed genes (DEGs)).

MakeGeneList("M")

MakeGeneList("F")

# LFC = 1 (highly differentially expressed genes (HDEGs)).

MakeGeneList("M", y = 1)

MakeGeneList("F", y = 1)
