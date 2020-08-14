#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Produce lists of DEGs from Kallisto and HISAT2/HTSeq data.
#-------------------------------------------------------------------------------
# Inputs:
# ddsKallisto and ddsHTSeq objects.

# Outputs:
# .csv files containing lists of DEGs for a given method/caste/time point.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# DESeq2 and tximport loaded in the scripts in the "Data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Kallisto Data ----

source("S4_11_kallisto_tximport_DESeq2.R")

# HISAT2/HTSeq Data ----

setwd("../../01_scripts")

source("S4_22_htseq_DESeq2.R")
sessionInfo()

# FUNCTION DEFINITIONS ----

MakeGeneList <- function (x, y = 0, z, method) {
  # Produces and saves a list of DEGs that meet a specific log-fold change (LFC)
  # for a given time point, caste and DEGs analysis method.
  #
  # Args:
  #   x = string denoting the time point ("E", "M" or "L").
  #   y = numeric denoting the LFC to use when filter DESeq2 results
  #       (default of 0 (i.e. all significant genes)).
  #   z = string denoting whether the gene list is DEGs expressed more in
  #       worker or queen larvae ("W" or "Q").
  #   method = string denoting which analysis method was used to produe the 
  #            DEG list ("Kallisto" or "HTSeq")
  
  # Returns:
  #   A .csv file of Bombus terrestris gene names for the significant 
  #   genes that meet the LFC threshold set.
  
  # Note: false discovery rate (FDR) of 0.05 used, but this can be changed by
  # changing the alpha level in the results() call below.
  
  # Set variables.
  
  if (x == "E") {
    a <- "EQ"
    b <- "EW"
  } else if (x == "M") {
    a <- "MQ"
    b <- "MW"
  } else if (x == "L") {
    a <- "LQ"
    b <- "LW"
  } else {
    stop('Argument x must equal "E", "M" or "L"')
  }
  
  # Extract results.
  
  if (method == "Kallisto") {
    dds <- ddsKallisto
  } else if (method == "HTSeq") {
    dds <- ddsHTSeq
  } else {
    stop('Argument x must be either "Kallisto" or "HTSeq".')
  }
  
  allExtractedResults <- results(dds, contrast = c("condition", a, b),
                                 alpha = 0.05, lfcThreshold = y)
  
  # Remove NAs from padj column (represent genes excluded from analysis as all 
  # counts were 0 or it contained an extreme count outlier).
  
  allExtractedResults <- allExtractedResults[!is.na(allExtractedResults$padj), ]
  
  # Index results based on argument z.
  
  if (z == "Q") {
    indexedResults <- 
      allExtractedResults[allExtractedResults[, "log2FoldChange"] > 0, ]
  } else if (z == "W") {
    indexedResults <- 
      allExtractedResults[allExtractedResults[, "log2FoldChange"] < 0, ]
  } else {
    stop('Argument z must equal "Q" or "W"')
  }
  
  # Index only significant (padj < 0.05) genes.
  
  sigResults <- as.data.frame(indexedResults[indexedResults[, "padj"] < 0.05, ])
  
  # Make column with gene names from row names.
  
  sigResults$Bter_gene_ID <- row.names(sigResults)
  
  # Write sigResultsNames to a csv file.
  
  nameOfFile <- paste0(method, "_results_", x, z, "_regulated_LFC", y, ".csv")
  
  write.csv(sigResults, file = nameOfFile, row.names = FALSE)
  
}

# EXECUTED STATEMENTS ----

# Kallisto DEG Lists ----

# Make a results directory.

dir.create("../02_outputs/50_Kallisto_DEG_lists")
# Will produce a warning if directory already exists.

setwd("../02_outputs/50_Kallisto_DEG_lists")

# Index DESeq2 results into individual lists.

# LFC = 0 (all differentially expressed genes (DEGs)).

MakeGeneList("E", 0, "Q", method = "Kallisto")

MakeGeneList("E", 0, "W", method = "Kallisto")

MakeGeneList("M", 0, "Q", method = "Kallisto")

MakeGeneList("M", 0, "W", method = "Kallisto")

MakeGeneList("L", 0, "Q", method = "Kallisto")

MakeGeneList("L", 0, "W", method = "Kallisto")

# HISAT2/HTSeq DEG Lists ----

# Make a results directory.

dir.create("../60_HISAT2_HTSeq_DEG_lists")
# Will produce a warning if directory already exists.

setwd("../60_HISAT2_HTSeq_DEG_lists")

# Index DESeq2 results into individual lists.

# LFC = 0 (all differentially expressed genes (DEGs)).

MakeGeneList("E", 0, "Q", method = "HTSeq")

MakeGeneList("E", 0, "W", method = "HTSeq")

MakeGeneList("M", 0, "Q", method = "HTSeq")

MakeGeneList("M", 0, "W", method = "HTSeq")

MakeGeneList("L", 0, "Q", method = "HTSeq")

MakeGeneList("L", 0, "W", method = "HTSeq")
