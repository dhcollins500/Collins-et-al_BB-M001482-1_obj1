#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1, 2 and 4
# Tasks: Produce lists of differentially expressed genes (DEGs) from Cameron 
# et al. (2013) and He et al. (2017) data.
#-------------------------------------------------------------------------------
# Inputs:
# ddsCameron and ddsHe objects.

# Outputs:
# .csv files containing lists of DEGs for a given study/caste/time point.
# .csv files containing all genes expressed in each study (background lists).
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# DESeq2 loaded in the script in the "Data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Cameron et al. (2013) and He et al. (2017) Data ----

source("S5_50_htseq_to_DESeq2.R")
sessionInfo()

# FUNCTION DEFINITIONS ----

MakeGeneList <- function (x, y = 0, z, day) {
  # Produces and saves a list of DEGs that meet a specific log-fold change (LFC)
  # for a given time point, caste and study.
  #
  # Args:
  #   x: string denoting whether the data loaded should be from Cameron
  #      et al. (2013) or He et al. (2017) ("cameron" or "he").
  #   y: numeric denoting the LFC to use when filter DESeq2 results
  #      (default of 0 (i.e. all significant genes)).
  #   z: string denoting whether the gene list is DEGs expressed more in
  #      worker or queen larvae ("W" or "Q").
  #   day: numeric denoting which day in the He et al. (2017) data should be
  #        used (2 or 4).
  #
  # Returns:
  #   A .csv file of Apis mellifera gene names for the significant 
  #   genes that meet the LFC threshold set.
  
  # Note: false discovery rate (FDR) of 0.05 used, but this can be changed by
  # changing the alpha level in the results() call below.
  
  # Set variable based on arguments.
  
  if (x == "cameron") {
    dds <- ddsCameron
    a <- "Q"  # condition (Queen)
    b <- "W"  # condition (Worker)
  } else if (x == "he" && day == 2) {
    dds <- ddsHe
    a <- "2d_Q" # condition (day 2 Queen)
    b <- "2d_W" # condition (day 2 Worker)
  } else if (x == "he" && day == 4) {
    dds <- ddsHe
    a <- "4d_Q" # condition (day 4 Queen)
    b <- "4d_W" # condition (day 4 Worker)
  } else {
    stop('Argument x must be either "cameron" or "he".')
  }
  
  # Extract results.
  
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
  
  sigResults$Amel_gene_ID <- row.names(sigResults)
  
  # Write sigResultsNames to a csv file.
  
  if (z == "Q") {
    timePhenotype <- a
  } else if (z == "W") {
    timePhenotype <- b
  }
  
  nameOfFile <- paste0(x, "_results_", timePhenotype, "_regulated_LFC", y, ".csv")
  
  write.csv(sigResults, file = nameOfFile, row.names = FALSE)
  
}

# EXECUTED STATEMENTS ----

# Cameron et al. (2013) DEG Lists ----

# Make a results directory.

dir.create("../02_outputs/50_Cameron_DEG_lists")
# Will produce a warning if directory already exists.

setwd("../02_outputs/50_Cameron_DEG_lists")

# Generate a list of all expressed genes for statistical test background list.

# Extract results from ddsCameron.

cameronResults <- results(ddsCameron)

# Generate data.frame of the row names (i.e. the expressed genes).

allCameronGenes <- as.data.frame(row.names(cameronResults))

# Correct column name.

colnames(allCameronGenes) <- "Gene"

# Save data.frame as .csv.

write.csv(allCameronGenes, "cameron_background_genelist.csv",
          row.names = FALSE)

# Index DESeq2 results into individual lists.

# LFC = 0 (all differentially expressed genes (DEGs)).

MakeGeneList("cameron", z = "Q")

MakeGeneList("cameron", z = "W")

# LFC = 1 (highly differentially expressed genes (HDEGs)).

MakeGeneList("cameron", y = 1, z = "Q")

MakeGeneList("cameron", y = 1, z = "W")

# He et al. (2017) DEG Lists ----

# Make a results directory.

dir.create("../60_He_DEG_lists")
# Will produce a warning if directory already exists.

setwd("../60_He_DEG_lists")

# Generate a list of all expressed genes for statistical test background list.

# Extract results from ddsCameron.

heResults <- results(ddsHe)

# Generate data.frame of the row names (i.e. the expressed genes).

allHeGenes <- as.data.frame(row.names(heResults))

# Correct column name.

colnames(allHeGenes) <- "Gene"

# Save data.frame as .csv.

write.csv(allHeGenes, "he_background_genelist.csv",
          row.names = FALSE)

# Index DESeq2 results into individual lists.

# LFC = 0 (all DEGs).

MakeGeneList("he", z = "Q", day = 2)

MakeGeneList("he", z = "W", day = 2)

MakeGeneList("he", z = "Q", day = 4)

MakeGeneList("he", z = "W", day = 4)

# LFC = 1 (HDEGs).

MakeGeneList("he", y = 1, z = "Q", day = 2)

MakeGeneList("he", y = 1, z = "W", day = 2)

MakeGeneList("he", y = 1, z = "Q", day = 4)

MakeGeneList("he", y = 1, z = "W", day = 4)
