#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Compare differentially expressed genes (DEGs) lists from current study,
# Kallisto pipeline and HISAT2/HTSeq pipeline.
#-------------------------------------------------------------------------------
# Inputs:
# DEGs lists from current study (character vectors), Kallisto pipeline (.csv 
# files) and HISAT2/HTSeq pipeline (.csv files).

# Outputs:
# .csv file summarising the results of the comparisons. 
# .csv files to generate Euler diagrams.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

source("S4_32_isolate_current_study_DEGs.R")

# FUNCTION DEFINITIONS ----

compare3Methods <- function(x, y = "DEGs") {
  # Compares the overlap in DEGs lists from the current study, Kallisto pipeline
  # and HISAT2/HTSeq pipeline.
  #
  # Args:
  #   x: a string denoting the time point and caste data to be used in the 
  #      function ("EQ", "EW", "MQ", "MW", "LQ" or "LW").
  #   y: a string denoting whether the genes of interest are DEGs (offset fold 
  #      change (OFC) greater than 0), or HDEGs (OFC > 1) ("DEGs" (default) or
  #      "HDEGs")).
  #
  # Returns:
  #   A data.frame of results of the comparison.
  
  # Check that the arguments supplied are correct, and set variables.
  
  if (x == "EQ") {
    currentStudyDEGs <- eqDEGs
    currentStudyHDEGs <- eqHDEGs
  } else if (x == "EW") {
    currentStudyDEGs <- ewDEGs
    currentStudyHDEGs <- ewHDEGs
  } else if (x == "MQ") {
    currentStudyDEGs <- mqDEGs
    currentStudyHDEGs <- mqHDEGs
  } else if (x == "MW") {
    currentStudyDEGs <- mwDEGs
    currentStudyHDEGs <- mwHDEGs
  } else if (x == "LQ") {
    currentStudyDEGs <- lqDEGs
    currentStudyHDEGs <- lqHDEGs
  } else if (x == "LW") {
    currentStudyDEGs <- lwDEGs
    currentStudyHDEGs <- lwHDEGs
  } else {
    stop("Argument x must be either 'EQ', 'EW', 'MQ', 'MW', 'LQ' or 'LW'.")
  }
  
  # Choose current study gene list based on y argument.
  
  if (y == "DEGs") {
    currentStudy <- currentStudyDEGs
    levelOfDEG <- "DEGs"
  } else if (y == "HDEGs") {
    currentStudy <- currentStudyHDEGs
    levelOfDEG <- "HDEGs"
  } else {
    stop("Argument y must be either 'DEGs' or 'HDEGs'.")
  }
  
  # Load Kallisto results.
  
  setwd("../50_Kallisto_DEG_lists")
  
  kallisto <- read.csv(paste0("Kallisto_results_", x, "_regulated_LFC0.csv"))
  
  # Load HISAT2/HTSeq_count results.
  
  setwd("../60_HISAT2_HTSeq_DEG_lists")
  
  htseq <- read.csv(paste0("HTSeq_results_", x, "_regulated_LFC0.csv"))
  
  # Isolate gene lists from loaded results.
  
  kallistoDEGs <- kallisto[, "Bter_gene_ID"]
  
  htseqDEGs <- htseq[, "Bter_gene_ID"]
  
  # Compare the 3 gene lists in pairs.
  
  htseqANDKallisto <- intersect(htseqDEGs, kallistoDEGs)
  
  currentStudyANDKallisto <- intersect(currentStudy, kallistoDEGs)
  
  currentStudyANDHtseq <- intersect(currentStudy, htseqDEGs)
  
  # Compare which genes are shared by all methods.
  
  sharedByAllThree <- intersect(currentStudyANDHtseq, currentStudyANDKallisto)
  
  # Determine genes shared by 2 out of 3 methods.
  
  htseqANDKallistoNOTCurrentStudy <- setdiff(htseqANDKallisto, sharedByAllThree)
  
  currentStudyANDKallistoNOTHtseq <- setdiff(currentStudyANDKallisto, sharedByAllThree)
  
  currentStudyANDHtseqNOTKallisto	<- setdiff(currentStudyANDHtseq, sharedByAllThree)
  
  # Determine genes only present in 1 of the methods.
  
  currentStudyOnly <- length(currentStudy) - length(currentStudyANDHtseqNOTKallisto) - length(currentStudyANDKallistoNOTHtseq) - length(sharedByAllThree)
  
  kallistoOnly <- length(kallistoDEGs) - length(htseqANDKallistoNOTCurrentStudy) - length(currentStudyANDKallistoNOTHtseq) - length(sharedByAllThree)
  
  htseqOnly <- length(htseqDEGs) - length(htseqANDKallistoNOTCurrentStudy) - length(currentStudyANDHtseqNOTKallisto) - length(sharedByAllThree)
  
  # Collect data for euler diagram and output to file.
  
  eulerData <- 
    data.frame("category" = c("current_study_only", "htseq_only",
                              "kallisto_only", "current_study&kallisto",
                              "current_study&htseq",
                              "htseq&kallisto",
                              "current_study&htseq&kallisto"),
               "gene" = c(currentStudyOnly, htseqOnly, kallistoOnly,
                 length(currentStudyANDKallistoNOTHtseq),
                 length(currentStudyANDHtseqNOTKallisto),
                 length(htseqANDKallistoNOTCurrentStudy),
                 length(sharedByAllThree)))
  
  setwd("../70_euler_diagram_data")
  
  write.csv(eulerData, 
            paste0(x, "_larvae_", levelOfDEG, "_euler_data.csv"), 
            row.names = FALSE)
  
  # Percentages of other methods detected by current study method.
  # (CS = current study)
  
  htseqCSPercent <- round((length(sharedByAllThree) + length(currentStudyANDHtseqNOTKallisto))/length(htseqDEGs)*100, digits = 1)
  
  kallistoCSPercent <- round((length(sharedByAllThree) + length(currentStudyANDKallistoNOTHtseq))/length(kallistoDEGs)*100, digits = 1)
  
  csHtseqPercent <- round((length(sharedByAllThree) + length(currentStudyANDHtseqNOTKallisto))/length(currentStudy)*100, digits = 1)
  
  csKallistoPercent <- round((length(sharedByAllThree) + length(currentStudyANDKallistoNOTHtseq))/length(currentStudy)*100, digits = 1)
  
  htseqKallistoPercent <- round((length(sharedByAllThree) + length(htseqANDKallistoNOTCurrentStudy))/length(htseqDEGs)*100, digits = 1)
  
  kallistoHtseqPrecent <- round((length(sharedByAllThree) + length(htseqANDKallistoNOTCurrentStudy))/length(kallistoDEGs)*100, digits = 1)
  
  # Results data.frame
  
  resultsDF <- data.frame("time_point_caste" = rep(x, times = 3),
                          "level_of_DEGs" = rep(levelOfDEG, times = 3),
                          "method1" = c("current_study", "current_study", "HTSeq"),
                          "method2" = c("HTSeq", "Kallisto", "Kallisto"),
                          "shared_genes" = rep(NA, times = 3),
                          "total_genes_method_1" = rep(NA, times = 3), 
                          "total_genes_method_2" = rep(NA, times = 3))
  
  # Fill in the NAs.
  
  resultsDF[1, 5] <- length(currentStudyANDHtseq)
  
  resultsDF[2, 5] <- length(currentStudyANDKallisto)
  
  resultsDF[3, 5] <- length(htseqANDKallisto)
  
  resultsDF[1, 6] <- length(currentStudy)
  
  resultsDF[2, 6] <- length(currentStudy)
  
  resultsDF[3, 6] <- length(htseqDEGs)
  
  resultsDF[1, 7] <- length(htseqDEGs)
  
  resultsDF[2, 7] <- length(kallistoDEGs)
  
  resultsDF[3, 7] <- length(kallistoDEGs)
  
  # Add new columns with percentages
  
  resultsDF[, "percent_method_1_in_method_2"] <- round(resultsDF[, 5]/resultsDF[, 7]*100, digits = 1)
  
  resultsDF[, "percent_method_2_in_method_1"] <- round(resultsDF[, 5]/resultsDF[, 6]*100, digits = 1)
  
  # Return results data.frame

  return(resultsDF)
  
}

# EXECUTED STATEMENTS ----

# Make a Result Directories ----

dir.create("../../02_outputs/70_euler_diagram_data")

dir.create("../../02_outputs/80_comparison_results")
# Will produce a warning if directory already exists.

# Compare DEGs ----

setwd("../../02_outputs/50_Kallisto_DEG_lists")

ewDEGsResults <- compare3Methods("EW")

eqDEGsResults <- compare3Methods("EQ")

mwDEGsResults <- compare3Methods("MW")

mqDEGsResults <- compare3Methods("MQ")

lqDEGsResults <- compare3Methods("LQ")

lwDEGsResults <- compare3Methods("LW")

# Compare HDEGs

ewHDEGsResults <- compare3Methods("EW", "HDEGs")

eqHDEGsResults <- compare3Methods("EQ", "HDEGs")

mwHDEGsResults <- compare3Methods("MW", "HDEGs")

mqHDEGsResults <- compare3Methods("MQ", "HDEGs")

lqHDEGsResults <- compare3Methods("LQ", "HDEGs")

lwHDEGsResults <- compare3Methods("LW", "HDEGs")

# Combine all results into one table.

combinedResults <- 
  rbind(ewDEGsResults, ewHDEGsResults, eqDEGsResults, eqHDEGsResults,
        mwDEGsResults, mwHDEGsResults, mqDEGsResults, mqHDEGsResults,
        lwDEGsResults, lwHDEGsResults, lqDEGsResults, lqHDEGsResults)

# Output the results to a csv file.

setwd("../80_comparison_results")

write.csv(combinedResults, "comparison_of_methods_results.csv",
          row.names = FALSE)
