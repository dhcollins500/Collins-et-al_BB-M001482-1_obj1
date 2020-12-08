#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 5
# Tasks: Conduct Fisher's Exact tests to compare overlap between differentially
# expressed genes (DEGs) from the current study and diapause genes from
# Amsalem et al. (2015).
#-------------------------------------------------------------------------------
# Inputs:
# DEGs lists from the current study and Amsalem et al. (2015).

# Outputs:
# .csv file with the results of the all the comparisons.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# Data loaded by function.

# FUNCTION DEFINITIONS ----

CompareDiapauseGeneListOverlap <- function (x, y = 0) {
  # Calculates the overlap in differentially expressed genes (DEGs) from the 
  # current study determined using the HISAT2/HTSeq pathway and diapause genes.
  #
  # Args:
  #   x: a string denoting the time point and caste data from the current study
  #      to be compared in the function ("EQ", "EW", "MQ", "MW", "LQ" or "LW").
  #   y: a number denoting the log fold change (LFC) to be used for the data to 
  #      be compared.
  #
  # Returns:
  #   A data.frame summarising the results of the comparison.   
  
  # Set variables based on arguments
  
  if (y == 0) {
    levelDEG <- "DEGs"
  } else if (y == 1) {
    levelDEG <- "HDEGs"
  }
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("This_study_list" = paste0(x, " ", levelDEG),
                          "Comparison_study" = "Amsalem et al. (2015)",
                          "Comparison_study_list" = "diapause-associated genes",
                          "DE_in_both" = NA,
                          "DE_only_in_this_study" = NA,
                          "DE_only_in_comparison_study" = NA,
                          "DE_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = 0.0083,
                          "Percentage_of_overlap" = NA)
  
  # Load and combine requested gene lists.
  
  # Load current study gene list based on argument x.
  
  currentStudyGenes <- 
    read.csv(paste0("HTSeq_results_", x, "_regulated_LFC", y, ".csv"))
  
  # Load Amsalem et al. 2015 diapause gene lists based on the arguments.
  
  diapauseUpreg <- read.csv(paste0("diapause_specific_genes_up-regulated_LFC",
                                   y, ".csv"))
  
  diapauseDownreg <- read.csv(paste0("diapause_specific_genes_down-regulated_LFC",
                                   y, ".csv"))
  
  # Combine diapause gene lists.
  
  diapauseGenes <- rbind(diapauseUpreg, diapauseDownreg)
  
  # Calculate overlap between current study and diapause lists.
  
  # Determine overlapping genes.
  
  overlappingGenes <- 
    currentStudyGenes[currentStudyGenes$Bter_gene_ID %in% diapauseGenes$Bter_Gene_ID, ]
  
  # Generate results for resultsDF
  
  # Number of overlapping DEGs.
  
  resultsDF$DE_in_both <- length(overlappingGenes$Bter_gene_ID)
  
  # Number of non-overlapping Bombus DEGs.
  
  resultsDF$DE_only_in_this_study <- 
    length(currentStudyGenes$Bter_gene_ID) - length(overlappingGenes$Bter_gene_ID)
  
  # Number of non-overlapping diapause DEGs.
  
  resultsDF$DE_only_in_comparison_study <-
    length(diapauseGenes$Bter_Gene_ID) - length(overlappingGenes$Bter_gene_ID)
  
  
  # Size of background list = all the genes in the Bter GTF used for 
  # HTSeq (12008).
  
  totalGeneListSize = 12008
  
  # Add results to resultsDF
  
  resultsDF$DE_in_neither <-
    totalGeneListSize - resultsDF$DE_in_both - resultsDF$DE_only_in_this_study - resultsDF$DE_only_in_comparison_study
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both lists, only 
  # current study list, only diapause list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$DE_in_both, resultsDF$DE_only_in_this_study,
                              resultsDF$DE_only_in_comparison_study, 
                              resultsDF$DE_in_neither))
  
  # Change the dimensions to 2 rows and 2 columns
  
  dim(contingencyTable) <- c(2,2)
  
  # Conduct two-tailed Fisher's exact test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance
  
  fisherResults <- fisher.test(contingencyTable)
  
  # Add results of statistical tests to resultsDF.
  
  resultsDF$p_value <- fisherResults$p.value
  
  resultsDF$Odds_ratio <- fisherResults$estimate[[1]]
  
  resultsDF$Percentage_of_overlap <- 
    resultsDF$DE_in_both/((resultsDF$DE_in_both + resultsDF$DE_only_in_this_study))*100
  
  # Return resultsDF
  
  return(resultsDF)
  
}

# EXECUTED STATEMENTS ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Set Working Directory ----

setwd("../02_outputs/90_diapause_comparison_DEG_lists")

# Comparison of Early-Instar Differentially Expressed Genes with Diapause 
# Genes ----

ewResults <- CompareDiapauseGeneListOverlap("EW")

eqResults <- CompareDiapauseGeneListOverlap("EQ")

# Comparison of Mid-Instar Differentially Expressed Genes with Diapause 
# Genes ----

mwResults <- CompareDiapauseGeneListOverlap("MW")

mqResults <- CompareDiapauseGeneListOverlap("MQ")

# Comparison of Late-Instar Differentially Expressed Genes with Diapause 
# Genes ----

lwResults <- CompareDiapauseGeneListOverlap("LW")

lqResults <- CompareDiapauseGeneListOverlap("LQ")

# Combine all Results Data.Frames ----

# Search global environment for all objects with "Results" in the name

resultsTables <- grep("Results", names(.GlobalEnv), value = TRUE)

# Combine all the objects into a list.

resultsTablesList <- do.call("list", mget(resultsTables))

# Perform rbind on the list of objects.
# (Rather than write out all the names).

combinedResults <- do.call("rbind", resultsTablesList)

# Write combinedResults to a .csv file.

write.csv(combinedResults, "supplementary_Table_S17.csv",
          row.names = FALSE)
