#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Conduct Fisher's exact test to compare overlap between differentially
# expressed genes (DEGs) from the current study and those of Cameron et al 2013 
# and He et al 2017, with and without diapause genes from Amsalem et al 2015 
# included.
#-------------------------------------------------------------------------------
# Inputs:
# DEGs lists from the current study, Cameron et al 2013, He et al 2017 and 
# Amsalem et al 2015.

# Outputs:
# .csv file with the results of the all the comparisons.
#-------------------------------------------------------------------------------

#First you must always reset R
rm(list=ls())

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script Sx_01

# Load and Refine Orthologs List ----

# Set working directory.

setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/HISAT2_comparison/")

# Load Bombus to Apis ortholog gene list.

bterToApisDFComplete <- read.csv("Table_S7.csv")

# Refine ortholog gene list.
# This comparison requires the "BEEBASE" and "Bter_gene_symbol" columns,
# and therefore the instances where one of these columns = NA will be removed.

bterToApisDFRefined <- bterToApisDFComplete[, c("BEEBASE", "Bter_gene_symbol")] 

bterToApisDFRefined <- na.omit(bterToApisDFRefined)  # 1109 rows omitted.

# Test for redundancy/duplication in the otholog list and remove if
# found.

if (length(unique(bterToApisDFRefined$BEEBASE)) != length(bterToApisDFRefined$BEEBASE)) {
  bterApisOrthos <- bterToApisDFRefined[!duplicated(bterToApisDFRefined), ]
} else if (length(unique(bterToApisDFRefined$Bter_gene_symbol)) != length(bterToApisDFRefined$Bter_gene_symbol)) {
  bterApisOrthos <- bterToApisDFRefined[!duplicated(bterToApisDFRefined), ]
} else {
  bterApisOrthos <- bterToApisDFRefined
}

# Rename columns to match gene lists.

colnames(bterApisOrthos) <- c("Amel_gene_ID", "Bter_gene_ID")

# Load Diapause Associated Gene Lists ----

# Load and combine log-fold change >0 lists.

# Load.

diapauseLFC0Upreg <- read.csv("diapause_specific_genes_up-regulated_LFC0.csv")

diapauseLFC0Downreg <- read.csv("diapause_specific_genes_down-regulated_LFC0.csv")

# Combine.

diapauseGenesLFC0 <- rbind(diapauseLFC0Upreg, diapauseLFC0Downreg)

# Load and combine log-fold change >1 lists.

# Load.

diapauseLFC1Upreg <- read.csv("diapause_specific_genes_up-regulated_LFC1.csv")

diapauseLFC1Downreg <- read.csv("diapause_specific_genes_down-regulated_LFC1.csv")

# Combine.

diapauseGenesLFC1 <- rbind(diapauseLFC1Upreg, diapauseLFC1Downreg)

# Rename columns

colnames(diapauseGenesLFC0) <- "Bter_gene_ID"

colnames(diapauseGenesLFC1) <- "Bter_gene_ID"

# Add Apis orthologs to the Bombus diapause genes.

diapauseGenesLFC0 <- merge(diapauseGenesLFC0, bterApisOrthos, all.x = TRUE)

diapauseGenesLFC1 <- merge(diapauseGenesLFC1, bterApisOrthos, all.x = TRUE)

# Remove unnecessary objects.

rm(diapauseLFC0Downreg, diapauseLFC0Upreg, diapauseLFC1Downreg, diapauseLFC1Upreg,
   bterToApisDFComplete, bterToApisDFRefined)

# FUNCTION DEFINITIONS ----

CompareGeneListOverlap <- function (x, y, z, a, diapause = "no") {
  # Calculates the overlap in differentially expressed genes (DEGs) determined 
  # using the HISAT2/HTSeq pathway, and optionally performs the analysis with 
  # Bombus terrestris diapause genes removed.
  #
  # Args:
  #   x: a string denoting the time point and caste data from the current study
  #      to be compared in the function ("EQ", "EW", "MQ", "MW", "LQ" or "LW").
  #   y: a string denoting the Apis mellifera study to be compared in the function
  #      ("he" or "cameron").
  #   z: a string denoting the Apis mellifera study time point and caste data 
  #      to be compared in the function ("Q", "W" (for Cameron et al 2013), 
  #      "2d_Q", "2d_W", "4d_Q", or "4d_W"). 
  #   a: a number denoting the log fold change (LFC) to be used for the data to 
  #      be compared.
  #   diapause: string denoting whether diapause-specific genes from Amsalem et
  #              al 2015 have been removed from the analysis. 
  #              ("yes" or "no" (default)).
  #
  # Returns:
  #   A data.frame summarising the results of the comparison.   
  
  # Set variables based on arguments
  
  if (y == "cameron") {
    compStudy <- "Cameron et al. 2013"
    timePoint <- "60 hr"
  } else if (y == "he") {
    compStudy <- "He et al. 2017"
  }
  
  if (substr(z, 1, 2) == "2d") {
    timePoint <- "2d"
  } else if (substr(z, 1, 2) == "4d") {
    timePoint <- "4d"
  }
  
  if (grepl("Q", z, fixed = TRUE)) {
    caste <- "queen"
  } else if (grepl("W", z, fixed = TRUE)) {
    caste <- "worker"
  }
  
  if (a == 0) {
    levelDEG <- "DEGs"
    diapauseList <- diapauseGenesLFC0
  } else if (a == 1) {
    levelDEG <- "HDEGs"
    diapauseList <- diapauseGenesLFC1
  }
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("This_study_list" = paste0(x, " ", levelDEG),
                          "Comparison_study" = compStudy,
                          "Comparison_study_list" = paste0(timePoint, " ", caste, " larvae"),
                          "Diapause_genes_excluded" = diapause,
                          "DE_in_both" = NA,
                          "DE_only_in_this_study" = NA,
                          "DE_only_in_comparison_study" = NA,
                          "DE_in_neither" = NA,
                          "This_study_genes_with_no_orthologs" = NA,
                          "Comparison_study_genes_with_no_orthologs" = NA,
                          "p_value" = NA,
                          "odds_ratio" = NA,
                          "alpha_value" = NA,
                          "Percentage_of_overlap" = NA)
  
  # Load requested gene lists.
  
  # Load Bombus gene list based on the arguments.
  
  bombusGeneList <- 
    read.csv(paste0("HTSeq_results_", x, "_regulated_LFC", a, ".csv"))
  
  # Load Apis gene list based on the arguments.
  
  apisGeneList <- 
    read.csv(paste0(y, "_results_", z, "_regulated_LFC", a, ".csv"))
  
  # Remove diapause genes from the gene list if the diapause option is "yes".
  
  if (diapause == "yes") {
    bombusGeneList <- bombusGeneList[!(bombusGeneList$Bter_gene_ID %in% diapauseList$Bter_Gene_ID), ]
    apisGeneList <- apisGeneList[!(apisGeneList$Amel_gene_ID %in% diapauseList$Amel_gene_ID), ]
  }
  
  # Determine orthologs in Bombus and Apis lists.
  
  # Bombus genes with Apis orthologs.
  
  bombusGenesWithApisOrthos <- 
    bombusGeneList[bombusGeneList$Bter_gene_ID %in% bterApisOrthos$Bter_gene_ID, ]
  
  # Apis genes with Bombus orthologs.
  
  apisGenesWithBombusOrthos <-
    apisGeneList[apisGeneList$Amel_gene_ID %in% bterApisOrthos$Amel_gene_ID, ]
  
  # And number of genes without orthologs.
  
  resultsDF$This_study_genes_with_no_orthologs <- 
    length(bombusGeneList$Bter_gene_ID) - length(bombusGenesWithApisOrthos$Bter_gene_ID)
  
  resultsDF$Comparison_study_genes_with_no_orthologs <- 
    length(apisGeneList$Amel_gene_ID) - length(apisGenesWithBombusOrthos$Amel_gene_ID)
  
  # Calculate overlap between Bombus and Apis lists.
  
  # Add Bter IDs to the Apis gene list.
  
  apisGenesWithBombusOrthosAdded <- merge(apisGenesWithBombusOrthos, bterApisOrthos, by = "Amel_gene_ID")
  
  # Determine overlapping genes.
  
  overlappingGenes <- 
    bombusGenesWithApisOrthos[bombusGenesWithApisOrthos$Bter_gene_ID %in% apisGenesWithBombusOrthosAdded$Bter_gene_ID, ]
  
  # Generate results for resultsDF
  
  # Number of overlapping DEGs.
  
  resultsDF$DE_in_both <- length(overlappingGenes$Bter_gene_ID)
  
  # Number of non-overlapping Bombus DEGs.
  
  resultsDF$DE_only_in_this_study <- 
    length(bombusGenesWithApisOrthos$Bter_gene_ID) - length(overlappingGenes$Bter_gene_ID)
  
  # Number of non-overlapping Apis DEGs.
  
  resultsDF$DE_only_in_comparison_study <-
    length(apisGenesWithBombusOrthos$Amel_gene_ID) - length(overlappingGenes$Bter_gene_ID)
  
  # Calculate the size of background list.
  # (i.e. The genes with orthologs in both species that were also expressed 
  # in both studies)
  
  # Load Bombus background list.
  
  bombusBackgroundList <- read.csv("collins_htseq_background_genelist.csv")
  
  # Load Apis background list.
  
  apisBackgroundList <- read.csv(paste0(y, "_background_genelist.csv"))
  
  # Remove diapause genes from the background gene lists if the diapause option is "yes".
  
  if (diapause == "yes") {
    bombusBackgroundList <- 
      as.data.frame(bombusBackgroundList[!(bombusBackgroundList$Gene %in% diapauseList$Bter_gene_ID), ])
    apisBackgroundList <- 
      as.data.frame(apisBackgroundList[!(apisBackgroundList$Gene %in% diapauseList$Amel_gene_ID), ])
  } 
  
  # Check that the background lists only contain unique gene IDs, and remove
  # any that are redundant.
  
  if (length(unique(bombusBackgroundList$Gene)) != length(bombusBackgroundList$Gene)) {
    bombusBackgroundList <- as.data.frame(unique(bombusBackgroundList$Gene))
  }
  
  if (length(unique(apisBackgroundList$Gene)) != length(apisBackgroundList$Gene)) {
    apisBackgroundList <- as.data.frame(unique(apisBackgroundList$Gene))
  }
  
  # Rename columns of background lists.
  
  colnames(bombusBackgroundList) <- "Gene"
  
  colnames(apisBackgroundList) <- "Gene"
  
  # Determine which of the pairs of orthologs between Bombus and Apis were
  # expressed in the Bombus and Apis studies.
  
  expressedInBombus <- 
    bterApisOrthos[bterApisOrthos$Bter_gene_ID %in% bombusBackgroundList$Gene, ]
  
  notExpressedInBombus <- 
    bterApisOrthos[!(bterApisOrthos$Bter_gene_ID %in% bombusBackgroundList$Gene), ]
  
  notBombusButApis <-
    notExpressedInBombus[notExpressedInBombus$Amel_gene_ID %in% apisBackgroundList$Gene, ]
  
  # Total number of ortholog pairs expressed in the two studies is: 
  
  totalGeneListSize <- 
    length(expressedInBombus$Amel_gene_ID) + length(notBombusButApis$Amel_gene_ID)
  
  # Add results to resultsDF
  
  resultsDF$DE_in_neither <-
    totalGeneListSize - resultsDF$DE_in_both - resultsDF$DE_only_in_this_study - resultsDF$DE_only_in_comparison_study
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both lists, only 
  # Bombus list, only Apis list, and neither list.
  
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
  
  resultsDF$odds_ratio <- fisherResults$estimate[[1]]
  
  resultsDF$Percentage_of_overlap <- 
    resultsDF$DE_in_both/((resultsDF$DE_in_both + resultsDF$DE_only_in_this_study))*100
  
  # Return resultsDF
  
  return(resultsDF)
  
}

# EXECUTED STATEMENTS ----

# Comparison with Cameron et al 2013 ----

# Early-instar larvae from the current study compared to the 60h timepoint 
# from Cameron et al. 2013.

# Log fold change > 0.

# With diapause genes included. 

eqCameronQLFC0Results <- CompareGeneListOverlap("EQ", "cameron", "Q", 0)
ewCameronWLFC0Results <- CompareGeneListOverlap("EW", "cameron", "W", 0)

# With diapause genes excluded.

eqCameronQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("EQ", "cameron", "Q", 0, "yes")
ewCameronWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("EW", "cameron", "W", 0, "yes")

# Mid-instar larvae from the current study compared to the 60h timepoint 
# from Cameron et al. 2013.

# Log fold change > 0.

# With diapause genes included. 

mqCameronQLFC0Results <- CompareGeneListOverlap("MQ", "cameron", "Q", 0)
mwCameronWLFC0Results <- CompareGeneListOverlap("MW", "cameron", "W", 0)

# With diapause genes excluded.

mqCameronQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MQ", "cameron", "Q", 0, "yes")
mwCameronWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MW", "cameron", "W", 0, "yes")

# Late-instar larvae from the current study compared to the 60h timepoint 
# from Cameron et al. 2013.

# Log fold change > 0.

# With diapause genes included. 

lqCameronQLFC0Results <- CompareGeneListOverlap("LQ", "cameron", "Q", 0)
lwCameronWLFC0Results <- CompareGeneListOverlap("LW", "cameron", "W", 0)

# With diapause genes excluded.

lqCameronQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LQ", "cameron", "Q", 0, "yes")
lwCameronWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LW", "cameron", "W", 0, "yes")

# Comparison with He et al. 2017 ----

# Early-instar larvae from the current study and 2d larvae from
# He et al. 2017.

# Log fold change > 0.

# With diapause genes included.

eqHe2dQLFC0Results <- CompareGeneListOverlap("EQ", "he", "2d_Q", 0)
ewHe2dWLFC0Results <- CompareGeneListOverlap("EW", "he", "2d_W", 0)

# With diapause genes excluded.

eqHe2dQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("EQ", "he", "2d_Q", 0, "yes")
ewHe2dWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("EW", "he", "2d_W", 0, "yes")

# Mid-instar larvae from the current study and 2d larvae from
# He et al. 2017.

# Log fold change > 0.

# With diapause genes included.

mqHe2dQLFC0Results <- CompareGeneListOverlap("MQ", "he", "2d_Q", 0)
mwHe2dWLFC0Results <- CompareGeneListOverlap("MW", "he", "2d_W", 0)

# With diapause genes excluded.

mqHe2dQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MQ", "he", "2d_Q", 0, "yes")
mwHe2dWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MW", "he", "2d_W", 0, "yes")

# Mid-instar larvae from the current study and 4d larvae from
# He et al. 2017.

# Log fold change > 0.

# With diapause genes included.

mqHe4dQLFC0Results <- CompareGeneListOverlap("MQ", "he", "4d_Q", 0)
mwHe4dWLFC0Results <- CompareGeneListOverlap("MW", "he", "4d_W", 0)

# With diapause genes excluded.

mqHe4dQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MQ", "he", "4d_Q", 0, "yes")
mwHe4dWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("MW", "he", "4d_W", 0, "yes")

# Late-instar larvae from the current study and 2d larvae from
# He et al. 2017.

# Log fold change > 0.

# With diapause genes included.

lqHe2dQLFC0Results <- CompareGeneListOverlap("LQ", "he", "2d_Q", 0)
lwHe2dWLFC0Results <- CompareGeneListOverlap("LW", "he", "2d_W", 0)

# With diapause genes excluded.

lqHe2dQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LQ", "he", "2d_Q", 0, "yes")
lwHe2dWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LW", "he", "2d_W", 0, "yes")

# Mid-instar larvae from the current study and 4d larvae from
# He et al. 2017.

# Log fold change > 0.

# With diapause genes included.

lqHe4dQLFC0Results <- CompareGeneListOverlap("LQ", "he", "4d_Q", 0)
lwHe4dWLFC0Results <- CompareGeneListOverlap("LW", "he", "4d_W", 0)

# With diapause genes excluded.

lqHe4dQLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LQ", "he", "4d_Q", 0, "yes")
lwHe4dWLFC0ResultsDiapauseExcluded <- CompareGeneListOverlap("LW", "he", "4d_W", 0, "yes")

# Combine all results data.frames.

# Search global environment for all objects with "Results" in the name

resultsTables <- grep("Results", names(.GlobalEnv), value = TRUE)

# Combine all the objects into a list.

resultsTablesList <- do.call("list", mget(resultsTables))

# Perform rbind on the list of objects.
# (Rather than write out all the names).

combinedResults <- do.call("rbind", resultsTablesList)

# Write combinedResults to a .csv file.

write.csv(combinedResults, "comparative_analysis_diapause_results_.csv",
          row.names = FALSE)
