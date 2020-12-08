#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 3, 4 and 5
# Tasks: Produce diapause list of differentially expressed genes (DEGs) from 
# Amsalem et al. (2015).
#-------------------------------------------------------------------------------
# Inputs:
# Lists of DEGs from Amsalem et al. (2015).

# Outputs:
# .csv files containing lists of DEGs specific to diapause from 
# Amsalem et al. (2015).
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# No data loaded, as loaded by function.

# FUNCTION DEFINITION ----

GenerateDiapauseGeneList <- function (x, y) {
  # Generates a list of diapause specific differentially expressed genes (DEGs)
  # by finding the intersection of DEGs list for diapause vs. mating and
  # diapause vs. founding (post mating) for a specific direction of differential
  # expression and log fold change (LFC).
  #
  # Args:
  #   x: string denoting which direction of differential expression is to
  #      be used for the gene lists (more expressed in diapause = up-regulated 
  #     ("up"), less expressed in diapause = down-regulated ("down")).
  #   y: numeric denoting the LFC to use when selecting the list of DEGs 
  #      (0 or 1).
  #
  # Returns:
  #   a .csv file containing the list of DEGs.
  
  # Load requested gene lists.
  
  vsMatingList <- 
    read.csv(paste0("amsalem_results_D_vs_M_Diapause_",
                    x, "-regulated_LFC", y, ".csv"))
  
  vsFoundingList <- 
    read.csv(paste0("amsalem_results_D_vs_F_Diapause_",
                    x, "-regulated_LFC", y, ".csv"))
  
  # Determine genes that overlap.
  
  overlappingGenes <- 
    vsMatingList[vsMatingList$Bter_gene_ID %in% vsFoundingList$Bter_gene_ID, ]
  
  # Remove unnecessary columns.
  
  overlappingGenes <- as.data.frame(overlappingGenes$Bter_gene_ID)
  
  colnames(overlappingGenes) <- "Bter_Gene_ID"
  
  # Write results to file.
  
  setwd("../80_diapause_DEG_lists")
  
  nameOfFile <- paste0("diapause_specific_genes_", x, "-regulated_LFC", y, ".csv")
  
  write.csv(overlappingGenes, nameOfFile, row.names = FALSE)
  
  setwd("../70_Amsalem_DEG_lists")
  
}

# EXECUTED STATEMENTS ----

# Amsalem et al. (2015) Diapause DEG List ----

# Note, following the approach of Amsalem et al. (2015), a gene is considered
# specific to diapause if it is differentially regulated in the same direction
# (up-regulated or down-regulated) with respect to both mating and founding 
# (post diapause).

# Make a results directory.

dir.create("../02_outputs/80_diapause_DEG_lists")
# Will produce a warning if directory already exists.

# Location of gene lists.

setwd("../02_outputs/70_Amsalem_DEG_lists")

# Generate lists.

GenerateDiapauseGeneList("up", 0)

GenerateDiapauseGeneList("up", 1)

GenerateDiapauseGeneList("down", 0)

GenerateDiapauseGeneList("down", 1)
