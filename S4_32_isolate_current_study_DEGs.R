#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Produce data.frames of differentially expressed genes (DEGs) from 
# current study Table S10.
#-------------------------------------------------------------------------------
# Inputs:
# Table S10 tabs in .csv format.

# Outputs:
# Character vectors of DEGs. 
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# No data loaded, as the function does it.

# FUNCTION DEFINITIONS ----

isolateDEGs <- function(x, y = "DEGs") {
  # Load a gene expression data.frame for a given time point and caste and 
  # isolate the annotated DEGs or highly differentially expressed genes (HDEGs).
  #
  # Args:
  #   x: a string denoting the time point and caste data to be used in the 
  #      function ("EQ", "EW", "MQ", "MW", "LQ" or "LW").
  #   y: a string denoting whether the genes of interest are DEGs (offset fold 
  #      change (OFC) greater than 0), or HDEGs (OFC > 1) ("DEGs" (default) or
  #      "HDEGs")).
  #
  # Returns:
  #   A character vector of the Bombus terrestris gene symbol for DEGs/HDEGs 
  #   of the specified time point and caste.
  
	# Check that x is one of "EQ", "EW", "MQ", "MW", "LQ" or "LW".
	
	if (x != "EQ" && x != "EW" && x != "MQ" && x != "MW" && x != "LQ" && x != "LW") {
		stop("Argument x must be either 'EQ', 'EW', 'MQ', 'MW', 'LQ' or 'LW'.")
	}
	
	# Load appropriate .csv file
	
	fileToLoad <- paste0("Table_S10_", x, ".csv")
	
	geneExpressionDF <- read.csv(fileToLoad)

	# Isolate columns of interest.

	geneExpressionDF <- geneExpressionDF[, c(21, 22, 23, 26, 36, 37)]

	# Filter to retain annotated genes only.
	
	annotatedGenesDF <- geneExpressionDF[geneExpressionDF[, "type"] == "gene", ]
	
	# Filter to retain differentially expressed genes (DEGs) only or highly DEGs 
	# (HDEGs) based on argument y.

	if (y == "DEGs") {
		degsDF <- annotatedGenesDF[annotatedGenesDF[, 1] == "YES", ]
	} else if (y == "HDEGs") {
		degsDF <- annotatedGenesDF[annotatedGenesDF[, 3] == "TRUE", ]
	} else {
		stop("Argument y must be either 'DEGs' or 'HDEGs'")
	}	
	
	# Isolate the gene symbol only.
	
	geneSymbolCV <- degsDF[, "Gene_symbol"]
	
	# Return character vector.

	return(geneSymbolCV)

}

# EXECUTED STATEMENTS ----

# Set Working Directory ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

setwd("../00_data/20_table_s10_csvs")

# Generate Current Study Gene Lists ----

# DEGs lists.

ewDEGs <- isolateDEGs("EW")

eqDEGs <- isolateDEGs("EQ")

mwDEGs <- isolateDEGs("MW")

mqDEGs <- isolateDEGs("MQ")

lwDEGs <- isolateDEGs("LW")

lqDEGs <- isolateDEGs("LQ")

# HDEGs lists.

ewHDEGs <- isolateDEGs("EW", "HDEGs")

eqHDEGs <- isolateDEGs("EQ", "HDEGs")

mwHDEGs <- isolateDEGs("MW", "HDEGs")

mqHDEGs <- isolateDEGs("MQ", "HDEGs")

lwHDEGs <- isolateDEGs("LW", "HDEGs")

lqHDEGs <- isolateDEGs("LQ", "HDEGs")
