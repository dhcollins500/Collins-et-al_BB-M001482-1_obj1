#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 22.08.2018 
# File last updated: 31.03.2020
# Project: NEL0067851 (Experiment 1a) 
# Analysis: Comparative analysis
# Subsection: Comparing Bombus and Apis studies
# Tasks: Create a reciprocal best hit BLAST data.frame to convert between
# Bombus and different Apis gene IDs.
# IMPORTANT NOTE: Need to change subdirectoryPath variable to be specific to
# the computer being used for the analysis.
#-------------------------------------------------------------------------------
# Inputs:
# Output (tabular with comment lines format) from BLASTp of Amel v4.5 cf Bombus 
# terrestris proteins (and vice versa). Apis mellifera v4.5 protein table and 
# gff file from NCBI. Apis to Apis RBH analysis. Bombus terrestris annotations.

# Outputs:
# .csv file with the results of the analysis, allowing conversion between 
# Bombus and different Apis mellifera gene IDs.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ape)  # read.gff()
library(tidyr) # separate()
sessionInfo()

# LOADING CUSTOM FUNCTIONS ----

source("f_CalculateRBH.R")

# LOADING DATA ----

# Loading Apis mellifera Files ----

# GFF file.

setwd("../02_exp1a_data")

apisGFF <- dir(pattern = "62")

apisGFFAnnotation <- read.gff(apisGFF)

apisGeneAnnotations <- subset(apisGFFAnnotation, type == "gene")  # genes

# Protein table.

apisProteinTable <- read.delim("63_Apis_ProteinTable48_22683.txt", 
                               stringsAsFactors = FALSE)

# Apis to Apis RBH file.

setwd("../04_analysis_output")

apisToApisRBH <- read.csv("61_Apis_cf_Apis_RBH.csv", stringsAsFactors = FALSE)

# Bombus terrestris Features Table ----

setwd("../02_exp1a_data")

bombusAnnotations <- 
  read.table("40_GCF_000214255.1_Bter_1.0_feature_table.tabular", sep = "\t")

colnames(bombusAnnotations) <- 
  c("feature",	"class", "assembly",	"assembly_unit",	"seq_type", "chromosome",
    "genomic_accession",	"start",	"end",	"strand", "product_accession",
    "non-redundant_refseq",	"related_accession",	"gene_description",	
    "gene_symbol", "GeneID", "locus_tag",	"feature_interval_length",
    "product_length",	"attributes")

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Calculating RBH BLAST Orthologs Between B. terrestris and A. mellifera 
# ----

# Running RBH Script ----

# Setting input values for the function.

x <- "70_BlastOutput_Amelv4_5_protein_on_bt_protein"  

y <- "71_BlastOutput_Bt_protein_on_amelv4_5_protein"

z <- "90_Amel_v4_5_to_Bter_RBH_Results"

a <- "U:/Documents/NEL0067581/10_Exp1a"
### IMPORTANT - change this value to the appropriate address for your computer #

# Call CalculateRBH function.

CalculateRBH(x, y, z, a)

# Parsing the RBH Results ----

# Loading the results of the RBH script

setwd(paste(subdirectoryPath, pathToOutput, sep = ""))

results <- dir(pattern = "90_Amel_v4_5")

RBHResults <- read.table(results, stringsAsFactors = FALSE)

colnames(RBHResults) <- c("Bter_protein_ID", "Amel_v4_5")

# Parsing the Apis GFF Gene Annotations ----

# Separating out the desired information from the attributes column.

separatedApisGeneAnnotations <- 
  separate(apisGeneAnnotations, attributes, 
           c("Apis_gene_number", "Apis_gene_ID", "Apis_gene_symbol"),
           sep = ";", remove = TRUE)

# Removing columns not needed.

separatedApisGeneAnnotations <- separatedApisGeneAnnotations[, c(10:11)]

# separating out the two columns to isolate the BEEBASE and gene symbol info.

sep1 <- 
  separate(separatedApisGeneAnnotations, Apis_gene_ID, c("to_delete", "BEEBASE"),
           sep = "BEEBASE:")

sep2 <- 
  separate(sep1, BEEBASE, c("BEEBASE", "to_delete2"), sep = ",")

sep3 <- 
  separate(sep2, Apis_gene_symbol, c("to_delete3", "Apis_gene_symbol"), sep = "=")

# remove unwanted columns.

sep3 <- sep3[, c(2, 5)]

# remove rows with NAs.

beebaseGeneSymbol <- na.omit(sep3)

colnames(beebaseGeneSymbol) <- c("BEEBASE", "Apis_locus")

# Parsing Protein Table ----

# Remove unwanted columns.

apisProteinTable <- apisProteinTable[, c(6:8, 10)]

# Adding extra column to turn gene IDs into Gene symbols.

apisProteinTable$gene_symbol <- 
  rep("LOC", times = length(apisProteinTable$GeneID))

# Pasting gene ID and gene symbol together.

apisProteinTable$gene_symbol <- 
  paste(apisProteinTable$gene_symbol, apisProteinTable$GeneID, sep = "")

# Renaming the columns.

colnames(apisProteinTable) <- c("Apis_gene_ID", "Apis_locus", "Amel_v4_5", 
                                "Apis_gene_description", "Apis_gene_symbol")

# Parsing the Apis to Apis RBH data.frame.

apisToApisRBH <- apisToApisRBH[, c(3, 6)]

colnames(apisToApisRBH) <- c("OGS_v1_1", "Apis_gene_symbol")

# Parsing Bombus Annotations ----

# Remove unnecessary columns from bombusAnnotations.

bombusAnnotations <- 
  bombusAnnotations[, c("feature", "product_accession", "related_accession",
                        "gene_description", "gene_symbol")]

# Subset bombusAnnotations by CDS (to isolate all the proteins).

bombusProteinAnnotations <- subset(bombusAnnotations, feature == "CDS")

# Remove feature column and rename the other columns.

bombusProteinAnnotations <- bombusProteinAnnotations[, 2:5]

colnames(bombusProteinAnnotations) <- 
  c("Bter_protein_ID", "Bter_gene_ID", "Bter_gene_description", 
    "Bter_gene_symbol")

# Merging the RBHResults with the Annotations ----

# Merging RBH and protein table.

RBHPlusProteinTable <- merge(RBHResults, apisProteinTable, 
                             by = "Amel_v4_5", all = TRUE)

# Merging RBH/protein table with BEEBASE info.

RBHPlusProteinPlusGFF <- merge(RBHPlusProteinTable, beebaseGeneSymbol, 
                               by = "Apis_locus", all = TRUE)

# Merging RBHPlusProteinPlusGFF with apisToApisRBH.

apisRBH <- merge(RBHPlusProteinPlusGFF, apisToApisRBH,
                 by = "Apis_gene_symbol", all = TRUE)

# Merging the Bter annotations with the RBH results/Apis annotations.

bombusApisRBHInitial <- merge(apisRBH, bombusProteinAnnotations,
                              by = "Bter_protein_ID", all = TRUE)

# Parse the BombusApisRBH to remove rows that are not needed (e.g. where either
# Bombus or all Apis IDs = NA).

bombusApisRBHremoval1 <- 
  bombusApisRBHInitial[!(is.na(bombusApisRBHInitial$Bter_gene_ID)), ]

bombusApisRBHremoval2 <- 
  bombusApisRBHremoval1[!(is.na(bombusApisRBHremoval1$Apis_gene_symbol) & is.na(bombusApisRBHremoval1$BEEBASE) & is.na(bombusApisRBHremoval1$OGS_v1_1)), ]

# Output the RBH data.frame to a file.

write.csv(bombusApisRBHremoval2, "91_Bombus_cf_Apis_RBH.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 25.07.18
# File last updated: 21.01.20
# Project: NEL0067851 (Experiment 1a) 
# Analysis: Enriched Gene Ontology (GO) term analysis.
# Used in script number: 30
#-------------------------------------------------------------------------------

CalculateRBH <- function (x, y, z, a) {
  # Creates a list of likely orthologs between protein sequences in two species
  # based on reciprocal best hit (RBH) BLAST results
  #
  # Args:
  #   x: A string denoting an output file (tabular format with comment lines) 
  #      from BLASTp, where the protein sequences of species 1 have been blasted 
  #      against a database made up of the protein sequences of species 2.
  #   y: A string denoting an output file (tabular format with comment lines) 
  #      from BLASTp, where the protein sequences of species 2 have been 
  #      blasted against a database made up of the protein sequences of species 
  #      1.
  #   z: A string denoting the name of the output file
  #   a: A string denoting the path to the directory where the analysis is 
  #      occuring.
  #
  # Returns:
  #   A .txt file (tab delimited) containing the pairs of likely RBH orthologs
  #   between species 1 and 2.
  
  # Set working directory to the location of my python executable.
  
  setwd("C:/Anaconda3/")
  
  # Designate the command for the system call (in the case "python").
  
  command = "python"
  
  # Designate the path to the RBH BLAST python script (taken from  
  # https://github.com/hongqin/Simple-reciprocal-best-blast-hit-pairs).
  
  pathToDirectories <- a 
  
  pathToScript <- paste(a, "/03_exp1a_scripts/31_NEL0067581_exp1a_RBHv1.py", 
                        sep = "")
  
  # Designate the paths to the input files.
  
  pathToInput <- "/02_exp1a_data/"
  
  # Create complete paths for input files.
  
  inputFile1 <- paste(a, pathToInput, x, sep = "")
  
  inputFile2 <- paste(a, pathToInput, y, sep = "")
  
  # Designate the path for the output file.
  
  pathToOutput <- "/04_analysis_output/"
  
  # Create complete path for output file.
  
  outputFile <- paste(a, pathToOutput, z, "_", Sys.Date(), ".txt", sep = "")
  
  # Combine all information about system arguments.
  
  allArgs <- c(pathToScript, inputFile1, inputFile2, outputFile)
  
  # System call to run the python script.
  
  system2(command, allArgs, stdout = TRUE)
}