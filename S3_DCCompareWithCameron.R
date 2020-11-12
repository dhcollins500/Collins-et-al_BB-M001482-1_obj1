#-----------------------------------------------------------------------------------
# David Collins
# Started 19/01/2018
# Last modified 12/11/2020
# This script is a modification of David Prince's script to compare the mRNAseq results from Bombus terrestris with a microarray
# on Apis mellifera published by Cameron et al. 2013
#-----------------------------------------------------------------------------------

# DESCRIPTION

# Purpose: to calculate whether the overlap between two gene lists is statistically significant

## Situation - the two lists to be compared are from different species
## AND
## one of the lists is a microarray study
### Necessary data conversion - convert the gene list from one of the species into a list of orthologs in the other species
### Necessary inputs: 
#### input1 = a list of genes from this study as .txt file - For example "DC_EQ_AllSig.txt" (10 files)
#### input2 = a list of genes from the study to be compared to as a .txt file - "Cameron_108_worker" and "Cameron_004_queen" (14 files)
#### input3 = the number of genes with orthologs in both species - "new_conversion_sheet_both_part2.csv"
#### input4 = a list of the total number of genes that appear in the study to be compared to as a .txt file -
####                                                        "Cameron_et_al_2013_microarray.csv"

## Output will be:
## the number of genes that overlap between the two lists
## a p-value for how significant the overlap between the two lists is
## an odds ratio (representing the strength and direction of the association)

# CODE ##########################################################################

# clear R's memory of previous inputs
rm(list=ls())

# get working directory (tells you were R is currently looking to get data from)
getwd()

# set working directory - tell R which folder the data to input is stored in
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/comparative_analyses/BombusApisCaste/") 
# change this address as necessary
# can copy and paste addresses from Windows Explorer but need to turn the \ into /

# double check that the working directory is set to the correct place
getwd()

# The dataset is from a microarray, therefore not all genes may be present on the array
# Load the list of genes present on the microarray
## File taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52288
microarray = read.csv("Cameron_et_al_2013_microarray.csv", header = TRUE) 

microarray <- read_csv("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/comparative_analyses/BombusApisCaste/Cameron_et_al_2013_microarray.csv")

# Isolate the Reporter.Comment column to find the names of the genes present on the microarray
gene.names = as.character(microarray$Beebase)

# Some genes on the microarray are represented by more than one spot
# Therefore identify the unique gene Ids (removing redundancy)
unique.gene.names = unique(gene.names)

# Calculate which of the unique.gene.names have an RBH BLAST ortholog in Bombus terrestris
# Input Bombus terrestris to Apis mellifera RBH BLAST conversion list
conversion.df = read.csv("new_conversion_sheet_both_part2.csv", header = TRUE) 

# The Cameron et al 2013 data set uses OGS_v1_1 Ids, therefore remove the other Apis Ids
conversion.df2 = conversion.df[,c(4,6)]
# and the instances where no RBH BLAST ortholog exists between OGS_v1_1 and Bter Gene_Ids
conversion.df3 = na.omit(conversion.df2)
conversion.df3 = unique(conversion.df3) #removes redundancy
# Calculate the overlap between the RBH BLAST ortholog list (the column OGS_v1_1 in the 'new conversion sheet 2') 
#and the microarray gene list (the values labelled 'unique.gene.names')
microarray.RBH = intersect(conversion.df3$OGS_v1_1, unique.gene.names)
# count the number of ortholog pairs between Bombus terrestris and Apis mellifera on the microarray
total.gene.list.size = length(microarray.RBH)
#Note the total number of orthologs is 7598

## Define a function for comparing the gene lists

compare.bt.amel.lists = function (x, y)
{
  # Import the necessary files
  bombus.list = readLines(x)
  apis.list = readLines(y)
  
  # Calculate the Bombus terrestris genes with Apis mellifera orthologs
  bombus.present.in.df = conversion.df3[match(bombus.list, conversion.df3$Bombus_Gene_ID),]
  bombus.present.in.df2 = na.omit(bombus.present.in.df)
  count.removed.NAs.bombus = length(bombus.present.in.df$Bombus_Gene_ID) - length(bombus.present.in.df2$Bombus_Gene_ID)
  print("Number of Bombus terrestris genes in the list with no Apis mellifera homolog:")
  print(count.removed.NAs.bombus)
  
  # Place the Bombus terrestris genes with orthologs in Apis mellifera into a list
  bombus.ortho.list = as.character(bombus.present.in.df2$Bombus_Gene_ID)
  
  # Calculate the Apis mellifera genes with Bombus terrestis orthologs
  apis.present.in.df = conversion.df3[match(apis.list, conversion.df3$OGS_v1_1),]
  apis.present.in.df2 = na.omit(apis.present.in.df)
  count.removed.NAs.apis = length(apis.present.in.df$OGS_v1_1) - length(apis.present.in.df2$OGS_v1_1)
  print("Number of Apis mellifera genes in the list with no Bombus terrestris homolog:")
  print(count.removed.NAs.apis)
  
  # Place the Apis mellifera genes with orthologs in Bombus terrestris into a list
  apis.ortho.list = as.character(apis.present.in.df2$Bombus_Gene_ID)
  
  # Determine the genes present in both list1 and list2, and store the number of genes separately
  both.lists = intersect(bombus.ortho.list, apis.ortho.list)
  print("Overlapping genes are:")
  print(both.lists)
  number.both.lists = length(both.lists)
  
  # Determine the genes present in list1 but not list2, and store the number of genes separately
  only.in.bombus = setdiff(bombus.ortho.list, apis.ortho.list)
  number.only.in.bombus = length(only.in.bombus)
  
  # Determine the genes present in list2 but not list1, and store the number of genes separately
  only.in.apis = setdiff(apis.ortho.list, bombus.ortho.list)
  number.only.in.apis = length(only.in.apis)
  
  # Determine the number of genes not present in either list
  neither.list = total.gene.list.size - (number.both.lists + number.only.in.bombus + number.only.in.apis)
  
  # Create a matrix representing the numbers of genes in both lists, only list1, only list2, and neither list
  contingency.table = matrix(c(number.both.lists, number.only.in.bombus, number.only.in.apis, neither.list))
  # Change the dimensions to 2 rows and 2 columns
  dim(contingency.table) <- c(2,2)
  
  # Print the matrix 
  print("Contingency table")
  print(contingency.table)
  
  # Conduct two-tailed Fisher's exact test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance
  print(fisher.test(contingency.table))
}

# Run comparisons between gene lists with the compare.bt.amel.lists function
# Form is comparelists("bombus list", "apis list")
# Gene list name needs to be a string


####Early instar significant (LFC>0) comparisons

#Queen
EQAllsig_006_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_006_queen.txt")
EQAllsig_012_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_012_queen.txt")
EQAllsig_036_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_036_queen.txt")
EQAllsig_060_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_060_queen.txt")

#Worker
EWAllsig_006_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_006_worker.txt")
EWAllsig_012_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_012_worker.txt")
EWAllsig_036_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_036_worker.txt")
EWAllsig_060_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_060_worker.txt")


####Mid instar significant (LFC>0) comparisons

#Queen
MQAllsig_036_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_036_queen.txt")
MQAllsig_060_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_060_queen.txt")
MQAllsig_084_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_084_queen.txt")

#Worker
MWAllsig_036_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_036_worker.txt")
MWAllsig_060_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_060_worker.txt")
MWAllsig_084_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_084_worker.txt")


####Late instar significant (LFC>0) comparisons

#Queen
LQAllsig_060_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_060_queen.txt")
LQAllsig_084_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_084_queen.txt")
LQAllsig_108_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_108_queen.txt")
LQAllsig_132_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_132_queen.txt")

#Worker
LWAllsig_060_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_060_worker.txt")
LWAllsig_084_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_084_worker.txt")
LWAllsig_108_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_108_worker.txt")
LWAllsig_132_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_132_worker.txt")



####Mid instar highly significant (LFC>1) comparisons

#Queen
MQHighsig_036_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_036_queen.txt")
MQHighsig_060_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_060_queen.txt")
MQHighsig_084_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_084_queen.txt")

#Worker
MWHighsig_036_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_036_worker.txt")
MWHighsig_060_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_060_worker.txt")
MWHighsig_084_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_084_worker.txt")


####Late instar highly significant (LFC>1) comparisons

#Queen

LQHighsig_060_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_060_queen.txt")
LQHighsig_084_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_084_queen.txt")
LQHighsig_108_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_108_queen.txt")
LQHighsig_132_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_132_queen.txt")

#Worker

LWHighsig_060_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_060_worker.txt")
LWHighsig_084_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_084_worker.txt")
LWHighsig_108_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_108_worker.txt")
LWHighsig_132_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_132_worker.txt")




