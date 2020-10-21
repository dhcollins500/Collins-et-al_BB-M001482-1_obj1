#-----------------------------------------------------------------------------------
# David Collins
# Started 19/01/2018
# Last modified 22/01/2018
# This script is a modification of David Prince's script to compare the mRNAseq results from He et al 2017 with the results
# from this study.
#-----------------------------------------------------------------------------------

# DESCRIPTION

# Purpose: to calculate whether the overlap between two gene lists is statistically significant

## Situation - the two lists to be compared are from different species
### Necessary data conversion - convert the gene list from one of the species into a list of orthologs in the other species
### Necessary inputs: 
#### input1 = a list of genes from this study as .txt file - Files labelled for example "DC_EQ_AllSig.txt" (10 files)
#### input2 = a list of genes from the study to be compared to as a .txt file - Files labelled for example "He_2d_queen","He_2d_worker" (4 files)
#### input3 = the number of genes with orthologs in both species - "new_conversion_sheet_both_part2.csv"
#### input4 = a list of the total number of genes that appear in the study to be compared to as a .txt file -
####                                                        "He_et_al_2017_All_genes_count.csv"

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
setwd("U:/Documents/RAnalysis/comparative_analyses/BombusApisCaste/") 
# change this address as necessary
# can copy and paste addresses from Windows Explorer but need to turn the \ into /

# double check that the working directory is set to the correct place
getwd()

# The dataset is from mRNAseq, but includes several unique identifiers that are not comparable between studies (for their novel genes)
# We will only look at the genes that are present in both studies.
# I recieved the total number of genes sequenced by mRNAseq from the authors themselves
mRNAseq = read.csv("He_et_al_2017_All_genes_count.csv", header = TRUE) 

# Isolate the BeeBase column to find the names of the genes present in the mRNAseq
gene.names = as.character(mRNAseq[,1])

# Some genes might be represented twice so remove these
# Therefore identify the unique gene Ids (removing redundancy)
unique.gene.names = unique(gene.names)

# Calculate which of the unique.gene.names have an RBH BLAST ortholog in Bombus terrestris
# Input Bombus terrestris to Apis mellifera RBH BLAST conversion list
conversion.df = read.csv("new_conversion_sheet_both_part2.csv", header = TRUE) 

# The He et al data is annotated using the OGS_v2 on beebase (i.e. after 2014)
#therefore match to beebase ids, to do this remove the other Apis Ids.
conversion.df2 = conversion.df[,c(3,6)]
# and the instances where no RBH BLAST ortholog exists between OGS_v1_1 and Bter Gene_Ids
conversion.df3 = na.omit(conversion.df2)
conversion.df3 = unique(conversion.df3) #removes redundancy
# Calculate the overlap between the RBH BLAST ortholog list (the column BEEBASE in the 'new conversion sheet 2') 
#and the gene list (the values labelled 'unique.gene.names')
mRNAseq.RBH = intersect(conversion.df3$BEEBASE, unique.gene.names)
# count the number of ortholog pairs between Bombus terrestris and Apis mellifera
total.gene.list.size = length(mRNAseq.RBH)
#Note the total number of orthologs is 6677

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
  apis.present.in.df = conversion.df3[match(apis.list, conversion.df3$BEEBASE),]
  apis.present.in.df2 = na.omit(apis.present.in.df)
  count.removed.NAs.apis = length(apis.present.in.df$BEEBASE) - length(apis.present.in.df2$BEEBASE)
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
# Bombus lists are available in "Supplementary datafile 4", Apis lists are located in the supplementary material for 
# He et al. 2017. Insect Sci. 00, 1-11, each gene list has been split into separate queen and worker lists at each time point (2d and 4d), and then stored
#as separate text files
## Gene list name needs to be a string

#### Early instar significant comparisons

#Within caste
EWsig_2dw = compare.bt.amel.lists("DC_EW_AllSig.txt","He_2d_worker.txt")
EQsig_2dq = compare.bt.amel.lists("DC_EQ_AllSig.txt","He_2d_queen.txt")

#Between caste
EWsig_2dq = compare.bt.amel.lists("DC_EW_AllSig.txt","He_2d_queen.txt")
EQsig_2dw = compare.bt.amel.lists("DC_EQ_AllSig.txt","He_2d_worker.txt")

#Random non-dif list
EWLR_2dq = compare.bt.amel.lists("DC_EW_LR.txt","He_2d_queen.txt")
EQLR_2dw = compare.bt.amel.lists("DC_EQ_LR.txt","He_2d_worker.txt")


##### Mid instar significant comparisons

#Within caste
MWsig_2dw = compare.bt.amel.lists("DC_MW_AllSig.txt","He_2d_worker.txt")
MQsig_2dq = compare.bt.amel.lists("DC_MQ_AllSig.txt","He_2d_queen.txt")
MWsig_4dw = compare.bt.amel.lists("DC_MW_AllSig.txt","He_4d_worker.txt")
MQsig_4dq = compare.bt.amel.lists("DC_MQ_AllSig.txt","He_4d_queen.txt")

#Between caste
MWsig_2dq = compare.bt.amel.lists("DC_MW_AllSig.txt","He_2d_queen.txt")
MQsig_2dw = compare.bt.amel.lists("DC_MQ_AllSig.txt","He_2d_worker.txt")
MWsig_4dq = compare.bt.amel.lists("DC_MW_AllSig.txt","He_4d_queen.txt")
MQsig_4dw = compare.bt.amel.lists("DC_MQ_AllSig.txt","He_4d_worker.txt")

#Random non-dif list
MWLR_2dq = compare.bt.amel.lists("DC_MW_LR.txt","He_2d_queen.txt")
MQLR_2dw = compare.bt.amel.lists("DC_MQ_LR.txt","He_2d_worker.txt")
MWLR_4dq = compare.bt.amel.lists("DC_MW_LR.txt","He_4d_queen.txt")
MQLR_4dw = compare.bt.amel.lists("DC_MQ_LR.txt","He_4d_worker.txt")


#### Late instar significant comparisons

#Within caste
LWsig_2dw = compare.bt.amel.lists("DC_LW_AllSig.txt","He_2d_worker.txt")
LQsig_2dq = compare.bt.amel.lists("DC_LQ_AllSig.txt","He_2d_queen.txt")
LWsig_4dw = compare.bt.amel.lists("DC_LW_AllSig.txt","He_4d_worker.txt")
LQsig_4dq = compare.bt.amel.lists("DC_LQ_AllSig.txt","He_4d_queen.txt")

#Between caste
LWsig_2dq = compare.bt.amel.lists("DC_LW_AllSig.txt","He_2d_queen.txt")
LQsig_2dw = compare.bt.amel.lists("DC_LQ_AllSig.txt","He_2d_worker.txt")
LWsig_4dq = compare.bt.amel.lists("DC_LW_AllSig.txt","He_4d_queen.txt")
LQsig_4dw = compare.bt.amel.lists("DC_LQ_AllSig.txt","He_4d_worker.txt")

#Random non-dif list
LWLR_2dq = compare.bt.amel.lists("DC_LW_LR.txt","He_2d_queen.txt")
LQLR_2dw = compare.bt.amel.lists("DC_LQ_LR.txt","He_2d_worker.txt")
LWLR_4dq = compare.bt.amel.lists("DC_LW_LR.txt","He_4d_queen.txt")
LQLR_4dw = compare.bt.amel.lists("DC_LQ_LR.txt","He_4d_worker.txt")


#### Medium instar highly significant (LFC>1) comparisons 

#Within caste
MWHighsig_2dw = compare.bt.amel.lists("DC_MW_HighSig.txt","He_2d_worker.txt")
MQHighsig_2dq = compare.bt.amel.lists("DC_MQ_HighSig.txt","He_2d_queen.txt")
MWHighsig_4dw = compare.bt.amel.lists("DC_MW_HighSig.txt","He_4d_worker.txt")
MQHighsig_4dq = compare.bt.amel.lists("DC_MQ_HighSig.txt","He_4d_queen.txt")

#Between caste
MWHighsig_2dq = compare.bt.amel.lists("DC_MW_HighSig.txt","He_2d_queen.txt")
MQHighsig_2dw = compare.bt.amel.lists("DC_MQ_HighSig.txt","He_2d_worker.txt")
MWHighsig_4dq = compare.bt.amel.lists("DC_MW_HighSig.txt","He_4d_queen.txt")
MQHighsig_4dw = compare.bt.amel.lists("DC_MQ_HighSig.txt","He_4d_worker.txt")

#Random non-dif list
MWSR_2dq = compare.bt.amel.lists("DC_MW_SR.txt","He_2d_queen.txt")
MQSR_2dw = compare.bt.amel.lists("DC_MQ_SR.txt","He_2d_worker.txt")
MWSR_4dq = compare.bt.amel.lists("DC_MW_SR.txt","He_4d_queen.txt")
MQSR_4dw = compare.bt.amel.lists("DC_MQ_SR.txt","He_4d_worker.txt")


####Late instar highly significant (LFC>1) comparisons

#Within caste
LWHighsig_2dw = compare.bt.amel.lists("DC_LW_HighSig.txt","He_2d_worker.txt")
LQHighsig_2dq = compare.bt.amel.lists("DC_LQ_HighSig.txt","He_2d_queen.txt")
LWHighsig_4dw = compare.bt.amel.lists("DC_LW_HighSig.txt","He_4d_worker.txt")
LQHighsig_4dq = compare.bt.amel.lists("DC_LQ_HighSig.txt","He_4d_queen.txt")

#Between caste
LWHighsig_2dq = compare.bt.amel.lists("DC_LW_HighSig.txt","He_2d_queen.txt")
LQHighsig_2dw = compare.bt.amel.lists("DC_LQ_HighSig.txt","He_2d_worker.txt")
LWHighsig_4dq = compare.bt.amel.lists("DC_LW_HighSig.txt","He_4d_queen.txt")
LQHighsig_4dw = compare.bt.amel.lists("DC_LQ_HighSig.txt","He_4d_worker.txt")

#Random non-dif list
LWSR_2dq = compare.bt.amel.lists("DC_LW_SR.txt","He_2d_queen.txt")
LQSR_2dw = compare.bt.amel.lists("DC_LQ_SR.txt","He_2d_worker.txt")
LWSR_4dq = compare.bt.amel.lists("DC_LW_SR.txt","He_4d_queen.txt")
LQSR_4dw = compare.bt.amel.lists("DC_LQ_SR.txt","He_4d_worker.txt")
