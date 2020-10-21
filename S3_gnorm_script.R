#------------------------------------------------------------------
# Author: David Prince, David Collins
# File started: 27.03.18
# File last updated: 16.07.20
# Tasks: Transform raw Cq values for candidate reference genes 
# into relative quantities (using the primer pair efficiency), conduct 
# stability analysis of expression of candidate reference genes using 
# geNORM plot the results and output the graphs. This script has been 
# modified by DC to work with his data
#------------------------------------------------------------------
# Inputs:
# Quality controlled raw Cq values for candidate reference genes.
#
# Outputs:
# Number of genes required for normalization and which genes are 
# most stably expressed, graphs showing the data
#------------------------------------------------------------------

#Reset R
rm(list=ls())

# LOADING PACKAGES AND DATA ----

# Packages ----

library(SLqPCR)  # relQuantPCR(), selectHKgenes()
library(ggplot2)  # ggplot
library(svglite)  # ggsave for svg files
sessionInfo()

# Loading Data ----

# Data are the means of the technical replicates of the quality 
# controlled, raw Cq values

# setwd("") fill in as appropriate
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/qRTPCRdata/normQBterdata/2020_analysis/")
cqsPostQC <- read.csv("ref_gene_mean_larvae_cqs.csv")

ampEff <- read.csv("ref_efficiency.csv")

# FUNCTION DEFINITIONS ----

# Plotting Stability of Expression Graphs

PlotStabilityGraph <- function(w, x, y, z, NF = FALSE) {
  # Plots a graph of the stability results against the gene name
  #
  # Args:
  #   w: The object from which the data is to be plotted
  #   x: w$x where x = the column name for the data for the x-axis 
  #      (the gene names in object w)
  #   y: w$y where y = the column name for the data for the y-axis 
  #      (the stability values in object w)
  #   z: A string containing the name of the y-axis label
  #   NF: Whether the data is from NormFinder
  # Returns:
  #   A plot of x vs. y from object w
  if (NF == TRUE) {
    newPlot <-
      ggplot(w, aes_(x = x, y = y)) + 
      geom_line(aes(group = 1)) +
      geom_point(size = 3, shape = 19) + 
      scale_x_discrete(limits = rev(x)) +
      xlab("Genes expressed with increasing stability") + 
      ylab(z)
  } else {
    newPlot <-
      ggplot(w, aes_(x = x, y = y)) + 
      geom_line(aes(group = 1)) +
      geom_point(size = 3, shape = 19) + 
      scale_x_discrete(limits = x) +
      xlab("Genes expressed with increasing stability") + 
      ylab(z)
  }
  newPlot <- 
    newPlot + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  return(newPlot)
}

# EXECUTED STATEMENTS ----

# Generating Relative Quantities from the Raw Data ----

# Rationale:
# It's important to take into account the primer pair efficiency in the
# candidate reference gene stability analysis (De Spiegelaere et al 2015, 
# PLoS One). Therefore the Cqs need to be converted into 
# relative quantities (RQs), with the formula:
# RQ = E^-(minCq-sampleCq)
# where E is the amplification efficiency of the primer pair:
# E = 10^(-1/slope)

# Select the first eight genes from the data frame, as these are the 
# candidate reference genes.

ampEff <- ampEff[1:8, ]

# Add an extra column where slope is used to calculate "E".

ampEff$E <- 10^(-1/ampEff$Slope)

# Converting to RQs
# RQ = E^(-(minCq-sampleCq))
# relQuantPCR in SLqPCR returns the same results as doing the RQ
# calculation by hand, and the code shows that it uses the 
# equation RQ = E^(-(minCq-sampleCq))

rqData <- data.frame(relQuantPCR(cqsPostQC$PepIsoA, E = ampEff$E[1]))
names(rqData) <- "PepisoA"
rqData$GAPDH2 <- (relQuantPCR(cqsPostQC$GAPDH2, E = ampEff$E[2]))
rqData$PP2A <- (relQuantPCR(cqsPostQC$PP2A, E = ampEff$E[3]))
rqData$TATA_box_binding_protein <-
  (relQuantPCR(cqsPostQC$TATA_box_binding_protein, E = ampEff$E[4]))
rqData$RPS5a <- (relQuantPCR(cqsPostQC$RPS5a, E = ampEff$E[5]))
rqData$RPS18 <- (relQuantPCR(cqsPostQC$RPS18, E = ampEff$E[6]))
rqData$Arginine_kinase <-
  (relQuantPCR(cqsPostQC$Arginine_kinase, E = ampEff$E[7]))
rqData$Tubulin_alpha_1_chain <-
  (relQuantPCR(cqsPostQC$Tubulin_alpha_1_chain, E = ampEff$E[8]))

# Stability of Reference Genes Analysis ----

# geNORM analysis
# Data already in correct format for the analysis.

genormResults <- selectHKgenes(rqData, geneSymbol = names(rqData))

# geNORM ranking of the genes from most to least stable

rankedGenormResults <- genormResults$ranking
rankedGenormResults

# Summary of Results ----

# geNORM analysis # also shows that two reference genes are necessary 
# for downstream analysis. PepisoA and Arginine_kinase are the two most 
# stably expressed reference genes in these data.

# Plotting the Data ----

# Setting the theme for the graphs

theme_set(theme_bw(base_size=10, base_family="Times") + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black")))


# geNORM results
# M (average expression stability) vs. gene name
# geNORM results put into data frame to plot with ggplot2

genormStabilityResults = data.frame(genormResults$meanM, attributes(genormResults$meanM))

# Combine the first two entries in rankedGenormResults so that the names can 
# be used as a legend for the plot

rankedGenormResultsPlotting <- unname(rankedGenormResults) 
rankedGenormResultsPlotting[1] <- 
  paste(rankedGenormResultsPlotting[1], rankedGenormResultsPlotting[2], sep="\n")
rankedGenormResultsPlotting <- rankedGenormResultsPlotting[-2]  
rankedGenormResultsPlotting <- rev(rankedGenormResultsPlotting) 
rankedGenormResultsPlotting

# Update genormStabilityResults with the new names

genormStabilityResults$gene_names <- rankedGenormResultsPlotting 
names(genormStabilityResults) <- c("meanM", "rank_of_gene", "gene_names") 
genormStabilityResults

# Plotting the graph
plotGenormMData <- 
  PlotStabilityGraph(genormStabilityResults, 
                     genormStabilityResults$gene_names,
                     genormStabilityResults$meanM,
                     "Average expression stability (M)")
plotGenormMData

# V (pairwise variation) vs. number of genes
# Would like to get the data into data frames so i can plot with ggplot2

genormVariationResults <- 
  data.frame(genormResults$variation, attributes(genormResults$variation))

# Plotting the graph

plotGenormVData <- 
  ggplot(genormVariationResults, aes(x = names, y = genormResults.variation)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.15, linetype = "dashed", size = 1) + 
  scale_x_discrete(limits = rev(genormVariationResults$names)) +
  ylim(0, 0.25) +
  xlab("The number of reference genes being compared") + 
  ylab("Pairwise variation (V)")
plotGenormVData

# Exporting the plots for figure construction ----

# setwd('') # Change as required

ggsave("DC_genormVData.svg", plot = plotGenormVData, width = 10.5, 
       height = 9.9, units = "cm")
ggsave("DC_genormMData.svg", plot = plotGenormMData, width = 10.5, 
       height = 9.9, units = "cm")

# Alternatively you could use the ggpubr package to combine them into one plot with
# labels etc, but I like the flexibility of putting them in inkscape and arranging/
# labelling with more control.
