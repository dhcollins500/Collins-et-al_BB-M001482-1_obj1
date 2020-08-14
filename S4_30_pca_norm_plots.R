#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. 2020
# Analysis: RNA-Seq method comparison
# Tasks: Produce principal component analysis and normalization boxplot
# from Kallisto and HISAT2/HTSeq data.
#-------------------------------------------------------------------------------
# Inputs:
# ddsKallisto and ddsHTSeq objects.

# Outputs:
# .svg figure of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggplot2)  # ggplot() and associated functions.
library(reshape2)  # melt()
library(ggpubr)  # ggarrange(), annotate_figure()
# DESeq2 and tximport loaded in the scripts in the "Data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Kallisto Data ----

source("S4_11_kallisto_tximport_DESeq2.R")

# HISAT2/HTSeq Data ----

setwd("../../01_scripts")

source("S4_22_htseq_DESeq2.R")
sessionInfo()

# FUNCTION DEFINITIONS ----

makeFigure <- function(x) {
  # Make a figure consisting of two plots: a PCA and a normalization
  # boxplot.
  #
  # Args:
  #   x = string denoting whether ddsKallisto or ddsHTSeq should be used
  #       to make the plots ("Kallisto" or "HTSeq").
  #
  # Returns:
  #   A .svg file of the two plots together as a figure
  
  # Set variables based on argument.
  
  if (x == "Kallisto") {
    dds <- ddsKallisto
    title <- "Kallisto pipeline"
    nameOfFile <- "pca_norm_figure_Kallisto.svg"
  } else if (x == "HTSeq") {
    dds <- ddsHTSeq
    title <- "HISAT2/HTSeq pipeline"
    nameOfFile <- "pca_norm_figure_HISAT2_HTSeq.svg"
  } else {
    stop('Argument x must be either "Kallisto" or "HTSeq".')
  }
  
  # Transform the data to stabilize the variance across the mean.
  # rlog transformation chosen over vst because there are 
  # relatively few samples (n < 30).
  
  ddsRlog <- rlog(dds)
  
  # Box plot of normalized, transformed counts.
  
  # Format data for plot.
  
  boxPlotData <- melt(assay(ddsRlog))
  boxPlotData$Var3 <- substr(boxPlotData$Var2, 0, 5)
  boxPlotData$Var2 <- NULL
  colnames(boxPlotData) <- c("gene", "rlog_transformed_counts", "sample")
  boxPlotData <- boxPlotData[, c("gene", "sample", "rlog_transformed_counts")]
  
  # Calculate median value of counts
  
  samplesMedian <- median(boxPlotData$rlog_transformed_counts)
  
  # Generate box plot.
  
  boxPlot <- ggplot(boxPlotData, aes(x = sample, y = rlog_transformed_counts)) + 
    geom_boxplot() +
    geom_hline(yintercept = samplesMedian) +  # Add median line
    ylab("rlog_transformed_counts")        
  
  # Format box plot.
  
  formattedBoxPlot <- boxPlot + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5))
  
  # Principle component analysis (PCA) of samples
  # (uses top 2000 genes by row variance).
  
  # Generate PCA (using 2000 genes with highest variance).
  
  pcaPlot <- plotPCA(ddsRlog, intgroup = "condition",
                     ntop = 2000)
  
  # Format PCA.
  
  formattedPCAPlot <- pcaPlot + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Group plots into a single plot and annotate with method name.
  
  figure <- ggarrange(formattedBoxPlot, formattedPCAPlot,
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
  
  figure <- annotate_figure(figure, top = title)
  
  # Save figure.
  
  ggsave(nameOfFile, figure, width = 29.7, height = 21, units = "cm")
  
}

# EXECUTED STATEMENTS ----

# Make a Results Directory ----

dir.create("../02_outputs/40_pca_norm_plots")
# Will produce a warning if directory already exists.

setwd("../02_outputs/40_pca_norm_plots")

# Make the Figures ----

makeFigure("Kallisto")

makeFigure("HTSeq")