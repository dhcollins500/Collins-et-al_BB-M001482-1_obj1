#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Generate Euler diagrams of the differentially expressed genes (DEGs) 
# comparison results.
#-------------------------------------------------------------------------------
# Inputs:
# .csv data files generated by script S4_40.

# Outputs:
# Two .svg figures (one for DEGs and one for HDEGs) of Euler diagrams.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(eulerr)  # euler()
library(ggplotify)  # as.ggplot()
library(ggplot2)
library(ggpubr)  # ggarrange()

# LOAD DATA ----

# No data loaded, as the function does it.

# FUNCTION DEFINITIONS ----

MakeEulerPlot <- function(x, y = "DEGs") {
  # Make a Euler plot of the results of comparing the DEGs lists from the 
  # current study, Kallisto pipeline and HISAT2/HTSeq pipeline.
  #
  # Args:
  #   x = string denoting time point and caste ("EW", "EQ", "MW", "MQ", 
  #       "LW" or "LQ").
  #   y = string denoting level of differentiated gene ("DEGs" (default)
  #       or "HDEGs").
  #
  # Returns:
  #   Euler diagram of the results.
  
  # Load relevant csv file based on arguments.
  
  resultsDF <- read.csv(paste0(x, "_larvae_", 
                               y, "_euler_data.csv"))
  
  # Fit euler diagram to data.
  
  eulerData <- euler(c(current_study = resultsDF[1, 2], HISAT2_HTSeq = resultsDF[2, 2], 
                       Kallisto = resultsDF[3, 2],
                       "current_study&Kallisto" = resultsDF[4,2],
                       "current_study&HISAT2_HTSeq" = resultsDF[5, 2],
                       "HISAT2_HTSeq&Kallisto" = resultsDF[6, 2],
                       "current_study&HISAT2_HTSeq&Kallisto" = resultsDF[7, 2]))
  
  # Plot euler diagram.
  
  ePlot <- as.ggplot(plot(eulerData, quantities = TRUE))
  
  # Return the plot.
  
  return(ePlot)
  
}

# EXECUTED STATEMENTS ----

# Set Working Directory ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

setwd("../02_outputs/70_euler_diagram_data")

# Make Plots for DEGs ----

ewPlot <- MakeEulerPlot("EW")

eqPlot <- MakeEulerPlot("EQ")

mwPlot <- MakeEulerPlot("MW")

mqPlot <- MakeEulerPlot("MQ")

lwPlot <- MakeEulerPlot("LW")

lqPlot <- MakeEulerPlot("LQ")

# Make Plots for HDEGs ----

mwPlotHDEGs <- MakeEulerPlot("MW", "HDEGs")

mqPlotHDEGs <- MakeEulerPlot("MQ", "HDEGs")

lwPlotHDEGs <- MakeEulerPlot("LW", "HDEGs")

lqPlotHDEGs <- MakeEulerPlot("LQ", "HDEGs")

# Combine the Plots into Figures ----

# DEGs.

degsCombined <- ggarrange(ewPlot, eqPlot, mwPlot, mqPlot, lwPlot, lqPlot,
                          labels = c("A", "B", "C", "D", "E", "F"),
                          ncol = 2, nrow = 3)

# HDEGs.

hdegsCombined <- ggarrange(mwPlotHDEGs, mqPlotHDEGs, lwPlotHDEGs, lqPlotHDEGs,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2)

# Save figures as SVG files so that the positioning of the labels
# can be adjusted manually for clarity.

setwd("../80_comparison_results/")

ggsave("supplementary_figure_S9_degs_overlap_euler_diagram.svg", 
       plot = degsCombined, height = 29.7, width = 21, units = "cm")

ggsave("supplementary_figure_S10_hdegs_overlap_euler_diagram.svg", 
       plot = hdegsCombined, height = 29.7, width = 21, units = "cm")
