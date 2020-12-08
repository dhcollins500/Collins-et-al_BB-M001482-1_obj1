#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Produce principal component analysis and normalization boxplot
# from Cameron et al. (2013), He et al. (2017) and Amsalem et al. (2015)
# HISAT2/HTSeq data.
#-------------------------------------------------------------------------------
# Inputs:
# ddsCameron, ddsHe and ddsAmsalem objects.

# Outputs:
# .svg figure of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggplot2)  # ggplot() and associated functions.
library(reshape2)  # melt()
library(ggpubr)  # ggarrange(), annotate_figure()
# DESeq2 loaded in the script in the "Data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Cameron et al. (2013), He et al. (2017) and Amsalem et al. (2015) Data ----

source("S5_50_HTSeq_to_DESeq2.R")
sessionInfo()

# FUNCTION DEFINITIONS ----

makeFigure <- function(x) {
  # Make a figure consisting of two plots: a PCA and a normalization
  # boxplot.
  #
  # Args:
  #   x: string denoting whether ddsCameron, ddsHe or ddsAmsalem should be used
  #      to make the plots ("cameron", "he" or "amsalem").
  #
  # Returns:
  #   A .svg file of the two plots together as a figure.
  
  # Set variables based on argument.
  
  if (x == "cameron") {
    dds <- ddsCameron
    title <- "Cameron et al. (2013)"
    nameOfFile <- "pca_norm_figure_Cameron.svg"
  } else if (x == "he") {
    dds <- ddsHe
    title <- "He et al. (2017)"
    nameOfFile <- "pca_norm_figure_He.svg"
  } else if (x == "amsalem") {
    dds <- ddsAmsalem
    title <- "Amsalem et al. (2015)"
    nameOfFile <- "pca_norm_figure_Amsalem.svg"
  } else {
    stop('Argument x must be either "cameron","he" or "amsalem".')
  }
  
  # Transform the data to stabilize the variance across the mean.
  # rlog transformation chosen over vst because there are 
  # relatively few samples (n < 30).
  
  ddsRlog <- rlog(dds)
  
  # Box plot of normalized, transformed counts.
  
  # Format data for plot.
  
  boxPlotData <- melt(assay(ddsRlog))
  colnames(boxPlotData) <- c("gene", "sample", "rlog_transformed_counts")
  
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
            # + geom_text(aes(label = name), vjust = 2, check_overlap = TRUE, size = 4)
  # Commented out code adds sample names to the PCA.
  
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

makeFigure("cameron")

makeFigure("he")

makeFigure("amsalem")

# Repeating Amsalem et al. (2015) Analysis ----

# Inspection of the Amsalem et al. (2015) PCA shows that sample D_line141 
# clusters with the mated and founding post-diapause samples rather than the
# other diapause samples (this can be seen easier by uncommenting the code in 
# line 104 above and re-running line 143).
# Therefore, this sample needs to be excluded and the analysis re-run.

# Deleting the HTSeq file for D_line141.

file.remove("../32_amsalem_htseq_count/D_line141_HTSeq.txt")

# Rename the current Amsalem et al. (2015) PCA plot so that it is not 
# overwritten.

file.rename("pca_norm_figure_Amsalem.svg", 
            "pca_norm_figure_Amsalem_all_samples.svg")

# Remove dds objects for re-run of analysis.

rm(list = c("ddsAmsalem", "ddsCameron", "ddsHe"))

# Load data again.

source("../../01_scripts/S5_50_HTSeq_to_DESeq2.R")

# Run makeFigure again.

makeFigure("amsalem")
