#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Gather previously generated differentially expressed gene (DEG) lists,
# background gene lists, and gene conversion lists between species into one 
# directory in preparation for comparative analysis.
#-------------------------------------------------------------------------------
# Inputs:
# DEGs lists from the current study, Cameron et al. (2013), He et al. (2017) and 
# Amsalem et al. (2015). Background gene lists from the current study, Cameron et 
# al. (2013) and He et al. (2017). Bombus to Apis conversion gene list
# (Table S7).

# Outputs:
# Directory containing all the lists.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# No data loaded.

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Gather the Necessary Gene Lists for the Comparison
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory created by script S5_01, and that the scripts for S4 and S5 have
# been run to generate the necessary gene lists.

# Make a directory to store gene lists.

dir.create("../02_outputs/90_diapause_comparison_DEG_lists")
# Will produce a warning if directory already exists.

# Copy the relevant lists into the new directory.

# Current study HISAT2/HTSeq DEGs lists.

# Set working directory.

setwd("../../S4_mRNA_seq_comparison/02_outputs/60_HISAT2_HTSeq_DEG_lists")

# Copy files to new directory.

file.copy(dir(), 
          "../../../S5_analyses/02_outputs/90_diapause_comparison_DEG_lists")

# Cameron et al. (2013) DEGs lists.

setwd("../../../S5_analyses/02_outputs/50_Cameron_DEG_lists/")

# Copy files to new directory.

file.copy(dir(), "../90_diapause_comparison_DEG_lists/")

# He et al. (2017) DEGs lists.

setwd("../60_He_DEG_lists/")

# Copy files to new directory.

file.copy(dir(), "../90_diapause_comparison_DEG_lists/")

# Amsalem et al. (20150 diapause DEGs lists.

setwd("../80_diapause_DEG_lists/")

# Copy files to new directory.

file.copy(dir(), "../90_diapause_comparison_DEG_lists/")

# Bombus to Apis conversion gene list (Table S7).

setwd("../../00_data/20_table_s7_csv")

# Copy file to new directory.

file.copy("Table_S7.csv", "../../02_outputs/90_diapause_comparison_DEG_lists")
