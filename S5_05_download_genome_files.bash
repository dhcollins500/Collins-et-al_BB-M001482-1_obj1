#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: Supplementary file S5
# Subsection: 1-5 (All analyses)
# Tasks: Download Apis mellifer (Amel) and Bombus terrestris (Bter) genome 
# files for analyses.
#-------------------------------------------------------------------------------
# Inputs:
# None.
#
# Outputs:
# Amel genome fasta file and GTF file.
# Bter genome fasta file and GTF file.
# CHECKSUMS files for the downloaded files.
#-------------------------------------------------------------------------------

# STEP 1: DOWNLOAD AND UNZIP AMEL GENOMNE FASTA FILE ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S5_01.

# Change working directory to correct data subdirectory.

cd ../00_data/01_fasta

# Download genome file.

# Args:
# -O: write file to following file path. 

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-48/fasta/apis_mellifera/dna/Apis_mellifera.Amel_4.5.dna.toplevel.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-48/fasta/apis_mellifera/dna/CHECKSUMS -O ../03_sums/Amel_dna_toplevel_checksums.txt

# Unzip genome file.

gunzip --stdout Apis_mellifera.Amel_4.5.dna.toplevel.fa.gz > Amel_v4_5_genome_DNA.fasta

# STEP 2: DOWNLOAD AND UNZIP BTER GENOME FILE ----

# Download Bombus terrestris fasta genome for HISAT2/HTSeq pipeline.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz

# Download associated CHECKSUMS.
# Args:
# -O: Name of output file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/CHECKSUMS -O ../03_sums/Bter_dna_toplevel_checksums.txt

# Unzip genome fasta file.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz > Bter_toplevel_dna.fasta

# STEP 3: DOWNLOAD and UNZIP AMEL GTF FILE ----

# Change working directory.

cd ../02_gtf

# Download GTF file.

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-48/gtf/apis_mellifera/Apis_mellifera.Amel_4.5.48.gtf.gz

wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-48/gtf/apis_mellifera/CHECKSUMS -O ../03_sums/Amel_gtf_checksums.txt

# Unzip GTF file.

gunzip --stdout Apis_mellifera.Amel_4.5.48.gtf.gz > Amel_v4_5_48.gtf

# STEP 4: BTER AND UNZIP GTF FILE ----

# Download Bombus terrestris GTF file for use with HTSeq.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/Bombus_terrestris.Bter_1.0.47.gtf.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/CHECKSUMS -O ../03_sums/Bter_gtf_checksums.txt

# Unzip GTF file.

gunzip --stdout Bombus_terrestris.Bter_1.0.47.gtf.gz > Bter.gtf
