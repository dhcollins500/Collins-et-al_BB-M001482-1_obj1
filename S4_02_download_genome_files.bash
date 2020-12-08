#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Download Bombus terrestris genome files for analysis.
#-------------------------------------------------------------------------------
# Inputs:
# None.
#
# Outputs:
# Bombus terrestris genome and transcriptome fasta files, and GTF file.
# CHECKSUMS files for the downloaded files.
#-------------------------------------------------------------------------------

# STEP 1: DOWNLOAD BOMBUS TERRESTRIS GENOME FILE ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Change working directory to correct data subdirectory.

cd ../00_data/01_fasta

# Download Bombus terrestris fasta genome for HISAT2/HTSeq pipeline.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz

# Download associated CHECKSUMS.
# Args:
# -O: Name of output file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/CHECKSUMS -O ../03_sums/Bter_dna_toplevel_checksums.txt

# STEP 2: BOMBUS TERRESTRIS TRANSCRIPTOME FILE ----

# Download Bombus terrestris transcriptome for Kallisto pipeline.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/cdna/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/cdna/CHECKSUMS -O ../03_sums/Bter_cdna_all_checksums.txt

# STEP 3: BOMBUS TERRESTRIS GTF FILE ----

# Change working directory to correct data subdirectory.

cd ../02_gtf

# Download Bombus terrestris GTF file for use with HTSeq.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/Bombus_terrestris.Bter_1.0.47.gtf.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/CHECKSUMS -O ../03_sums/Bter_gtf_checksums.txt
