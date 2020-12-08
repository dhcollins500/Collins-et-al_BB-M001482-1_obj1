#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: Collins et al. (2020)
# Analysis: mRNA-seq method comparison
# Tasks: Select relevant entries from checksums files and compare with the
# checksums generated from the downloaded files.
#-------------------------------------------------------------------------------
# Inputs:
# Bombus terrestris fasta and annotation files, and associated checksum files, downloaded
# from Ensembl using script 02.
#
# Outputs:
# .txt file stating whether the checksums match between the downloaded files
# and the files present on Ensembl.
#-------------------------------------------------------------------------------

# STEP 1: SELECT RELEVANT CHECKSUMS ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory
# created by script S4_01.

# Print the relevant checksums to a new document so that they can be compared.

# Change Directory ----

cd ../00_data/03_sums

# Print the Relevant Checksums ----

# Genome checksums.

awk '/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz/' Bter_dna_toplevel_checksums.txt > original_ensembl_checksums.txt

# cDNA checksums. 

awk '/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz/' Bter_cdna_all_checksums.txt >> original_ensembl_checksums.txt

# GTF checksums.

awk '/Bombus_terrestris.Bter_1.0.47.gtf.gz/' Bter_gtf_checksums.txt >> original_ensembl_checksums.txt

# STEP 2: CONDUCT CHECKSUMS ON DOWNLOADED FILES ----
# Print the checksums results and the file name to a new document so that they can be compared.

# Genome checksums

echo "$(sum ../01_fasta/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz) $(echo Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz)" > downloaded_ensembl_checksums.txt

# cDNA checksums

echo "$(sum ../01_fasta/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz) $(echo Bombus_terrestris.Bter_1.0.cdna.all.fa.gz)" >> downloaded_ensembl_checksums.txt

# GTF checksums

echo "$(sum ../02_gtf/Bombus_terrestris.Bter_1.0.47.gtf.gz) $(echo Bombus_terrestris.Bter_1.0.47.gtf.gz)" >> downloaded_ensembl_checksums.txt

# STEP 3: COMPARE THE ORIGINAL AND DOWNLOADED CHECKSUMS ----
# Check that the checksums of the downloaded file match those provided with the file
# by comparing the files.

# Make subdirectory in outputs directory

mkdir ../../02_outputs/01_checksums

# Make empty file for results output

touch ../../02_outputs/01_checksums/checksums_ensembl_comparison_results.txt

# Compare results 

# grep arguments:
# -q: quiet, exits with zero status (mapped to true in the "if" statement) if any matches
# found, and exits with non-zero status (i.e. "false") if no matches found.
# -f fileName: obtains pattern(s) for comparison from fileName. 

while IFS= read -r FILE
do
if grep -q "$FILE" original_ensembl_checksums.txt
then
	RESULTS="$FILE \t Checksums match, therefore downloaded file is likely to be complete"
else 
	RESULTS="$FILE \t Checksums do not match, therefore downloaded file is likely to be incomplete"
fi
echo -e $RESULTS >> ../../02_outputs/01_checksums/checksums_ensembl_comparison_results.txt
done < "downloaded_ensembl_checksums.txt"
