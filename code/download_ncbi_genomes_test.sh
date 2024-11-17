#!/bin/bash

# Path to the file with genome accessions
accession_file="../ncbi_accessions_batch2.csv"

# Directory to save downloaded genomes
output_dir="."

# Loop through each accession in the file
while IFS=, read -r accession; do
    
    # Download the genome data for the accession
    datasets download genome accession "$accession" --filename "$output_dir/$accession.zip" --include genome
    
done < "$accession_file"