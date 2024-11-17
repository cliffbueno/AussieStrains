#!/bin/bash

# Path to the file with genome accessions
accession_file="../ncbi_accessions_batch4.txt"

# Directory to save downloaded genomes
output_dir="."

# Loop through each accession in the file
while IFS=$'\t' read -r accession; do
    
    # Download the genome data for the accession
    datasets download genome accession "$accession" --filename "$output_dir/$accession.zip" --include genome
    unzip "$accession.zip"
    cd "ncbi_dataset/data/$accession" 
    mv *.fna ../../../
    cd "../../../"
    rm -r "$accession.zip" ncbi_dataset *.md *.txt  
    
done < "$accession_file"q