#!/bin/bash

# Path to the file with focal genus genome accessions
accession_file="../brady_accessions.csv"

# Directory with genome fastas
input_dir="./Ref_genomes"

# Directory to save downloaded genomes
output_dir="./brady"

# Loop through each accession in the file
while IFS=, read -r accession; do
    
    # Copy the genomes of that genus into its own folder for coverm
    cp "$input_dir/$accession" "$output_dir/$accession"
    
done < "$accession_file"