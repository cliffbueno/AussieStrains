#!/bin/bash

# Path to the file with genome accessions
accession_file="/scratch/alpine/clbd1748/Australia/dRep_filenames.txt"

# Define source and destination directories
source_dir="/scratch/alpine/clbd1748/Australia/Ref_genomes"
destination_dir="/scratch/alpine/clbd1748/Australia/Ref_genomes_dRep"

# Loop through each line in the text file (assuming each line is a filename)
while IFS=, read -r filename; do
  # Check if the file exists in the source directory
  if [ -e "$source_dir/$filename" ]; then
    # Copy the file to the destination directory
    cp "$source_dir/$filename" "$destination_dir/"
  else
    echo "File $filename does not exist in $source_dir"
  fi
done < "$accession_file"
