#!/bin/bash

    #metadata
    metadata=../ncbi_accessions_batch2.csv
    #
    Red="$(tput setaf 1)"
    Green="$(tput setaf 2)"
    Bold=$(tput bold)
    reset=`tput sgr0` # turns off all atribute
    while IFS=, read -r field1  

    do 
        datasets download genome accession "${field1}" --filename "${field1}".zip --include genome
    done < ${metadata}