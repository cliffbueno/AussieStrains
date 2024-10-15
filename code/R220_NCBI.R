# Download entire GTDB database

## Info:
# 596,859 total genomes
# Bacteria: 584,382. Archaea: 12,477
# 113104 species clusters (available for download)

# GTDB website only provides download of the representative species clusters sequences.
# Need to get others from NCBI
# Get NCBI names from GTDB metadata file
setwd("~/Documents/GitHub/AussieStrains/")

arcGT <- read.delim("~/Desktop/Strains/ar53_metadata_r220.tsv") # fast
bacGT <- read.delim("~/Desktop/Strains/bac120_metadata_r220.tsv") # takes a while to load

gtdb220_full <- rbind(arcGT, bacGT)
ncbi_names <- gtdb220_full$ncbi_assembly_name
#write.table(ncbi_names, "data/ncbi_names.txt")

# After deciding on focal taxa, just select those and get all of those genomes