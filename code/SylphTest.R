# Use Sylph to identify reference genomes of the focal taxa

setwd("~/Documents/GitHub/AussieStrains/")
library(readr)
library(dplyr)

#### _Test ####
# Test Sylph on 4 Australian metagenomes
# Test how many genomes it outputs based on different AI
# Try the profile command and the query command

# Profiling, with ANI > 95%
sp95 <- read_tsv("data/profiling.tsv")
sp95_gPerS <- sp95 %>%
  group_by(Sample_file) %>%
  summarise(n_genomes = n())
sp95_focal <- sp95 %>%
  filter(grepl("Acidothermus|Bacillus|Bradyrhizobium|Bryobacter|Haliangium|Mycobacterium|Solibacter|Udaeobacter",
               Contig_name))
# Bacillus, Brady, Udaeo found
# Bacillus and Brady different genomes in different samples

# Querying, with ANI > 90%
sq90 <- read_tsv("data/querying90.tsv")
sq90_gPerS <- sq90 %>%
  group_by(Sample_file) %>%
  summarise(n_genomes = n())
sq90_focal <- sq90 %>%
  filter(grepl("Acidothermus|Bacillus|Bradyrhizobium|Bryobacter|Haliangium|Mycobacterium|Solibacter|Udaeobacter", 
               Contig_name))
sq90_focal_count <- sq90_focal %>%
  group_by(Contig_name) %>%
  summarise(n_samples = n())

# Querying, with ANI > 85%
sq85 <- read_tsv("data/querying85.tsv")
sq85_gPerS <- sq85 %>%
  group_by(Sample_file) %>%
  summarise(n_genomes = n())
sq85_focal <- sq85 %>%
  filter(grepl("Acidothermus|Bacillus|Bradyrhizobium|Bryobacter|Haliangium|Mycobacterium|Solibacter|Udaeobacter", 
               Contig_name))
sq85_focal_count <- sq85_focal %>%
  group_by(Contig_name) %>%
  summarise(n_samples = n())

# Based on tests with the Zymo kit, probably want to use 90% minimum ANI for Sylph for identifying reference genomes.