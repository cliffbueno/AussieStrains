# Select 2023 NEON Samples to use

#### 1. Setup #####
# Libraries
library(tidyverse)
library(ggside)
library(ozmaps)
library(sf)
library(scales)
library(arules)

# Functions
`%notin%` <- Negate(`%in%`)

# Working directory
setwd("~/Documents/GitHub/AussieStrains/")

#### 2. Explore ####
# Metadata sent by Hugh Cross Jan 29th, 2025 (wasn't publicly available yet)
# Total without filtering: 263 from 27 sites
# Non-wetland: 242 from 27 sites
# Non-wetland and <= 0.10: 25 from 9 sites
# Not-wetland and < 0.11
d <- read.delim("data/neon_soil_metagenome_samples_2023.tsv") %>%
  drop_na(Genome.Size....assembled) %>% # Remove those not sequenced yet
  separate(depth..meters, remove = F, sep = " - ", into = c("topDepth", "botDepth")) %>%
  filter(ecosystem_subtype != "Wetlands") %>%
  filter(soil.horizon == "O horizon") %>%
  filter(botDepth <= 0.11) %>%
  mutate(Site = substr(sample.name, start = 0, stop = 4))
names(d)
table(d$depth..meters)
table(d$Site)
table(d$ecosystem_subtype)
hist(d$Genome.Size....assembled)
