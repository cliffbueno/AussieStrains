# Select LUCAS Samples to use

# At each LUCAS point, there are N, S, E, W subsamples taken
# Soils are 0 - 20 cm
# n = 885 from 2018 for
# n = 658 actually on NCBI

# Workflow:
# Downloaded "Soil", "0-0.1m depth", "Conservation and Natural Environments"
# This yields 376 samples
# 346 have define upland vegetation types (no mangrove, marsh, other, NA)
# 339 have AI (aridity index) data
# 331 have > 10 million reads and < 40 million reads sequencing depth
# 104 spreading aridity gradient. 4 bins with 26 MG each. 
# Chose 26 samples with greatest depth for 3 arid bins, used all 26 of the least arid bin.

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

d <- read.csv("~/Desktop/Fierer/Strains/LUCAS/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018.csv") %>%
  filter(Depth == "0-20 cm") %>% # 0-20 cm
  filter(Elev >= 0 & Elev <= 5642) %>% # Above sea level and below highest point in Europe
  filter(LC0_Desc %notin% c("Cropland", "Artificial land", "Water", "Wetlands")) %>%
  filter(LU1_Desc == "Semi-natural and natural areas not in use") %>%
  filter(LC1_Desc %notin% c("Lichens and Moss", "Other bare soil", "Rocks and stones",
                            "Sand"))
names(d)
table(d$LC0_Desc)
table(d$LC1_Desc)
table(d$LU1_Desc)
