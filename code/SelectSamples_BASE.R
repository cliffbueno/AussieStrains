# Select BASE Samples to use

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

#### 2. Explore ####
d <- read.csv("data/contextual_nounits.csv")
names(d)
hist(d$water_content)
hist(d$ph)
hist(d$organic_carbon)
hist(d$ammonium_nitrogen_wt)
hist(d$nitrate_nitrogen)
hist(d$total_nitrogen)
plot(d$phosphorus_colwell)
plot(d$potassium_colwell)
hist(d$clay)
hist(d$sand)
hist(d$silt)
hist(d$conductivity)
plot(d$dtpa_zinc)
plot(d$dtpa_copper)
plot(d$dtpa_iron)
plot(d$dtpa_manganese)
plot(table(d$vegetation_type))
table(d$vegetation_type)
plot(d$water_content, d$organic_carbon)



#### 3. Veg. Filter ####
veg <- c("Dune", "Forest", "Grassland", "Heathland", "Savannah", "Shrubland", "Woodland")
d <- read.csv("data/contextual_nounits.csv") %>%
  filter(vegetation_type %in% veg) %>%
  separate(sample_id, into = c("Prefix", "sampleID"), sep = "\\/")
table(d$vegetation_type)



#### 4. AI Filter ####
coords346 <- d %>%
  select(sampleID, longitude, latitude)
#write.csv(coords346, "coords346.csv")

# Per docs: Global-AI values need to be multiplied by 0.0001 to retrieve the values in the correct units.
# Note these are AI from version 3 database, Zomer et al. 2022
a <- read.csv("data/coords346_wAI.csv") %>%
  mutate(AI = ai_v3_yr * 0.0001) %>%
  mutate(Dataset = "BASE") %>%
  filter(is.na(et0_v3_yr) == FALSE) # 7 NA
plot(a$AI)
hist(a$AI)
pdf("InitialFigs/BASE_n339_AI.pdf", width = 7, height = 5)
ggplot(a, aes(x = Dataset, y = AI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = "Dataset", y = "Aridity index") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme_ggside_minimal()
dev.off()

coords_trans <- st_as_sf(a, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
ozmap()
sf_oz <- ozmap("states")
pdf("InitialFigs/Mapn339.pdf", width = 7, height = 5)
ggplot(data = sf_oz) + 
  geom_sf(fill = "grey80", color = "white") +
  geom_sf(data = coords_trans,
          aes(color = AI)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal()
dev.off()



#### 5. Seq Filter ####
# Read in files for every sample
setwd("data/ND17/")
list.files()
length(list.files())
stat.files <- list()
for (i in 1:length(list.files())) {
  dir_name <- list.files()[i]
  setwd(dir_name)
  dir_name2 <- list.dirs()[2]
  setwd(dir_name2)
  file_name <- list.files()[3]
  stat.files[[i]] <- read.delim(file_name)
  setwd("~/Documents/GitHub/AussieStrains/data/ND17/")
}

# Reset wd
setwd("~/Documents/GitHub/AussieStrains/")

# Combine files into one dataframe
seq_depth <- stat.files[[1]]
for (i in 2:length(stat.files)) {
  seq_depth <- rbind(seq_depth, stat.files[[i]])
}

# Get R1 and R2 and check they're the same
seq_depth <- seq_depth %>%
  filter(grepl("R1p|R2p", Sample))
length(unique(seq_depth$FastQC_mqc.generalstats.fastqc.total_sequences))

# Get sampleID, filter to R1 and 339
seq_depth <- seq_depth %>%
  filter(grepl("R1p", Sample)) %>%
  separate(Sample, into = c("sampleID", "Other1", "Other2"), remove = F, sep = "_") %>%
  dplyr::select(-Other1, -Other2) %>%
  filter(sampleID %in% a$sampleID) %>%
  set_names(gsub("FastQC_mqc.generalstats.fastqc.", "", names(.))) %>%
  mutate(Dataset = "BASE")

ggplot(seq_depth, aes(x = Dataset, y = total_sequences)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = "Dataset", y = "# Reads") +
  scale_y_continuous(labels = label_comma()) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme_ggside_minimal()

# Check 4 outliers (8160, 8158, 8154, 8268) with insane depth and remove
d4 <- d %>%
  filter(sampleID %in% c("8160", "8158", "8154", "8268"))
# 3 grassland and 1 woodland
a4 <- a %>%
  filter(sampleID %in% c("8160", "8158", "8154", "8268"))
# 3 in central and 1 in southeast

seq_depth <- seq_depth %>%
  filter(sampleID %notin% c("8160", "8158", "8154", "8268"))

ggplot(seq_depth, aes(x = Dataset, y = total_sequences)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = "Dataset", y = "# Reads") +
  scale_y_continuous(labels = label_comma()) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme_ggside_minimal()

# Remove 4 below 10 million
seq_depth <- seq_depth %>%
  filter(total_sequences > 10000000)

pdf("InitialFigs/BASE_n331_SeqDepth.pdf", width = 7, height = 5)
ggplot(seq_depth, aes(x = Dataset, y = total_sequences)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = "Dataset", y = "# Reads") +
  scale_x_discrete(labels = c("BASE Soil Metagenomes\n0 - 10 cm\nNatural and Conservation Areas\nn = 331")) +
  scale_y_continuous(labels = label_comma()) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme_ggside_void()
dev.off()



#### 6. Final Selection ####
# Divide up AI gradient
# Within AI bins, select the metagenomes with the highest sequencing depth
# AIM for around 100 samples
a_toMerge <- a %>%
  dplyr::select(sampleID, bio12, bio1, et0_v3_yr, ai_v3_yr, AI) %>%
  mutate(sampleID = as.character(sampleID))
s_toMerge <- seq_depth %>%
  dplyr::select(sampleID, percent_duplicates, percent_gc, avg_sequence_length,
                percent_fails, total_sequences) %>%
  mutate(sampleID = as.character(sampleID))
d <- d %>%
  filter(sampleID %in% s_toMerge$sampleID) %>%
  left_join(., a_toMerge, by = "sampleID") %>%
  left_join(., s_toMerge, by = "sampleID")
# Save these 331 metadata!
#write.csv(d, "data/metadata_331.csv")

# Make 10 AI bins
max(d$AI)
min(d$AI)
table(discretize(d$AI, breaks = 13))
d$AI_bin <- cut(d$AI,
                breaks = c(0.0576, 0.483, 0.635, 0.943, 2.55),
                labels = c('Bin1', 'Bin2', 'Bin3', 'Bin4'))
table(d$AI_bin)

# For each of the first 3 AI bin, select the 26 samples with the highest seq depth
# For the wet bin, use all 26
d_sub <- d %>%
  arrange(AI_bin, desc(total_sequences)) %>%
  group_by(AI_bin) %>%
  slice_head(n = 26)
table(d_sub$AI_bin)

check <- d_sub %>%
  dplyr::select(sampleID, AI_bin, AI, latitude, longitude, vegetation_type, total_sequences)

plot(d_sub$AI_bin, d_sub$total_sequences)
plot(d_sub$AI, d_sub$total_sequences)

# Make sure veg. and geography are pretty well distributed
table(d_sub$vegetation_type)

coords_trans <- st_as_sf(d_sub, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
ozmap()
sf_oz <- ozmap("states")
pdf("InitialFigs/Mapn104.pdf", width = 7, height = 5)
ggplot(data = sf_oz) + 
  geom_sf(fill = "grey80", color = "white") +
  geom_sf(data = coords_trans,
          aes(color = AI)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal()
dev.off()

d_sub$Dataset <- "BASE"
pdf("InitialFigs/BASE_n104_SeqDepth_AI.pdf", width = 7, height = 5)
ggplot(d_sub, aes(x = Dataset, y = total_sequences)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75, aes(colour = AI)) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = NULL, y = "# Reads") +
  scale_color_gradient(low = "red", high = "blue") +
  scale_x_discrete(labels = c("BASE Soil Metagenomes\n0 - 10 cm\nNatural and Conservation Areas\nn = 104")) +
  scale_y_continuous(labels = label_comma()) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme_ggside_void()
dev.off()

plot(d_sub$AI)
hist(d_sub$AI)
pdf("InitialFigs/BASE_n104_AI.pdf", width = 7, height = 5)
ggplot(d_sub, aes(x = Dataset, y = AI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  geom_ysidedensity(aes(x = after_stat(density))) +
  labs(x = NULL, y = "Aridity index") +
  scale_x_discrete(labels = c("BASE Soil Metagenomes\n0 - 10 cm\nNatural and Conservation Areas\nn = 104")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme_ggside_void()
dev.off()

#write.csv(d_sub, "data/metadata104.csv")

sort(as.numeric(d_sub$sampleID))

# Go to BASE and select these 104 sampleIDs
# Download metagenomic data
# Download 55 metagenomic files for 104 samples
# Note: BASE uses the SqueezeMeta pipeline
# Data requested 9/23/24
# Thank you for your data request
# Your request id is 18, lodged at 2024-09-24T00:36:15.176Z
# We will be in touch. Please contact am-data-requests@bioplatforms.com for more information.

# Feb 1 24, decided we want all 331 of these! Send sampleID list to Matt Smith
# But since we already have 104, remove those, send list of 227
table(d$vegetation_type)
d_227 <- d %>%
  filter(sampleID %notin% d_sub$sampleID) %>%
  dplyr::select(sampleID)
writexl::write_xlsx(d_227, "~/Desktop/Cliff_BASE_samples227.xlsx")
