# Australia Preliminary Analysis
# Check Metadata and Hannah's PhyloFlash Data

# BioClim
# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO12 = Annual Precipitation
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# Due to issues with packages, do this in QGIS

#### 1. Setup ####
library(plyr)
library(tidyverse)
library(ozmaps)   
library(sf)
library(mctoolsr)
library(vegan)
library(FSA)
library(plotly)
library(ggside)
#library(ggstatsplot)
library(cowplot)

# Working Directory
setwd("~/Documents/GitHub/AussieStrains/")

# Metadata
moist <- read.csv("data/MoistureIndexALA.csv") %>%
  select(decimalLatitude, 
         #decimalLongitude, 
         Moisture.Index...annual.mean..Bio28.) %>%
  rename(latitude = decimalLatitude,
         #longitude = decimalLongitude,
         MoistureIndex = Moisture.Index...annual.mean..Bio28.) %>%
  mutate(latjoin = round(latitude, digits = 4)) %>%
  select(-latitude)
tp <- read.csv("data/coords_wTP.csv") %>%
  mutate(sampleID = as.character(sampleID)) %>%
  mutate(latjoin = round(latitude, digits = 4)) %>%
  left_join(., moist, by = "latjoin") %>%
  rename(MAP_bc = MAP,
         MAT_bc = MAT) %>%
  select(sampleID, MAP_bc, MAT_bc, MoistureIndex)

d <- read.delim("data/first_50_samples_data.txt") %>%
  mutate(sampleID = substr(Name2, start = 1, stop = 5)) %>%
  mutate(sampleID = gsub("_", "", sampleID)) %>%
  mutate(sampleID2 = gsub("_R1_001.fastq.gz", "", Name2)) %>%
  mutate(sampleID2 = gsub("_R2_001.fastq.gz", "", sampleID2)) %>%
  mutate(sampleID2 = gsub("_R1.fastq.gz", "", sampleID2)) %>%
  mutate(sampleID2 = gsub("_R2.fastq.gz", "", sampleID2)) %>%
  mutate(sampleID2 = gsub("-", ".", sampleID2)) %>%
  mutate(sampleID2 = paste("Sample_", sampleID2, sep = "")) %>%
  left_join(., tp, by = "sampleID")

uniq <- d %>%
  group_by(sampleID) %>%
  slice_head(n = 1)

# uniqL12 <- d %>%
#   group_by(sampleID2) %>%
#   slice_head(n = 1) %>%
#   select(sampleID2, 9, 10, 13, 16:21, 23:80, 85, 86) %>%
#   rename(sampleID = sampleID2)
# write.table(uniqL12, "data/metadata.txt", sep = "\t", row.names = F)

# From Hannah
metadata <- read.table("data/50_random_samples.txt", header = TRUE, sep = "\t") %>%
  rename(ammonium_nitrogen = ammonium_nitrogen..mg.kg., 
         boron = boron_hot_cacl2..mg.kg., 
         conductivity = conductivity_dsm..ds.m., 
         depth = depth..m., 
         copper = dtpa_copper..mg.kg., 
         iron = dtpa_iron..mg.kg., 
         manganese = dtpa_manganese..mg.kg., 
         zinc = dtpa_zinc..mg.kg.,
         extractable_aluminium = exc_aluminium..meq.100g., 
         extractable_calcium = exc_calcium..meq.100g., 
         extractable_magnesium = exc_magnesium..meq.100g., 
         extractable_potassium = exc_potassium..meq.100g., 
         extractable_sodium = exc_sodium..meq.100g.,
         nitrate_nitrogen = nitrate_nitrogen..mg.kg.,
         ph_level_cacl2 = ph_level_cacl2..cacl2., 
         ph_level_h2o = ph_level_h2o..h2o., 
         phosphorus_colwell = phosphorus_colwell..mg.kg., 
         potassium_colwell = potassium_colwell..mg.kg.,
         slope_aspect = slope_aspect..direction.or.degrees..e.g...nw.or.315.., 
         sulfur = sulphur..mg.kg.) %>%
  # Create new variables
  mutate(Inorganic_Nitrogen = nitrate_nitrogen + ammonium_nitrogen, 
         texture_2 = clay_percent + silt_percent) %>%
  # Filter out duplicated sample 7061_2
  filter(sample_extraction_id != "7061_2") %>%
  mutate(SampleID = paste0("Sample_", gsub("102.100.100/", "", sample_id))) %>%
  mutate(rowname = SampleID) %>%
  column_to_rownames("rowname")

extra <- read.csv("data/extra_metadata_base_samples.csv") %>%
  filter(catalogID %in% uniq$sampleID) %>%
  arrange(`Aridity.index...annual.mean`) %>%
  mutate(SampleID = paste0("Sample_", catalogID)) %>%
  filter(SampleID %in% metadata$SampleID) %>% 
  select(catalogID, SampleID, WorldClim..Temperature...annual.mean, # units = degC*10, worldclim
         Precipitation...annual.mean, # CSIRO ecosystem sciences; annual mean precip in mm; https://spatial.ala.org.au/ws/layers/view/more/rain_ann
         Aridity.index...annual.mean, # CSIRO ecosystem sciences; annual mean AI, ratio of precp to pot  evap. ; https://spatial.ala.org.au/ws/layers/view/more/arid_mean
         NPP.Mean, # tonnes/ha/yr; Modis derived from CSIRO
         Carbon...organic # % soil organic carbon flux; from ASRIS - modeled mean store kg/ha in soil
  ) %>%
  rename(MAT = WorldClim..Temperature...annual.mean,
         MAP = Precipitation...annual.mean, 
         AridityIndex = Aridity.index...annual.mean,
         NPP = NPP.Mean, 
         organic_Carbon_extrametadata = Carbon...organic) %>%
  mutate(MAT = MAT/10) # Correct WorldClim temperature units
head(extra$catalogID)
tail(extra$catalogID)
median(extra$AridityIndex)

# Huge outlier in carbon.
ggplot(extra, aes(AridityIndex, organic_Carbon_extrametadata)) +
  geom_point()

# Test my QGIS bioclim conferral vs. Hannah's
extra$catalogID <- as.character(extra$catalogID)
check_TP <- tp %>%
  left_join(., extra, by = c("sampleID" = "catalogID"))
plot(check_TP$MAP_bc, check_TP$MAP)
plot(check_TP$MAT_bc, check_TP$MAT)

# Map
coords <- uniq %>%
  select(sampleID, "latitude..decimal.degrees.", "longitude..decimal.degrees.") %>%
  set_names(c("sampleID", "latitude", "longitude"))
#write.csv(coords, "coords49.csv")
coords_trans <- st_as_sf(coords, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
ozmap()
sf_oz <- ozmap("states")
ggplot(data = sf_oz) + 
  geom_sf() +
  geom_sf(data = coords_trans) +
  theme_minimal()

# Hannah map
#base_map <- readRDS("data/base_map.RDS")
#base_map

# PhyloFlash output
pf <- read.delim("data/Australia_NTU_table.tsv") %>%
  rename(taxonomy = Taxonomy)
#sum(uniqL12$sampleID %in% names(pf))
#out_fp <- "seqtab_wTax_mctoolsr.txt"
#write("#Exported for mctoolsr", out_fp)
#suppressWarnings(write.table(pf, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Import mctoolsr
tax_table_fp <- "data/seqtab_wTax_mctoolsr.txt"
map_fp <- "data/metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 75 samples loaded

# Filter Chloroplast, Mitochondria, Domain NA
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 147 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 165 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 0 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 1681 removed

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 5542.8
se(colSums(input_filt$data_loaded)) # 492.6347
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 1015) # n = 75 still
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, index = "shannon", MARGIN = 2)

# Save
#saveRDS(input_filt_rare, "data/input_filt_rare_AusPF.rds")



#### _Combine Samples ####
# Combine samples with the same name that have L1, L2
# PhyloFlash output
pf_tax <- read.delim("data/Australia_NTU_table.tsv") %>%
  rename(taxonomy = Taxonomy)
pf <- read.delim("data/Australia_NTU_table.tsv") %>%
  rename(taxonomy = Taxonomy) %>%
  select(-taxonomy) %>%
  column_to_rownames(var = "NTU") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  separate(SampleID, into = c("SampleID", "Other"), sep = "_PE_") %>%
  select(-Other) %>%
  mutate(SampleID = substr(SampleID, start = 1, stop = nchar(SampleID) - 2)) %>%
  group_by(SampleID) %>%
  summarise_each(list(sum)) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  as.data.frame() %>%
  mutate(taxonomy = pf_tax$taxonomy) %>%
  rownames_to_column(var = "NTU")
# out_fp <- "data/seqtab_wTax_mctoolsr_comb.txt"
# write("#Exported for mctoolsr", out_fp)
# suppressWarnings(write.table(pf, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

metadata <- metadata %>%
  left_join(., extra, by = "SampleID") %>%
  select(SampleID, everything())
#write.table(metadata, "data/metadata_hannah.txt", sep = "\t", row.names = F)

tax_table_fp <- "data/seqtab_wTax_mctoolsr_comb.txt"
map_fp <- "data/metadata_hannah.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 49 samples loaded

# Filter Chloroplast, Mitochondria, Domain NA
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 147 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 165 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 0 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 1681 removed

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 8483.878
se(colSums(input_filt$data_loaded)) # 997.5413
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 2046) # n = 49
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, index = "shannon", MARGIN = 2)


# Save
#saveRDS(input_filt_rare, "data/input_filt_rare_AusPF_comb.rds")



#### _Start Here ####
input_filt_rare <- readRDS("data/input_filt_rare_AusPF_comb.rds")



#### 2. Focal Taxa ####
# Plot the 8 focal genera across the MAP gradient
# Bryobacter, Solibacter, Acidothermus, Mycobacterium, Bacillus, Haliangium, Bradyrhizobium, Udaeobacter
# Also get data from Hannah which included aridity index
View(input_filt_rare$taxonomy_loaded)
tax_sum_gen <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 6, 
                                  report_higher_tax = F)
plot_taxa_bars(tax_sum_gen, 
               input_filt_rare$map_loaded, 
               type_header = "vegetation_type", 
               num_taxa = 12)

tax_sum_gen_focal <- as.data.frame(t(tax_sum_gen)) %>%
  select(Bryobacter, `Candidatus Solibacter`, Acidothermus, Mycobacterium, Bacillus, 
         Haliangium, Bradyrhizobium, `Candidatus Udaeobacter`) %>%
  mutate(sampleID = rownames(.))

d_focal <- input_filt_rare$map_loaded %>%
  mutate(sampleID = rownames(.),
         catalogID = as.character(catalogID)) %>%
  left_join(., tax_sum_gen_focal, by = "sampleID") %>%
  left_join(., tp, by = c("catalogID" = "sampleID")) %>%
  arrange(MAP)
head(d_focal$sampleID) # Driest
tail(d_focal$sampleID) # Wettest

ggplot(d_focal, aes(MAP, Bryobacter)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, `Candidatus Solibacter`)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, Acidothermus)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, Mycobacterium)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, Bacillus)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, Haliangium)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, Bradyrhizobium)) +
  geom_point() +
  geom_smooth()
ggplot(d_focal, aes(MAP, `Candidatus Udaeobacter`)) +
  geom_point() +
  geom_smooth()

d_focal_long <- d_focal %>%
  select(-name) %>%
  pivot_longer(cols = c(Bryobacter, `Candidatus Solibacter`, Acidothermus, Mycobacterium, Bacillus, 
                        Haliangium, Bradyrhizobium, `Candidatus Udaeobacter`))
ggplot(d_focal_long, aes(MAP, value)) +
  geom_point() +
  labs(x = "MAP (cm)",
       y = "Relative abundance") +
  facet_wrap(~ name, ncol = 4) +
  theme_bw()
ggplotly(ggplot(d_focal_long, aes(MAP, value)) +
           geom_point() +
           labs(x = "MAP (cm)",
                y = "Relative abundance") +
           facet_wrap(~ name, ncol = 4) +
           theme_bw())

m <- lm(Bryobacter ~ MAP, data = d_focal)
summary(m) # NS
m <- lm(`Candidatus Solibacter` ~ MAP, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(Acidothermus ~ MAP, data = d_focal)
summary(m) # NS
m <- lm(Mycobacterium ~ MAP, data = d_focal)
summary(m) # NS
m <- lm(Bacillus ~ MAP, data = d_focal)
summary(m) # NS
m <- lm(Haliangium ~ MAP, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(Bradyrhizobium ~ MAP, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(`Candidatus Udaeobacter` ~ MAP, data = d_focal)
summary(m) # NS

m <- lm(Bryobacter ~ AridityIndex, data = d_focal)
summary(m) # NS
m <- lm(`Candidatus Solibacter` ~ AridityIndex, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(Acidothermus ~ AridityIndex, data = d_focal)
summary(m) # Marg. Sig. Pos.
m <- lm(Mycobacterium ~ AridityIndex, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(Bacillus ~ AridityIndex, data = d_focal)
summary(m) # NS
m <- lm(Haliangium ~ AridityIndex, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(Bradyrhizobium ~ AridityIndex, data = d_focal)
summary(m) # Sig. Pos.
m <- lm(`Candidatus Udaeobacter` ~ AridityIndex, data = d_focal)
summary(m) # Sig. Pos.
sig_tax <- c("Candidatus Solibacter", "Mycobacterium", 
             "Haliangium", "Bradyrhizobium", "Candidatus Udaeobacter")

plot(extra$MAP, extra$AridityIndex) # Positive relationship - not like some other aridity indices!
ggplot(d_focal, aes(AridityIndex, MoistureIndex)) +
  geom_point()

pdf("InitialFigs/AridityRelAbund.pdf", width = 7, height = 5)
ggplot(d_focal_long, aes(AridityIndex, value)) +
  geom_point() +
  geom_smooth(data = subset(d_focal_long, name %in% sig_tax),
              method = "lm") +
  labs(x = "Aridity Index",
       y = "Relative abundance") +
  facet_wrap(~ name, ncol = 4) +
  theme_bw()
dev.off()
