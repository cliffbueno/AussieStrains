# Compare BASE and NEON
# Assess sequencing depth and aridity gradients

# Libraries
library(tidyverse)
library(scales)
library(plotly)
library(readxl)
library(usmap)
library(ggrepel)
library(ggspatial)
library(ggstance)

# WD
setwd("~/Desktop/Strains/")


#### Depth ####
#d <- read.csv("CheckSeqDepth.csv")
#t.test(d$Reads ~ d$Dataset)

# pdf("BASEvsNEON_Reads.pdf", width = 7, height = 5)
# ggplot(d, aes(x = Dataset, y = Reads)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(size = 3, width = 0.25, alpha = 0.75) +
#   labs(x = "Dataset", y = "# Reads") +
#   scale_y_continuous(labels = label_comma()) +
#   ggtitle(label = NULL, subtitle = "BASE: 50 sites, 76 metagenomes (1 per site, some sequenced twice)\nNEON: 2019 data, 24 sites, 72 metagenomes (3 randomly selected from each)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12))
# dev.off()

# Add 2017 NEON data
d <- read.csv("CheckSeqDepth.csv")
m <- aov(d$Reads ~ d$Dataset)
summary(m)
TukeyHSD(m)
high <- d %>%
  filter(Reads > 10000000) %>%
  group_by(Dataset) %>%
  summarise(num = n())
text_df <- data.frame(Dataset = c("BASE", "NEON_2017", "NEON_2019"),
                      y = c(50000000, 50000000, 50000000),
                      label = c("33 over 10 mill", "35 over 10 mill", "1 over 10 mill"))

pdf("BASEvsNEONyears_Reads.pdf", width = 7, height = 5)
ggplot(d, aes(x = Dataset, y = Reads)) +
  geom_hline(yintercept = 10000000, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.25, alpha = 0.75) +
  geom_text(data = text_df, aes(Dataset, y, label = label)) +
  labs(x = "Dataset", y = "# Reads") +
  scale_y_continuous(labels = label_comma()) +
  ggtitle(label = NULL, subtitle = "BASE: 50 sites shown, 76 metagenomes (1 per site, some sequenced twice)\nNEON 2017: max 36 sites, 108 metagenomes (3 randomly selected from each)\nNEON 2019: max 24 sites, 72 metagenomes (3 randomly selected from each)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()



#### MAP ####
neon_sites <- read_xlsx("NEON/NEON_Soil_MG.xlsx") %>%
  filter(`2017` == 1)
neon_meta <- read.csv("NEON/NEON_Field_Site_Metadata_20240827.csv") %>%
  filter(field_site_id %in% neon_sites$Site)

# Map of 36 2017 MG sites
coords <- neon_meta %>%
  dplyr::select(field_site_id, field_longitude, field_latitude) %>%
  rename(Site = field_site_id,
         Latitude = field_latitude,
         Longitude = field_longitude)
#write.csv(coords, "NEON/neon_2017_site_coords.csv")
test_data <- data.frame(lon = coords$Longitude, lat = coords$Latitude)
transformed_data <- usmap_transform(test_data)
coords <- coords %>%
  left_join(., transformed_data, by = c("Latitude" = "lat"))
pdf("NEON/Map.pdf", width = 8, height = 6)
map <- plot_usmap(exclude = "HI",
                  color = "white",
                  fill = "grey80",
                  size = 0.3) +
  geom_point(data = coords, 
             aes(x = x, y = y),
             fill = "red",
             color = "black",
             size = 4,
             shape = 21) +
  geom_text_repel(data = coords,
                  aes(x = x, y = y, label = Site),
                  size = 4) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr") +
  theme(legend.position = "none")
map
dev.off()
ggplot(neon_meta, aes(field_mean_annual_precipitation_mm)) +
  geom_histogram() +
  labs(x = "NEON Site MAP (mm)",
       y = "Count") +
  theme_bw()

min(neon_meta$field_mean_annual_precipitation_mm)
max(neon_meta$field_mean_annual_precipitation_mm)
mean(neon_meta$field_mean_annual_precipitation_mm)
sd(neon_meta$field_mean_annual_precipitation_mm)

moist <- read.csv("Australia/MoistureIndexALA.csv") %>%
  select(decimalLatitude, 
         #decimalLongitude, 
         Moisture.Index...annual.mean..Bio28.) %>%
  rename(latitude = decimalLatitude,
         #longitude = decimalLongitude,
         MoistureIndex = Moisture.Index...annual.mean..Bio28.) %>%
  mutate(latjoin = round(latitude, digits = 4)) %>%
  select(-latitude)
tp <- read.csv("Australia/coords_wTP.csv") %>%
  mutate(sampleID = as.character(sampleID)) %>%
  mutate(latjoin = round(latitude, digits = 4)) %>%
  left_join(., moist, by = "latjoin") %>%
  rename(MAP_bc = MAP,
         MAT_bc = MAT) %>%
  select(sampleID, MAP_bc, MAT_bc, MoistureIndex)

min(tp$MAP_bc)
max(tp$MAP_bc)
mean(tp$MAP_bc)
sd(tp$MAP_bc)

ggplot(tp, aes(MAP_bc)) +
  geom_histogram() +
  labs(x = "BASE Site MAP (mm)",
       y = "Count") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900,
                                1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
                                1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
                                2600, 2700)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine
neon_map <- neon_meta %>%
  select(field_site_id, field_mean_annual_precipitation_mm) %>%
  rename(ID = field_site_id,
         MAP = `field_mean_annual_precipitation_mm`) %>%
  mutate(Dataset = "NEON_2017")
base_map <- tp %>%
  select(sampleID, MAP_bc) %>%
  rename(ID = sampleID,
         MAP = MAP_bc) %>%
  mutate(Dataset = "BASE")
comb <- rbind(neon_map, base_map)

pdf("BASEvsNEON_Precip.pdf", width = 7, height = 5)
ggplot(comb, aes(x = Dataset, y = MAP)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.25, alpha = 0.75) +
  #geom_text(data = text_df, aes(Dataset, y, label = label)) +
  labs(x = "Dataset", y = "MAP (mm)") +
  ggtitle(label = NULL, subtitle = "BASE: n = 49 shown, min MAP = 193, max MAP = 2634\nNEON 2017: max n = 36, min MAP = 105, max MAP = 2451") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()


#### Aridity ####
# Loaded 2.5 min BIOCLIM 12 and eto into QGIS
# Point sampling tool
# Aridity = MAP/PET
neon_aridity <- read.csv("NEON/neon_2017_site_coords_eto_map.csv") %>%
  mutate(aridity = wc2.1_2.5m/et0_yr) %>%
  mutate(Dataset = "NEON_2017 (n = 36)")
base_aridity <- read.csv("Australia/base_49_site_coords_eto_map.csv") %>%
  mutate(aridity = wc2.1_2.5m/et0_yr) %>%
  mutate(Dataset = "BASE_Hannah (n = 49)") %>%
  rename(Site = sampleID,
         Latitude = latitude,
         Longitude = longitude)
a <- rbind(neon_aridity, base_aridity)
ggplot(a, aes(x = Dataset, y = aridity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.75) +
  #geom_density(aes(y = Dataset), inherit.aes = FALSE)
  #geom_text(data = text_df, aes(Dataset, y, label = label)) +
  labs(x = "Dataset", y = "Aridity index (MAP/PET)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

pdf("BASEvsNEON_Aridity.pdf", width = 7, height = 5)
ggplot(a, aes(x = aridity, y = -0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxploth(outlier.shape = NA) +
  geom_jitter(size = 3, height = 0.1, alpha = 0.75) +
  geom_density(aes(x = aridity), inherit.aes = FALSE) +
  scale_y_continuous(labels = c("", "0", "0.5", "1"),
                     breaks = c(-0.5, 0, 0.5, 1)) +
  labs(x = "Aridity index (MAP/PET)",
       y = "Density") +
  facet_grid(Dataset ~ .) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.ticks.y = element_blank())
dev.off()

ggplot(a, aes(x = aridity, y = -0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxploth(outlier.shape = NA, aes(fill = Dataset)) +
  geom_jitter(size = 3, height = 0.1, alpha = 0.75, aes(colour = Dataset, group = Dataset)) +
  geom_density(aes(x = aridity, fill = Dataset), inherit.aes = FALSE) +
  scale_y_continuous(labels = c("", "0", "0.5", "1"),
                     breaks = c(-0.5, 0, 0.5, 1)) +
  labs(x = "Aridity index (MAP/PET)",
       y = "Density") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.ticks.y = element_blank()) +
  coord_flip()

# test
ggplot(diamonds, aes(x = carat, y = -0.5)) +
  
  # horizontal box plot
  geom_boxploth(aes(fill = cut)) +
  
  # normal density plot
  geom_density(aes(x = carat), inherit.aes = FALSE) +
  
  # vertical lines at Q1 / Q2 / Q3
  stat_boxploth(geom = "vline", aes(xintercept = ..xlower..)) +
  stat_boxploth(geom = "vline", aes(xintercept = ..xmiddle..)) +
  stat_boxploth(geom = "vline", aes(xintercept = ..xupper..)) +
  
  facet_grid(cut ~ .) +
  
  # reproduce original chart's color scale (o/w ordered factors will result
  # in viridis scale by default, using the current version of ggplot2)
  scale_fill_discrete()
  
  