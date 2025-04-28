# Analysis for soil bacterial strain manuscript
# By Cliff Bueno de Mesquita, Fierer Lab, Fall 2024, Spring 2025

#### 1. Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(ggside)
library(ozmaps)
library(sf)
library(ggspatial)
library(scales)
library(arules)
library(vegan)
library(corrplot)
library(ggrepel)
library(mctoolsr)
library(FSA)
library(RColorBrewer)
library(broom)
library(readxl)
library(pheatmap)
library(writexl)
library(geosphere)
library(ComplexUpset)
library(UpSetR)
library(rcompanion)
library(ape)
library(ggtree)
library(picante)
library(ggnewscale)
library(bestglm)
library(phyloseq)
library(phytools)
library(gdm)
library(treeio)
library(cowplot)
library(sp)

# Functions
`%notin%` <- Negate(`%in%`)
source("~/Documents/GitHub/SunflowerGxE/code/cliffplot_taxa_bars.R")
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
move_to_lower_triangle <- function(mat) {
  n <- nrow(mat)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # If there's a non-NA in the upper triangle and space in the lower triangle
      if (!is.na(mat[i, j]) && is.na(mat[j, i])) {
        mat[j, i] <- mat[i, j]  # Move to the lower triangle
        mat[i, j] <- NA         # Set the original position to NA
      }
    }
  }
  return(mat)
}
replace_lower_triangle_na <- function(mat) {
  n <- nrow(mat)
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      if (is.na(mat[i, j])) {
        mat[i, j] <- 0
      }
    }
  }
  return(mat)
}

# Working directory
setwd("~/Documents/GitHub/AussieStrains/")

d <- read.csv("data/metadata104.csv") %>%
  dplyr::select(-X, -Prefix) %>%
  mutate(sampleID = as.character(sampleID))

# Output just sampleID and coordinates
c <- d %>%
  dplyr::select(sampleID, latitude, longitude)
#write.csv(c, "data/coords104.csv")

# Elevation, slope, aspect were missing, so conferred for all using DEM
#dem <- read.csv()

meta <- d %>%
  #left_join(., dem, by = "sampleID")
  mutate(Sample = paste("X", sampleID, sep = "")) %>%
  dplyr::select(Sample, everything())

# Save metadata as text for mctoolsr import
length(unique(meta$Sample))
#write.table(meta, "data/metadata104wDEM.txt", row.names = F, sep = "\t")



#### 2. Environment ####
# Assess environmental variation in the 104 samples
names(d)
table(d$vegetation_type)
table(d$local_class)
env_vars <- c("boron_hot_cacl2", "clay", "conductivity",
              "dtpa_copper", "dtpa_iron", "dtpa_manganese", "dtpa_zinc",
              #"elev",
              "exc_aluminium", "exc_calcium", "exc_magnesium",
              "exc_potassium", "exc_sodium", "latitude",
              "longitude", "nitrate_nitrogen", "organic_carbon",
              "ph", "phosphorus_colwell",
              "sand", "silt", "sulphur", 
              "water_content", # 51 NA
              "bio1", "bio12", "AI")
d_AI <- d %>%
  dplyr::select(sampleID, AI) %>%
  mutate(sampleID = as.character(sampleID))
d_ai_sort <- d %>%
  arrange(AI)
d_env <- d %>%
  dplyr::select(all_of(env_vars))
hist(d_env$organic_carbon)
hist(d_env$nitrate_nitrogen)
hist(d_env$ph)
hist(d_env$bio1)
hist(d_env$bio12)
hist(d_env$clay)
#hist(d_env$water_content)
hist(d_env$conductivity)
#hist(d_env$elev)
hist(d_env$phosphorus_colwell)
hist(d_env$exc_potassium)
# plot(d_env$potassium_colwell, d_env$exc_potassium) # just use one

# Check NA by variable
n_na <- c()
for (i in 1:ncol(d_env)) {
  n_na[i] <- sum(is.na(d_env[,i]))
}
n_na
# Most with 27, 37 NA
# 8 with no NAs

# 40 samples with data for all 19 of these columns! no NA.
d_env_nona <- d_env %>%
  filter_at(vars(1:24), all_vars(!is.na(.))) %>%
  set_names(c("B", "Clay", "Conductivity", "Cu", "Fe", "Mn", "Zn", "Al", "Ca", "Mg",
              "K", "Na", "Lat", "Long", "NO3", "Org. C", "pH", "Colwell P",
              "Sand", "Silt", "S", "H2O", "Temp", "Precip", "Aridity")) %>%
  dplyr::select(Clay, Silt, Sand, `Org. C`, NO3, `Colwell P`, pH, Conductivity,
                B, Cu, Fe, Mn, Zn, Al, Ca, Mg, K, Na, S,
                Lat, Long, H2O, Temp, Precip, Aridity)

# Corrplot
m <- cor(d_env_nona)
pdf("InitialFigs/Env_Corrplot.pdf", width = 7, height = 5)
corrplot(m, 
         method = "square",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5)
dev.off()

# PCA with vectors
d.pcx <- prcomp(d_env_nona)
set.seed(100)
ef <- envfit(d.pcx, d_env_nona, permutations = 999, na.rm = TRUE)
ef
ordiplot(d.pcx)
plot(ef, p.max = 0.05, cex = 0.5)
manual_factor <- 0.2
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(PC1 = PC1 * manual_factor,
         PC2 = PC2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = variables)
d.mvar <- sum(d.pcx$sdev^2)
PC1 <- paste("PC1: ", round((sum(d.pcx$sdev[1]^2)/d.mvar)*100, 1), "%")
PC2 <- paste("PC2: ", round((sum(d.pcx$sdev[2]^2)/d.mvar)*100, 1), "%")
d_env_nona$Axis01 <- vegan::scores(d.pcx)[,1]
d_env_nona$Axis02 <- vegan::scores(d.pcx)[,2]
d_env_nona$Axis01 <- d_env_nona$Axis01/sqrt(sum((d_env_nona$Axis01 - mean(d_env_nona$Axis01))^2))
d_env_nona$Axis02 <- d_env_nona$Axis02/sqrt(sum((d_env_nona$Axis02 - mean(d_env_nona$Axis02))^2))
pdf("InitialFigs/Env_PCA.pdf", width = 7, height = 5)
ggplot(d_env_nona, aes(Axis01, Axis02)) +
  geom_point(size = 3, alpha = 1, pch = 16) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "red") +
  labs(x = PC1, 
       y = PC2) +
  theme_bw()
dev.off()



#### 3. Taxa ####
# Use mTAGs to select abundant and ubiquitous soil bacteria

#### _Setup ####
# mTAGs output. Format for mctoolsr
# Need to separate prokaryotes and eukaryotes
mt <- read.delim("data/merged_profile.otu.tsv") %>%
  rename(taxonomy = X.taxpath) %>%
  filter(., !grepl("Eukaryota", taxonomy)) %>%
  mutate(taxonomy = gsub("root__Root;", "", taxonomy)) %>%
  mutate(taxonomy = gsub("domain__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("phylum__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("class__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("order__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("family__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("genus__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("otu__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("silva_138_complink_cons_", "", taxonomy)) %>%
  mutate(taxonomy = gsub("unknown", "NA", taxonomy)) %>%
  mutate(taxonomy = gsub("otu_", "NA;otu_", taxonomy)) %>%
  filter(taxonomy != "Unassigned") %>%
  filter(taxonomy != "Unaligned") %>%
  separate(taxonomy, into = c("a","s","d","f","g","h","j","otu"), remove = F, sep = ";") %>%
  dplyr::select(-c(a,s,d,f,g,h,j)) %>%
  dplyr::select(otu, everything()) %>%
  dplyr::select(-taxonomy, taxonomy)
n <- data.frame(name = names(mt)) %>%
  separate(name, into = c("Sample", "Junk"), sep = "_")
mt <- mt %>%
  set_names(n$Sample)
# out_fp <- "data/seqtab_wTax_mctoolsr_mtagProk.txt"
# names(mt)[1] = "#OTU_ID"
# write("#Exported for mctoolsr", out_fp)
# suppressWarnings(write.table(mt, 
#                              out_fp, 
#                              sep = "\t", 
#                              row.names = FALSE, 
#                              append = TRUE))
sum(n$Sample %in% meta$Sample)

# Import mctoolsr
tax_table_fp <- "data/seqtab_wTax_mctoolsr_mtagProk.txt"
map_fp <- "data/metadata104wDEM.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 104 samples loaded

# Filter Chloroplast, Mitochondria, Domain NA
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 92 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 53 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 0 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 10291 removed

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 4824.74
se(colSums(input_filt$data_loaded)) # 149.4816
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 1849) # n = 104 still
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", MARGIN = 2)

# Save
#saveRDS(input_filt_rare, "data/input_filt_rare_mtags.rds")



#### _Analysis ####
input_filt_rare <- readRDS("data/input_filt_rare_mtags.rds")

# Check alpha diversity and aridity
m1 <- lm(rich ~ AI, data = input_filt_rare$map_loaded)
shapiro.test(m1$residuals)
summary(m1) # R2 = 0.11, p < 0.001
m2 <- lm(shannon ~ AI, data = input_filt_rare$map_loaded)
shapiro.test(m2$residuals)
summary(m2) # R2 = 0.11, p < 0.001
pdf("InitialFigs/Alpha_Rich_AI.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, alpha = 0.75, pch = 16, aes(colour = vegetation_type)) +
  geom_smooth(method = "lm") +
  labs(x = "Aridity index",
       y = "Richness (# mTAG OTUs)",
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  theme_bw()
dev.off()
ggplot(input_filt_rare$map_loaded, aes(AI, shannon)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  geom_smooth(method = "lm") +
  labs(x = "Aridity index",
       y = "Shannon diversity") +
  theme_bw()

# Check beta diversity
bc <- calc_dm(input_filt_rare$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
d_env <- input_filt_rare$map_loaded %>%
  dplyr::select(all_of(env_vars)) %>%
  set_names(c("B", "Clay", "Conductivity", "Cu", "Fe", "Mn", "Zn", "Al", "Ca", "Mg",
              "K", "Na", "Lat", "Long", "NO3", "Org. C", "pH", "Colwell P",
              "Sand", "Silt", "S", "Temp", "Precip", "Aridity")) %>%
  dplyr::select(Clay, Silt, Sand, `Org. C`, NO3, `Colwell P`, pH, Conductivity,
                B, Cu, Fe, Mn, Zn, Al, Ca, Mg, K, Na, S,
                Lat, Long, Temp, Precip, Aridity)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = variables)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_filt_rare$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/Beta_Bray_PCoA.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = vegetation_type)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Check taxa
input_filt_rare$map_loaded$Sample = paste("X", input_filt_rare$map_loaded$sampleID, sep = "")
cliffplot_taxa_bars(input_filt_rare, 2, variable = "Sample")
cliffplot_taxa_bars(input_filt_rare, 3, variable = "Sample")
cliffplot_taxa_bars(input_filt_rare, 4, variable = "Sample")
cliffplot_taxa_bars(input_filt_rare, 5, variable = "Sample")
cliffplot_taxa_bars(input_filt_rare, 6, variable = "Sample")

# Identify ubiquitous and abundant genera
# Note: based on Hannah's 48 samples, we were thinking:
# Acidothermus, Bacillus, Bradyrhizobium, Bryobacter, 
# Haliangium, Mycobacterium, Solibacter, Udaeobacter
# New additions >95% prev and 0.3% abund are Conexibacter, Singulisphaera, Nocardioides 
View(input_filt_rare$taxonomy_loaded)
tax_sum_gen <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 6, 
                                  report_higher_tax = T,
                                  relative = T)
toBin <- tax_sum_gen %>%
  filter(grepl("Bradyrhizobium|Streptomyces|Udaeobacter|Mycobacterium|Acidothermus",
               rownames(.))) %>%
  t() %>%
  as.data.frame()
# Brady: 401644
# Strepto: 138530
# Udaeo: 12818
# Myco: 401610
# Acido: 401654
gen_abund <- data.frame("Genus" = rownames(tax_sum_gen),
                        "MeanPercRelAbund" = rowMeans(tax_sum_gen)*100) %>%
  arrange(desc(MeanPercRelAbund))
gen_prev <- data.frame("Genus" = rownames(tax_sum_gen),
                       "Absent" = rowSums(tax_sum_gen==0)) %>%
  mutate(Present_n = ncol(tax_sum_gen) - Absent) %>%
  mutate(Present_Perc = Present_n/ncol(tax_sum_gen)*100) %>%
  arrange(desc(Present_Perc))
# Note: There were 1067 genera including NA and uncultured. 690 with assigned names.
gen_prev_abund <- gen_abund %>%
  left_join(., gen_prev, by = "Genus") %>%
  separate(Genus, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
          sep = "; ") %>%
  filter(Genus != "NA") %>%
  filter(Genus != "uncultured") %>%
  dplyr::select(-Absent, -Present_n)
# Export this for Noah
#writexl::write_xlsx(gen_prev_abund, "data/gen_prev_abund.xlsx", format_headers = F)
gen_90_point1 <- gen_prev_abund %>%
  filter(MeanPercRelAbund > 0.1) %>%
  filter(Present_Perc > 90) # 7/8
gen_85_point5 <- gen_prev_abund %>%
  filter(MeanPercRelAbund > 0.5) %>%
  filter(Present_Perc > 85) %>%
  mutate(Genus2 = gsub("Candidatus", "", Genus)) %>%
  mutate(Genus3 = gsub("Candidatus", "Candidatus ", Genus))

# At this step need to be less conservative. Want to input more to Sylph, then use that to refine.
# Top 10% of the named genera is about 0.1% rel abund. And try 75% prevalence.
# 40 genera
gen_75_point1 <- gen_prev_abund %>%
  filter(MeanPercRelAbund > 0.1) %>%
  filter(Present_Perc > 75) %>%
  mutate(Genus2 = gsub("Candidatus", "", Genus)) %>%
  mutate(Genus3 = gsub("Candidatus", "Candidatus ", Genus))

# Check if there are dominant OTUs for each focal genus
tax_sum_otu <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 8, 
                                  report_higher_tax = T,
                                  relative = T)
otu_abund <- data.frame("OTU" = rownames(tax_sum_otu),
                        "MeanPercRelAbund" = rowSums(tax_sum_otu)) %>%
  arrange(desc(MeanPercRelAbund))
otu_prev <- data.frame("OTU" = rownames(tax_sum_otu),
                       "Absent" = rowSums(tax_sum_otu==0)) %>%
  mutate(Present_n = ncol(tax_sum_otu) - Absent) %>%
  mutate(Present_Perc = Present_n/ncol(tax_sum_otu)*100) %>%
  arrange(desc(Present_Perc))
otu_prev_abund <- otu_abund %>%
  left_join(., otu_prev, by = "OTU") %>%
  separate(OTU, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus",
                         "Species", "OTU"),
           sep = "; ") %>%
  filter(Genus != "NA") %>%
  filter(Genus != "uncultured") %>%
  filter(Genus %in% gen_75_point1$Genus) %>%
  dplyr::select(-Absent, -Present_n)
genus_otus <- otu_prev_abund %>%
  group_by(Genus) %>%
  summarise(n_OTUs = n()) %>%
  arrange(desc(n_OTUs))
sum(genus_otus$n_OTUs)
# 2781 OTUs across those 15 genera! Most Bryobacter.

# Plot ubiq/abund genera over aridity gradient
tax_sum_gen <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 6, 
                                  report_higher_tax = F,
                                  relative = T)
tax_sum_gen_focal <- as.data.frame(t(tax_sum_gen)) %>%
  select(gen_75_point1$Genus) %>%
  mutate(sampleID = rownames(.)) %>%
  rename("Solibacter" = "CandidatusSolibacter",
         "Udaeobacter" = "CandidatusUdaeobacter",
         "Xiphinematobacter" = "CandidatusXiphinematobacter",
         "Koribacter" = "CandidatusKoribacter")
d_focal <- input_filt_rare$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., tax_sum_gen_focal, by = "sampleID")
d_focal_long <- d_focal %>%
  pivot_longer(cols = gen_75_point1$Genus2)
# Test linear and quadratic models
results <- as.data.frame(matrix(NA, nrow(gen_75_point1), 6)) %>%
  set_names(c("Genus", "LinearR2", "LinearP", "QuadraticR2", "QuadraticP", "Compare"))
test_data <- d_focal %>%
  dplyr::select(gen_75_point1$Genus2, AI)
for (i in 1:nrow(gen_75_point1)) {
  results$Genus[i] <- gen_75_point1$Genus2[i]
  m1 <- lm(test_data[,i] ~ test_data$AI)
  results$LinearR2[i] <- summary(m1)$r.squared
  results$LinearP[i] <- summary(m1)$coefficients[2,4]
  m2 <- lm(test_data[,i] ~ test_data$AI + I(test_data$AI^2))
  results$QuadraticR2[i] <- summary(m2)$r.squared
  results$QuadraticP[i] <- glance(m2)$p.value
  results$Compare[i] <- anova(m1, m2)$'Pr(>F)'[2]
}
# Decide best model or N.S.
results <- results %>%
  mutate(Model = ifelse(LinearP > 0.05 & QuadraticP > 0.05, "NS",
                        ifelse(Compare < 0.05, "Quadratic", "Linear")))
quad <- results %>%
  filter(Model == "Quadratic")
lin <- results %>%
  filter(Model == "Linear")
ggplot(d_focal_long, aes(AI, value)) +
  geom_point() +
  geom_smooth(data = subset(d_focal_long, name %in% lin$Genus),
              method = "lm") +
  geom_smooth(data = subset(d_focal_long, name %in% quad$Genus),
              method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Aridity index",
       y = "Relative abundance") +
  facet_wrap(~ name, ncol = 8, scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 6))

# Subset to the 33 actually used (need to make genus_genomes_2 first)
genus_info <- read_xlsx("data/genus_info_75point1.xlsx")
genus_genomes_2 <- genus_info %>%
  filter(n_GTDB_genomes >= 2)
d_focal_long_used <- d_focal_long %>%
  filter(name %in% genus_genomes_2$Genus) %>%
  droplevels()
pdf("InitialFigs/TopGenera_AI_75point1_used33.pdf", width = 12, height = 8)
ggplot(d_focal_long_used, aes(AI, value)) +
  geom_point() +
  geom_smooth(data = subset(d_focal_long_used, name %in% lin$Genus),
              method = "lm") +
  geom_smooth(data = subset(d_focal_long_used, name %in% quad$Genus),
              method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Aridity index",
       y = "Relative abundance") +
  facet_wrap(~ name, ncol = 7, scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 6))
dev.off()

# Looping back here after running Sylph
# Want to check OTU distributions for Bradyrhizobium, Streptomyces, Udaeobacter, Mycobacterium
mtag_brady <- filter_taxa_from_input(input_filt_rare,
                                     taxa_to_keep = "Bradyrhizobium",
                                     at_spec_level = 6)
mtag_brady$map_loaded$sampleID <- as.character(mtag_brady$map_loaded$sampleID)
mtag_brady_otu <- summarize_taxonomy(input = mtag_brady, 
                                     level = 8, 
                                     report_higher_tax = F,
                                     relative = T) %>%
  replace(is.na(.), 0)
mtag_brady_taxaSort <- data.frame(meanAbund = rowMeans(mtag_brady_otu),
                                  otu = rownames(mtag_brady_otu)) %>%
  arrange(meanAbund)
mtag_brady_aiSort <- mtag_brady$map_loaded %>%
  arrange(desc(AI))
g <- plot_ts_heatmap(mtag_brady_otu, 
                metadata_map = mtag_brady$map_loaded,
                type_header = "sampleID",
                min_rel_abund = 0,
                custom_sample_order = mtag_brady_aiSort$sampleID,
                custom_taxa_order = mtag_brady_taxaSort$otu,
                remove_other = T) +
  coord_flip() +
  ggtitle("Bradyrhizobium") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
g$layers[[2]] <- NULL
pdf("InitialFigs/otus_brady.pdf", width = 8, height = 8)
g
dev.off()

mtag_strepto <- filter_taxa_from_input(input_filt_rare,
                                     taxa_to_keep = "Streptomyces",
                                     at_spec_level = 6)
mtag_strepto$map_loaded$sampleID <- as.character(mtag_strepto$map_loaded$sampleID)
mtag_strepto_otu <- summarize_taxonomy(input = mtag_strepto, 
                                     level = 8, 
                                     report_higher_tax = F,
                                     relative = T) %>%
  replace(is.na(.), 0)
mtag_strepto_taxaSort <- data.frame(meanAbund = rowMeans(mtag_strepto_otu),
                                  otu = rownames(mtag_strepto_otu)) %>%
  arrange(meanAbund)
mtag_strepto_aiSort <- mtag_strepto$map_loaded %>%
  arrange(desc(AI))
g <- plot_ts_heatmap(mtag_strepto_otu, 
                     metadata_map = mtag_strepto$map_loaded,
                     type_header = "sampleID",
                     min_rel_abund = 0,
                     custom_sample_order = mtag_strepto_aiSort$sampleID,
                     custom_taxa_order = mtag_strepto_taxaSort$otu,
                     remove_other = T) +
  coord_flip() +
  ggtitle("Streptomyces") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
g$layers[[2]] <- NULL
pdf("InitialFigs/otus_strepto.pdf", width = 8, height = 8)
g
dev.off()

mtag_udaeo <- filter_taxa_from_input(input_filt_rare,
                                       taxa_to_keep = "Udaeobacter",
                                       at_spec_level = 6)
mtag_udaeo$map_loaded$sampleID <- as.character(mtag_udaeo$map_loaded$sampleID)
mtag_udaeo_otu <- summarize_taxonomy(input = mtag_udaeo, 
                                       level = 8, 
                                       report_higher_tax = F,
                                       relative = T) %>%
  replace(is.na(.), 0)
mtag_udaeo_taxaSort <- data.frame(meanAbund = rowMeans(mtag_udaeo_otu),
                                    otu = rownames(mtag_udaeo_otu)) %>%
  arrange(meanAbund)
mtag_udaeo_aiSort <- mtag_udaeo$map_loaded %>%
  arrange(desc(AI))
g <- plot_ts_heatmap(mtag_udaeo_otu, 
                     metadata_map = mtag_udaeo$map_loaded,
                     type_header = "sampleID",
                     min_rel_abund = 0,
                     custom_sample_order = mtag_udaeo_aiSort$sampleID,
                     custom_taxa_order = mtag_udaeo_taxaSort$otu,
                     remove_other = T) +
  coord_flip() +
  ggtitle("Udaeobacter") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
g$layers[[2]] <- NULL
pdf("InitialFigs/otus_udaeo.pdf", width = 8, height = 8)
g
dev.off()

mtag_myco <- filter_taxa_from_input(input_filt_rare,
                                     taxa_to_keep = "Mycobacterium",
                                     at_spec_level = 6)
mtag_myco$map_loaded$sampleID <- as.character(mtag_myco$map_loaded$sampleID)
mtag_myco_otu <- summarize_taxonomy(input = mtag_myco, 
                                     level = 8, 
                                     report_higher_tax = F,
                                     relative = T) %>%
  replace(is.na(.), 0)
mtag_myco_taxaSort <- data.frame(meanAbund = rowMeans(mtag_myco_otu),
                                  otu = rownames(mtag_myco_otu)) %>%
  arrange(meanAbund)
mtag_myco_aiSort <- mtag_myco$map_loaded %>%
  arrange(desc(AI))
g <- plot_ts_heatmap(mtag_myco_otu, 
                     metadata_map = mtag_myco$map_loaded,
                     type_header = "sampleID",
                     min_rel_abund = 0,
                     custom_sample_order = mtag_myco_aiSort$sampleID,
                     custom_taxa_order = mtag_myco_taxaSort$otu,
                     remove_other = T) +
  coord_flip() +
  ggtitle("Mycobacterium") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
g$layers[[2]] <- NULL
pdf("InitialFigs/otus_myco.pdf", width = 8, height = 8)
g
dev.off()

# Need to get the actual 16S sequence for the top rel abund OTU for those 4
# Download an mTAGs fasta for a sample with all of them
# Then blast Sylph genomes to see if the Sylph genome is same as the mTAGs OTU
top4_otu <- tax_sum_otu %>%
  filter(grepl("otu_24062|otu_187776|otu_97930|otu_41087", rownames(.)))
# 401570, 401586, 401590

#### _Strepto ####
tax_sum_fam <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 5, 
                                  report_higher_tax = F,
                                  relative = T)
fam_abund <- data.frame("Family" = rownames(tax_sum_fam),
                        "MeanPercRelAbund" = rowMeans(tax_sum_fam)*100) %>%
  arrange(desc(MeanPercRelAbund))
# Streptosporangiaceae is the 68th most abundant family across the dataset
# Mean % rel abund is 0.08788534 (cutoff was 0.1%)
fam_prev <- data.frame("Family" = rownames(tax_sum_fam),
                       "Absent" = rowSums(tax_sum_fam==0)) %>%
  mutate(Present_n = ncol(tax_sum_fam) - Absent) %>%
  mutate(Present_Perc = Present_n/ncol(tax_sum_fam)*100) %>%
  arrange(desc(Present_Perc))
# Streptosporangiaceae is the 37th most prevalent family. Present in 79/104 or 76% (cutoff was 75%)
# Make OTU heatmap as for other genera of interest

# Check Streptosporangiaceae OTUs
tax_sum_otu <- summarize_taxonomy(input = input_filt_rare, 
                                  level = 8, 
                                  report_higher_tax = T,
                                  relative = T)
otu_abund <- data.frame("OTU" = rownames(tax_sum_otu),
                        "MeanPercRelAbund" = rowSums(tax_sum_otu)) %>%
  arrange(desc(MeanPercRelAbund))
otu_prev <- data.frame("OTU" = rownames(tax_sum_otu),
                       "Absent" = rowSums(tax_sum_otu==0)) %>%
  mutate(Present_n = ncol(tax_sum_otu) - Absent) %>%
  mutate(Present_Perc = Present_n/ncol(tax_sum_otu)*100) %>%
  arrange(desc(Present_Perc))
otu_prev_abund <- otu_abund %>%
  left_join(., otu_prev, by = "OTU") %>%
  separate(OTU, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus",
                         "Species", "OTU"),
           sep = "; ") %>%
  filter(Family == "Streptosporangiaceae")
family_otus <- otu_prev_abund %>%
  group_by(Family) %>%
  summarise(n_OTUs = n()) %>%
  arrange(desc(n_OTUs))
sum(family_otus$n_OTUs) # 44 Streptosporangiaceae OTUs

mtag_sporang <- filter_taxa_from_input(input_filt_rare,
                                       taxa_to_keep = "Streptosporangiaceae",
                                       at_spec_level = 5)
mtag_sporang$map_loaded$sampleID <- as.character(mtag_sporang$map_loaded$sampleID)
mtag_sporang_otu <- summarize_taxonomy(input = mtag_sporang, 
                                       level = 8, 
                                       report_higher_tax = F,
                                       relative = T) %>%
  replace(is.na(.), 0)
mtag_sporang_taxaSort <- data.frame(meanAbund = rowMeans(mtag_sporang_otu),
                                    otu = rownames(mtag_sporang_otu)) %>%
  arrange(meanAbund)
mtag_sporang_aiSort <- mtag_sporang$map_loaded %>%
  arrange(desc(AI))
g <- plot_ts_heatmap(mtag_sporang_otu, 
                     metadata_map = mtag_sporang$map_loaded,
                     type_header = "sampleID",
                     min_rel_abund = 0,
                     custom_sample_order = mtag_sporang_aiSort$sampleID,
                     custom_taxa_order = mtag_sporang_taxaSort$otu,
                     remove_other = T) +
  coord_flip() +
  ggtitle("Streptosporangiaceae") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
g$layers[[2]] <- NULL
pdf("InitialFigs/otus_sporang.pdf", width = 8, height = 8)
g
dev.off()



#### _Brady ####
# Subset to the 53 samples used for StrainFinder
# Check which have the most Bradyrhizobium
# Get MAGs from those
# Goal: get a high quality MAG to compare to StrainFinder for validation
nrow(input_filt_rare$map_loaded)
nrow(d_brady)
class(input_filt_rare$map_loaded$sampleID)
input_filt_rare$map_loaded$sampleID_char <- as.character(input_filt_rare$map_loaded$sampleID)
class(d_brady$sampleID)
sum(d_brady$sampleID %in% input_filt_rare$map_loaded$sampleID_char)
sum(input_filt_rare$map_loaded$sampleID_char %in% d_brady$sampleID)
brady_mtags <- filter_data(input_filt_rare,
                           filter_cat = "sampleID_char",
                           keep_vals = d_brady$sampleID)
tax_sum_gen <- summarize_taxonomy(input = brady_mtags, 
                                  level = 6,
                                  report_higher_tax = F,
                                  relative = T)
brady_abund53 <- as.data.frame(t(tax_sum_gen)) %>%
  dplyr::select(Bradyrhizobium) %>%
  mutate(sampleID = rownames(.)) %>%
  arrange(desc(Bradyrhizobium))
  
  

#### 4. Ref Genomes ####
# Subset the GTDB full database to the genera of interest
bacGT <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") # takes a while to load
gtdb_focal <- bacGT %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";") %>%
  filter(Genus %in% gen_75_point1$Genus)
length(unique(gtdb_focal$Genus))
gen_75_point1$Genus
gen_75_point1$Genus %notin% gtdb_focal$Genus # Missing Cand. Udaeo, Cand. Soli, Acidibacter
# Knew Candidatus would cause issue. Check other taxonomy for Acidibacter. SILVA taxon.
# Tundrisphaera, Burkholderia-Caballeronia-Paraburkholderia, Actinomadura, P3OB-42
sum(grepl("Acidibacter", bacGT$gtdb_taxonomy)) # 0
sum(grepl("Acidibacter", bacGT$ncbi_taxonomy)) # 1
sum(grepl("Acidibacter", bacGT$ssu_gg_taxonomy)) # 0
sum(grepl("Acidibacter", bacGT$ssu_silva_taxonomy)) # 84
gtdb_acidibacter <- bacGT %>%
  filter(grepl("Acidibacter", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# Acidibacter is a mess. Several different genera in GTDB. 0 of 84 Acidibacter. Don't use.
gtdb_udaeobacter <- bacGT %>%
  filter(grepl("Udaeobacter", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# Udaeobacter isn't candidatus anymore. 34 of 82 Udaeobacter
gtdb_solibacter <- bacGT %>%
  filter(grepl("Solibacter", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# Solibacter isn't candidatus anymore. 3 of 51 Solibacter
gtdb_Tundrisphaera <- bacGT %>%
  filter(grepl("Tundrisphaera", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# Tundrisphaera has been reclassified to 2 genera
gtdb_bcp <- bacGT %>%
  filter(grepl("Burkholderia-Caballeronia-Paraburkholderia", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# These are split into 3 or more
gtdb_Actinomadura <- bacGT %>%
  filter(grepl("Actinomadura", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# These are split into Actinomadura, Spirillospora, Thermomonospora
gtdb_P3OB <- bacGT %>%
  filter(grepl("P3OB-42", ssu_silva_taxonomy)) %>%
  dplyr::select(gtdb_taxonomy, ncbi_taxonomy, ssu_gg_taxonomy, ssu_silva_taxonomy)
# Messy

# Actually should check how many were in SILVA for comparison
silva_focal <- bacGT %>%
  separate(ssu_silva_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus %in% gen_75_point1$Genus3)
genus_genomes_silva <- silva_focal %>%
  group_by(Genus) %>%
  summarise(n_Genomes_SILVA = n()) %>%
  arrange(desc(n_Genomes_SILVA))
  
# Subset to focal with Genus2
gtdb_focal <- bacGT %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus %in% gen_75_point1$Genus2)
length(unique(gtdb_focal$Genus)) 
# 35/40, missing Acidibacter, Tundrisphaera, Burkholderia-Caballeronia-Paraburkholderia, Actinomadura, P3OB-42

# First check how many genomes per focal genus
genus_genomes <- gtdb_focal %>%
  group_by(Genus) %>%
  summarise(n_Genomes = n()) %>%
  add_row(Genus = "Acidibacter", n_Genomes = 0) %>%
  add_row(Genus = "Tundrisphaera", n_Genomes = 0) %>%
  add_row(Genus = "Burkholderia-Caballeronia-Paraburkholderia", n_Genomes = 0) %>%
  add_row(Genus = "Actinomadura", n_Genomes = 0) %>%
  add_row(Genus = "P3OB-42", n_Genomes = 0) %>%
  arrange(desc(n_Genomes))

# Merge all of the info
genus_info <- gen_75_point1 %>%
  left_join(., genus_otus, by = "Genus") %>%
  left_join(., genus_genomes, by = c("Genus2" = "Genus")) %>%
  left_join(., genus_genomes_silva, by = c("Genus3" = "Genus")) %>%
  dplyr::select(-Genus, -Genus3) %>%
  rename(Genus = Genus2,
         n_GTDB_genomes = n_Genomes,
         n_GTDB_genomes_SILVAtax = n_Genomes_SILVA) %>%
  dplyr::select(Domain, Phylum, Class, Order, Family, Genus, 
                MeanPercRelAbund, Present_Perc, n_OTUs, n_GTDB_genomes, n_GTDB_genomes_SILVAtax)
# Export this for Noah
#writexl::write_xlsx(genus_info, "data/genus_info.xlsx", format_headers = F)
#writexl::write_xlsx(genus_info, "data/genus_info_75point1.xlsx", format_headers = F)

# Now, need to get a list of those NCBI accessions to download those genomes.
# For those 40 genera (35 with valid GTDB genus name), there are 21473 genomes.
# Need only 33 genera with valid GTDB genus name, and >= 2 GTDB genomes (n = 21471)
genus_info <- read_xlsx("data/genus_info_75point1.xlsx")
genus_genomes_2 <- genus_info %>%
  filter(n_GTDB_genomes >= 2)
gtdb_focal <- bacGT %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus %in% genus_genomes_2$Genus)
ncbi_accessions <- gtdb_focal %>%
  dplyr::select(ncbi_genbank_assembly_accession)
sum(is.na(ncbi_accessions$ncbi_genbank_assembly_accession)) # No NA.
sum(ncbi_accessions$ncbi_genbank_assembly_accession == "NA") # No NA string.
length(unique(ncbi_accessions$ncbi_genbank_assembly_accession)) # 21471
#write.table(ncbi_accessions, "data/ncbi_accessions.txt", sep = "\t", row.names = F, col.names = F)

# Met with Mads Albertson and apparently Acidothermus could actually Streptosporangiaceae
# So, download Streptosporangiaceae accessions and rerun Sylph
gtdb_strep <- bacGT %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Family == "Streptosporangiaceae")
ncbi_strep <- gtdb_strep %>%
  dplyr::select(ncbi_genbank_assembly_accession)
sum(is.na(ncbi_strep$ncbi_genbank_assembly_accession)) # No NA.
sum(ncbi_strep$ncbi_genbank_assembly_accession == "NA") # No NA string.
length(unique(ncbi_strep$ncbi_genbank_assembly_accession)) # 573
#write.table(ncbi_strep, "data/ncbi_strep.txt", sep = "\t", row.names = F, col.names = F)



#### _NCBI Download ####
# Take this list, and use the NCBI datasets command line interface to download genomes
# Looped through this:
# datasets download genome accession "${field1}" --filename "${field1}".zip --include genome
# Should have 21471 genomes
# First try downloaded 18670 genomes. What happenend to the other 2801!?
batch1 <- read.delim("data/batch1.txt", header = F) %>%
  mutate(accession = gsub(".zip", "", V1))
batch2 <- ncbi_accessions %>%
  filter(ncbi_genbank_assembly_accession %notin% batch1$accession)
#write.table(batch2, "data/ncbi_accessions_batch2.csv", sep = ",", row.names = F, col.names = F)
# Tried to download those with ncbi datasets and got this error:
# Error: invalid or unsupported assembly accession:
# But downloading them individually works!
# Try fixing the file sed -i 's/\r$//' ncbi_accessions_batch2.csv. Didn't work.
# Try fixing the file dos2unix sample_file.txt. -Worked!!!
# Still missing one genome. Which one!
batch1and2 <- read.delim("data/batch1and2.txt", header = F) %>%
  mutate(accession = gsub(".zip", "", V1))
batch3 <- ncbi_accessions %>%
  filter(ncbi_genbank_assembly_accession %notin% batch1and2$accession)
# Downloaded manually
# datasets download genome accession GCA_001969565.1 --filename GCA_001969565.1.zip --include genome

# Unzipped, then extracted fastas, only got 13877 .fna files
# What happened to the other 7594!?
fasta_batch1 <- read.delim("data/fasta_batch1.txt", header = F) %>%
  mutate(accession = substr(V1, start = 1, stop = 15))
batch3$ncbi_genbank_assembly_accession %in% fasta_batch1$accession # Got that one
batch4 <- ncbi_accessions %>%
  filter(ncbi_genbank_assembly_accession %notin% fasta_batch1$accession)
# Those accession have .fnas. Not sure what happened. Redownload, unzip, move .fna
# Must write as .csv, then save as .txt (tab-delim), then dos2unix for it to work!
#write.table(batch4, "data/ncbi_accessions_batch4.csv", sep = ",", row.names = F, col.names = F)

# Unzipped, then extracted fastas, got up to 21467 .fna files
# Need to get the last 5! Do manually.
fasta_batch1and2 <- read.delim("data/fasta_batch1and2.txt", header = F) %>%
  mutate(accession = substr(V1, start = 1, stop = 15))
last4 <- ncbi_accessions %>%
  filter(ncbi_genbank_assembly_accession %notin% fasta_batch1and2$accession)
# Downloaded, but these 4 did not have .fna files!!
# Checked online - those 4 were removed:
# Status: RefSeq GCF_030845235.1 is suppressed
# Status: RefSeq GCF_025379945.1 is suppressed
# Status: RefSeq GCF_026712145.1 is suppressed
# Status: RefSeq GCF_026278285.1 is suppressed
# Need to move on without those 4. Final number 21467.
# See which taxa?
gtdb_focal_removed <- gtdb_focal %>%
  filter(ncbi_genbank_assembly_accession %in% last4$ncbi_genbank_assembly_accession)
# 1 Streptomyces, 2 Gemmata, 1 Bacillus
# Final table
gtdb_focal_final <- gtdb_focal %>%
  filter(ncbi_genbank_assembly_accession %in% fasta_batch1and2$accession)
#saveRDS(gtdb_focal_final, "data/gtdb_focal_final.rds")
genus_info <- genus_info %>%
  mutate(n_Genomes_downloaded = ifelse(Genus == "Streptomyces",
                                       n_GTDB_genomes - 1,
                                       ifelse(Genus == "Gemmata",
                                              n_GTDB_genomes - 2,
                                              ifelse(Genus == "Bacillus",
                                                     n_GTDB_genomes - 1,
                                                     n_GTDB_genomes))))
#writexl::write_xlsx(genus_info, "data/genus_info_75point1.xlsx", format_headers = F)
# All 573 Streptosporangiaceae genomes downloaded fine



#### __dRep ####
# Also run Sylph those 33 genera from the dereplicated version of GTDB
# See if results are the same, especially for Bacillus
# This will be a subset of genomes already downloaded
# Just need list of accessions to copy them into a new directory to sketch
# 113104 total but 107235 bacteria
bacGT_dRep <-bacGT %>%
  filter(gtdb_representative == "t")
gtdb_focal_dRep <- bacGT_dRep %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus %in% genus_genomes_2$Genus) %>%
  filter(ncbi_genbank_assembly_accession %notin% last4$ncbi_genbank_assembly_accession)
ncbi_accessions_dRep <- gtdb_focal_dRep %>%
  dplyr::select(ncbi_genbank_assembly_accession)
# That won't work though, because there is other info in the downloaded genome names
# Get a list of the actual genome file names, then filter to the dRep
# ls Ref_genomes > genome_names.txt
genome_names <- read.delim("data/genome_names.txt", header = F) %>%
  mutate(Accession = substr(V1, start = 1, stop = 15))
genome_names_dRep <- genome_names %>%
  filter(Accession %in% ncbi_accessions_dRep$ncbi_genbank_assembly_accession)
#write.csv(genome_names_dRep$V1, "data/dRep_filenames.csv", row.names = F)
#write.table(genome_names_dRep$V1, "data/dRep_filenames.txt", sep = "\t", row.names = F, col.names = F)
# Don't forget dos2unix and to remove quotations!!
# Ran the copy script, looks like 1 is missing
dRep_moved <- read.delim("data/dRep_moved.txt", header = F)
missing <- genome_names_dRep %>%
  filter(V1 %notin% dRep_moved$V1)
# Moved it manually. Proceed to Sylph on the dRep set.
genus_genomes_dRep <- gtdb_focal_dRep %>%
  group_by(Genus) %>%
  summarise(n_GTDB_Genomes_dRep = n())



#### _Sylph ####
# Use Sylph to select reference genomes
# In terminal run Sylph on the R1 and R2 QC reads
# Sylph can do a profile at 95% ANI, or a query with 90% ANI. Try both.
# Note: Sylph cannot find super low abundance genomes. Sylph requires > 0.01-0.05x coverage
# Note: Sylph cannot reliably find genomes at genus level or higher (if not present at species level)
# Note: Sequence abundance per sample adds up to 100% just for the genomes hit.
# Note: Eff_cov is the effective coverage
# Need to isolate accession to map to GTDB taxonomy
# Need to get # genomes per sample and per focal genus
gtdb_focal_final <- readRDS("data/gtdb_focal_final.rds")
gtdb_ncbi_map <- gtdb_focal_final %>%
  dplyr::select(Domain, Phylum, Class, Order, Family, Genus, Species, 
                ncbi_genbank_assembly_accession)
gtdb_ncbi_map_strep <- rbind(gtdb_focal_final, gtdb_strep) %>%
  dplyr::select(Domain, Phylum, Class, Order, Family, Genus, Species, 
                ncbi_genbank_assembly_accession)
meta_key <- meta %>%
  dplyr::select(sampleID, vegetation_type, AI) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(-sampleID)
meta_AI <- meta_key %>%
  dplyr::select(-vegetation_type)



#### __95% ANI ####
sylph_profile <- read.delim("data/sylph_profile_ani95.tsv")
nrow(sylph_profile) # 714 total rows
length(unique(sylph_profile$Contig_name)) # 223 uniques
length(unique(sylph_profile$Genome_file)) # 223 uniques
length(unique(sylph_profile$Sample_file)) # 101 uniques - so 3 samples with none of these at 95% ANI.

sylph_profile <- read.delim("data/sylph_profile_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(Accession = substr(Genome_file, start = 48, stop = 62)) %>%
  left_join(., gtdb_ncbi_map, by = c("Accession" = "ncbi_genbank_assembly_accession")) %>%
  dplyr::select(-Sample_file, -Genome_file) %>%
  dplyr::select(SampleID, Accession, Domain, Phylum, Class, Order, Family, Genus, Species,
                everything())
length(unique(sylph_profile$Genus)) # 15 genera
unique(sylph_profile$Genus)
sylph_profile_samples <- sylph_profile %>%
  group_by(SampleID) %>%
  summarise(n_Genomes = n()) %>%
  left_join(., meta_key, by = "SampleID")
ggplot(sylph_profile_samples, aes(AI, n_Genomes)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
sylph_profile_genomes <- sylph_profile %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n())
sylph_profile_genera <- sylph_profile %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Genomes_Sylph95 = n())
sylph_profile_genera_samples <- sylph_profile %>%
  group_by(Genus, SampleID) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Samples_Sylph95 = n())
# Add this info
genus_info_toPlot <- read_xlsx("data/genus_info_75point1.xlsx") %>%
  left_join(., sylph_profile_genera, by = "Genus") %>%
  left_join(., sylph_profile_genera_samples, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph95 = 0,
                  n_Samples_Sylph95 = 0))
ggplot(genus_info_toPlot, aes(n_Genomes_downloaded, n_Genomes_Sylph95)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()



#### ___Select ####
# Need to select genomes for coverM and eventually StrainFinder
# Could just use most prevalent Sylph for those 4 most prevalent genera
# But should also check completeness and contamination
checkM <- gtdb_focal_final %>%
  dplyr::select(ncbi_genbank_assembly_accession, 
                checkm2_completeness, checkm2_contamination,
                checkm_completeness, checkm_contamination)
ggplot(checkM, aes(checkm_completeness, checkm2_completeness)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, pch = 16, alpha = 0.5) +
  geom_smooth(method = "lm") +
  scale_y_continuous(limits = c(25, 100)) +
  scale_x_continuous(limits = c(25, 100)) +
  theme_bw()
ggplot(checkM, aes(checkm_contamination, checkm2_contamination)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, pch = 16, alpha = 0.5) +
  geom_smooth(method = "lm") +
  scale_y_continuous(limits = c(0, 30)) +
  scale_x_continuous(limits = c(0, 30)) +
  theme_bw()
sylph_profile_genomes_select <- sylph_profile %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n()) %>%
  arrange(desc(n_Samples)) %>%
  filter(Genus %in% c("Bradyrhizobium", "Streptomyces", "Udaeobacter", "Mycobacterium")) %>%
  ungroup() %>%
  group_by(Genus) %>%
  slice_max(order_by = n_Samples, n = 4, with_ties = FALSE) %>%
  left_join(., checkM, by = c("Accession" = "ncbi_genbank_assembly_accession"))
# OK, the good thing is that the most prevalent is also best or almost best/comparable checkM
# Should also check competitive coverM for these 16 though

# Are those in the dRep set?
sum(sylph_profile_genomes_select$Accession %in% ncbi_accessions_dRep$ncbi_genbank_assembly_accession)
# 9. So 7 of those 16 are "strain level". I guess good to have both?
# Which ones?
sylph_profile_genomes_select_dRep <- sylph_profile_genomes_select %>%
  filter(Accession %in% ncbi_accessions_dRep$ncbi_genbank_assembly_accession)
# Brady, Udaeo, Myco

# File name list
genome_names_top16 <- genome_names %>%
  filter(Accession %in% sylph_profile_genomes_select$Accession)
#write.csv(genome_names_top16$V1, "data/top16_filenames.csv", row.names = F)



#### ___Heatmaps ####
# Now need to see if different genomes were hit for each genus across the gradient.
# Or if genomes present in multiple samples spanned a gradient
# Bradyrhizobium, Streptomyces, Udaeobacter, Mycobacterium present in >50 samples
# Bradyrhizobium
sylph_brady <- sylph_profile %>%
  dplyr::filter(Genus == "Bradyrhizobium")
no_brady <- meta %>%
  filter(sampleID %notin% sylph_brady$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_brady <- sylph_profile %>%
  filter(Genus == "Bradyrhizobium") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_brady) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_brady_mat <- sylph_brady %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_brady_ai <- sylph_brady %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_brady_mat), 
                       AI = sylph_brady_ai$AI)
pheatmap(sylph_brady_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Bradyrhizobium",
         width = 8,
         height = 8,
         filename = "InitialFigs/Sylph_Brady.png")
dev.off()
dev.set(dev.next())

# Streptomyces
sylph_strepto <- sylph_profile %>%
  filter(Genus == "Streptomyces")
no_strepto <- meta %>%
  filter(sampleID %notin% sylph_strepto$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_strepto <- sylph_profile %>%
  filter(Genus == "Streptomyces") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_strepto) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_strepto_mat <- sylph_strepto %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_strepto_ai <- sylph_strepto %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_strepto_mat), 
                       AI = sylph_strepto_ai$AI)
pheatmap(sylph_strepto_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Streptomyces",
         width = 8,
         height = 8,
         filename = "InitialFigs/Sylph_Strepto.png")
dev.off()
dev.set(dev.next())

# Udaeobacter
sylph_udaeo <- sylph_profile %>%
  filter(Genus == "Udaeobacter")
no_udaeo <- meta %>%
  filter(sampleID %notin% sylph_udaeo$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_udaeo <- sylph_profile %>%
  filter(Genus == "Udaeobacter") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_udaeo) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_udaeo_mat <- sylph_udaeo %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_udaeo_ai <- sylph_udaeo %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_udaeo_mat), 
                       AI = sylph_udaeo_ai$AI)
pheatmap(sylph_udaeo_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Udaeobacter",
         width = 8,
         height = 8,
         filename = "InitialFigs/Sylph_Udaeo.png")
dev.off()
dev.set(dev.next())

# Mycobacterium
sylph_myco <- sylph_profile %>%
  filter(Genus == "Mycobacterium")
no_myco <- meta %>%
  filter(sampleID %notin% sylph_myco$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_myco <- sylph_profile %>%
  filter(Genus == "Mycobacterium") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_myco) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_myco_mat <- sylph_myco %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_myco_ai <- sylph_myco %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_myco_mat), 
                       AI = sylph_myco_ai$AI)
pheatmap(sylph_myco_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Mycobacterium",
         width = 8,
         height = 8,
         filename = "InitialFigs/Sylph_Myco.png")
dev.off()
dev.set(dev.next())



#### __95% ANI dRep ####
# Do the same as above but with the Sylph output from the dereplicated GTDB set (n = 3281)
# Need to compare results from Sylph all vs. Sylph dRep
sylph_profile_dRep <- read.delim("data/sylph_profile_ani95_dRep.tsv")
length(unique(sylph_profile_dRep$Contig_name)) # 175 uniques
length(unique(sylph_profile_dRep$Genome_file)) # 175 uniques
length(unique(sylph_profile_dRep$Sample_file)) # 101 uniques - so 3 samples with none of these at 95% ANI.

sylph_profile_dRep <- read.delim("data/sylph_profile_ani95_dRep.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(Accession = substr(Genome_file, start = 53, stop = 67)) %>%
  left_join(., gtdb_ncbi_map, by = c("Accession" = "ncbi_genbank_assembly_accession")) %>%
  dplyr::select(-Sample_file, -Genome_file) %>%
  dplyr::select(SampleID, Accession, Domain, Phylum, Class, Order, Family, Genus, Species,
                everything())
length(unique(sylph_profile_dRep$Genus)) # 14 genera (1 less)
sylph_profile_dRep_samples <- sylph_profile_dRep %>%
  group_by(SampleID) %>%
  summarise(n_Genomes = n()) %>%
  left_join(., meta_key, by = "SampleID")
ggplot(sylph_profile_dRep_samples, aes(AI, n_Genomes)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
sylph_profile_dRep_genomes <- sylph_profile_dRep %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples_dRep = n())
sylph_profile_dRep_genera <- sylph_profile_dRep %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Genomes_Sylph95_dRep = n())
sylph_profile_dRep_genera_samples <- sylph_profile_dRep %>%
  group_by(Genus, SampleID) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Samples_Sylph95_dRep = n())
# Add this info
genus_info_dRep <- read_xlsx("data/genus_info_75point1.xlsx") %>%
  left_join(., sylph_profile_dRep_genera, by = "Genus") %>%
  left_join(., sylph_profile_dRep_genera_samples, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph95_dRep = 0,
                  n_Samples_Sylph95_dRep = 0))
ggplot(genus_info_dRep, aes(n_Genomes_downloaded, n_Genomes_Sylph95)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
genus_info_dRep <- read_xlsx("data/genus_info_75point1.xlsx") %>%
  left_join(., sylph_profile_dRep_genera, by = "Genus") %>%
  left_join(., sylph_profile_dRep_genera_samples, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph95_dRep = 0,
                  n_Samples_Sylph95_dRep = 0)) %>%
  dplyr::select(Genus, n_Genomes_Sylph95_dRep, n_Samples_Sylph95_dRep)

# Merge Rep and dRep. Compare
sylph_compare <- genus_info %>%
  left_join(., genus_info_dRep, by = "Genus") %>%
  left_join(., genus_genomes_dRep, by = "Genus")
ggplot(sylph_compare, aes(n_GTDB_genomes, n_GTDB_Genomes_dRep)) +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth(method = "lm") +
  theme_bw()
ggplot(sylph_compare, aes(n_Genomes_Sylph95, n_Genomes_Sylph95_dRep)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth(method = "lm") +
  theme_bw()
ggplot(sylph_compare, aes(n_Samples_Sylph95, n_Samples_Sylph95_dRep)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, pch = 16, alpha = 0.75) +
  geom_smooth(method = "lm") +
  theme_bw()

# Write file
sylph_compare_file <- sylph_compare %>%
  filter(n_GTDB_genomes >= 2) %>%
  dplyr::select(Genus, n_GTDB_genomes, n_GTDB_Genomes_dRep,
                n_Samples_Sylph95, n_Samples_Sylph95_dRep,
                n_Genomes_Sylph95, n_Genomes_Sylph95_dRep) %>%
  arrange(desc(n_Samples_Sylph95))
#writexl::write_xlsx(sylph_compare_file, "data/sylph_compare_drep.xlsx", format_headers = F)



#### __95% ANI w Strep ####
# Sylph output with the Streptosporangiaceae genomes included
sylph_profile <- read.delim("data/sylph_profile_strep_ani95.tsv")
nrow(sylph_profile) # 771 total rows
length(unique(sylph_profile$Contig_name)) # 245 uniques
length(unique(sylph_profile$Genome_file)) # 245 uniques
length(unique(sylph_profile$Sample_file)) # 101 uniques - so 3 samples with none of these at 95% ANI.

sylph_profile <- read.delim("data/sylph_profile_strep_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(Accession = substr(Genome_file, start = 21, stop = 35)) %>%
  left_join(., gtdb_ncbi_map_strep, by = c("Accession" = "ncbi_genbank_assembly_accession")) %>%
  dplyr::select(-Sample_file, -Genome_file) %>%
  dplyr::select(SampleID, Accession, Domain, Phylum, Class, Order, Family, Genus, Species,
                everything())
length(unique(sylph_profile$Genus)) # 26 genera
unique(sylph_profile$Genus)
length(unique(sylph_profile$Family)) # 14 families
sylph_profile_strep <- sylph_profile %>%
  filter(Family == "Streptosporangiaceae")
nrow(sylph_profile_strep) # 57 total rows
length(unique(sylph_profile_strep$Accession)) # 22 unique Streptosporangiaceae genomes found
length(unique(sylph_profile_strep$SampleID)) # 38 unique samples
sylph_profile_samples <- sylph_profile_strep %>%
  group_by(SampleID) %>%
  summarise(n_Genomes = n()) %>%
  left_join(., meta_key, by = "SampleID")
ggplot(sylph_profile_samples, aes(AI, n_Genomes)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
sylph_profile_genomes <- sylph_profile_strep %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n())
# The most samples per individual accession is 9 for a Trebonia and a Palsa-504
sylph_profile_genera <- sylph_profile_strep %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Genomes_Sylph95 = n())
# There are 11 different Streptosporangiaceae genera represented!
sylph_profile_genera_samples <- sylph_profile_strep %>%
  group_by(Genus, SampleID) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Samples_Sylph95 = n())
# Palsa-504 in 14 samples, Trebonia in 11 samples
# However, neither of these nor any Streptosporangiaceae were in the abund ubiq genera
# Add this info
genus_info <- read_xlsx("data/genus_info_75point1.xlsx") %>%
  left_join(., sylph_profile_genera, by = "Genus") %>%
  left_join(., sylph_profile_genera_samples, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph95 = 0,
                  n_Samples_Sylph95 = 0))
ggplot(genus_info, aes(n_Genomes_downloaded, n_Genomes_Sylph95)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()



#### __90% ANI ####
sylph_query <- read.delim("data/sylph_query_ani90.tsv")
length(unique(sylph_query$Contig_name)) # 16791 uniques
length(unique(sylph_query$Genome_file)) # 16791 uniques
length(unique(sylph_query$Sample_file)) # 104 uniques

sylph_query <- read.delim("data/sylph_query_ani90.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(Accession = substr(Genome_file, start = 48, stop = 62)) %>%
  left_join(., gtdb_ncbi_map, by = c("Accession" = "ncbi_genbank_assembly_accession")) %>%
  dplyr::select(-Sample_file, -Genome_file) %>%
  dplyr::select(SampleID, Accession, Domain, Phylum, Class, Order, Family, Genus, Species,
                everything())
length(unique(sylph_query$Genus)) # 29 genera
sylph_query_samples <- sylph_query %>%
  group_by(SampleID) %>%
  summarise(n_Genomes = n()) %>%
  left_join(., meta_key, by = "SampleID")
ggplot(sylph_query_samples, aes(AI, n_Genomes)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
sylph_query_genomes <- sylph_query %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n())
sylph_query_genera <- sylph_query %>%
  group_by(Genus, Accession) %>%
  summarise(n_Samples = n()) %>%
  ungroup() %>%
  group_by(Genus) %>%
  summarise(n_Genomes_Sylph90 = n())
# Add this info
genus_info <- genus_info %>%
  left_join(., sylph_query_genera, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph90 = 0))
ggplot(genus_info, aes(n_Genomes_downloaded, n_Genomes_Sylph90)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
ggplot(genus_info, aes(n_Genomes_Sylph90, n_Genomes_Sylph95)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
ggplot(genus_info, aes(n_OTUs, n_GTDB_genomes)) +
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  theme_bw()
genus_info_top4 <- genus_info %>%
  filter(Genus %in% c("Bradyrhizobium", "Streptomyces", "Udaeobacter", "Mycobacterium")) %>%
  mutate(Present_Perc_OTU = Present_Perc,
         Present_Perc_Sylph95 = n_Samples_Sylph95/104*100) %>% 
  dplyr::select(Phylum, Genus, n_OTUs, Present_Perc_OTU, n_GTDB_genomes, 
                n_Genomes_Sylph95, Present_Perc_Sylph95)



#### ___Heatmaps ####
# Now need to see if different genomes were hit for each genus across the gradient.
# Or if genomes present in multiple samples spanned a gradient
# Bradyrhizobium, Streptomyces, Udaeobacter, Mycobacterium present in >50 samples
# Bradyrhizobium
sylph_brady_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Bradyrhizobium")
no_brady <- meta %>%
  filter(sampleID %notin% sylph_brady_dRep$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_brady_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Bradyrhizobium") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_brady) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_brady_dRep_mat <- sylph_brady_dRep %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_brady_dRep_ai <- sylph_brady_dRep %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_brady_dRep_mat), 
                       AI = sylph_brady_dRep_ai$AI)
pheatmap(sylph_brady_dRep_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Bradyrhizobium",
         width = 8,
         height = 8,
         filename = "InitialFigs/sylph_brady_dRep.png")
dev.off()
dev.set(dev.next())

# Streptomyces
sylph_strepto_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Streptomyces")
no_strepto <- meta %>%
  filter(sampleID %notin% sylph_strepto_dRep$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_strepto_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Streptomyces") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_strepto) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_strepto_dRep_mat <- sylph_strepto_dRep %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_strepto_dRep_ai <- sylph_strepto_dRep %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_strepto_dRep_mat), 
                       AI = sylph_strepto_dRep_ai$AI)
pheatmap(sylph_strepto_dRep_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Streptomyces",
         width = 8,
         height = 8,
         filename = "InitialFigs/sylph_strepto_dRep.png")
dev.off()
dev.set(dev.next())

# Udaeobacter
sylph_udaeo_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Udaeobacter")
no_udaeo <- meta %>%
  filter(sampleID %notin% sylph_udaeo_dRep$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_udaeo_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Udaeobacter") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_udaeo) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_udaeo_dRep_mat <- sylph_udaeo_dRep %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_udaeo_dRep_ai <- sylph_udaeo_dRep %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_udaeo_dRep_mat), 
                       AI = sylph_udaeo_dRep_ai$AI)
pheatmap(sylph_udaeo_dRep_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Udaeobacter",
         width = 8,
         height = 8,
         filename = "InitialFigs/sylph_udaeo_dRep.png")
dev.off()
dev.set(dev.next())

# Mycobacterium
sylph_myco_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Mycobacterium")
no_myco <- meta %>%
  filter(sampleID %notin% sylph_myco_dRep$SampleID) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  mutate(Accession = NA,
         Eff_cov = 0) %>%
  dplyr::select(SampleID, Accession, Eff_cov)
sylph_myco_dRep <- sylph_profile_dRep %>%
  filter(Genus == "Mycobacterium") %>%
  dplyr::select(SampleID, Accession, Eff_cov) %>%
  rbind(., no_myco) %>%
  pivot_wider(names_from = Accession, values_from = Eff_cov) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-`NA`) %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI)
sylph_myco_dRep_mat <- sylph_myco_dRep %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
sylph_myco_dRep_ai <- sylph_myco_dRep %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(sylph_myco_dRep_mat), 
                       AI = sylph_myco_dRep_ai$AI)
pheatmap(sylph_myco_dRep_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         main = "Mycobacterium",
         width = 8,
         height = 8,
         filename = "InitialFigs/sylph_myco_dRep.png")
dev.off()
dev.set(dev.next())


#### _FastANI ####
# For top 4 genera, check ANI of top 4 genomes
genus_accession <- gtdb_focal_final %>%
  dplyr::select(Genus, ncbi_genbank_assembly_accession) %>%
  filter(Genus %in% c("Bradyrhizobium", "Streptomyces", "Udaeobacter", "Mycobacterium"))
ani <- read_xlsx("data/ani_top16.xlsx") %>%
  mutate(Query = substr(Query, start = 1, stop = 15),
         Reference = substr(Reference, start = 1, stop = 15)) %>%
  left_join(., genus_accession, by = c("Query" = "ncbi_genbank_assembly_accession")) %>%
  rename(Genus_Query = Genus) %>%
  left_join(., genus_accession, by = c("Reference" = "ncbi_genbank_assembly_accession")) %>%
  rename(Genus_Reference = Genus) %>%
  mutate(Comparison = ifelse(Genus_Query == Genus_Reference,
                             "Same",
                             "Different")) %>%
  filter(Comparison == "Same") %>%
  arrange(Genus_Query, desc(`ANI Estimate`)) %>%
  mutate(Genus = Genus_Query) %>%
  dplyr::select(Genus, Query, Reference, `ANI Estimate`)
ani$pair_id <- apply(ani, 1, function(row) paste(sort(row), collapse = "_"))
ani$pair_id <- substr(ani$pair_id, start = 9, stop = nchar(ani$pair_id))
ani <- ani[!duplicated(ani$pair_id), ]
ani <- ani %>%
  dplyr::select(Genus, Query, Reference, `ANI Estimate`)
length(unique(ani$`ANI Estimate`))
length(unique(ani$Query))
length(unique(ani$Reference))
length(unique(c(ani$Query, ani$Reference)))



#### _coverM ####
# Look at coverM mean coverage results for different genomes
# Need to loop through the output, merge into one, then merge other data
# Top4 (most prevalent of the top 4 genera)
# Compare those to Sylph Effective Coverage
# (and I guess Sylph dRep Effective Coverage)

#### __Mean cov ####
list.files("data/coverm_output/")
cm_top4_dfs <- list()
for (i in 1:length(list.files("data/coverm_output/"))) {
  filename <- list.files("data/coverm_output/")[i]
  cm_top4_dfs[[i]] <- read.delim(paste("data/coverm_output/", filename, sep = ""))
}
cm_top4 <- cm_top4_dfs[[1]]
for (i in 2:length(list.files("data/coverm_output/"))) {
  cm_top4 <- cbind(cm_top4, cm_top4_dfs[[i]][, 2, drop = FALSE])
}
cm_n_Samples <- cm_top4 %>%
  column_to_rownames(var = "Genome") %>%
  mutate(n_Samples_coverM = rowSums(. != 0)) %>%
  rownames_to_column(var = "Genome") %>%
  dplyr::select(Genome, n_Samples_coverM)
cm_top4_table <- cm_top4 %>%
  mutate(Genome = substr(Genome, start = 1, stop = 15)) %>%
  column_to_rownames(var = "Genome") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(SampleID = gsub("X", "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  set_names(paste(names(.), "_coverM", sep = "")) %>%
  rename(SampleID = SampleID_coverM)
sylph_brady_top <- sylph_brady %>%
  dplyr::select(SampleID, GCA_016616885.1) %>%
  set_names(c("SampleID", "GCA_016616885.1_Sylph"))
sylph_strepto_top <- sylph_strepto %>%
  dplyr::select(SampleID, GCA_012922115.1) %>%
  set_names(c("SampleID", "GCA_012922115.1_Sylph"))
sylph_udaeo_top <- sylph_udaeo %>%
  dplyr::select(SampleID, GCA_013372845.1) %>%
  set_names(c("SampleID", "GCA_013372845.1_Sylph"))
sylph_myco_top <- sylph_myco %>%
  dplyr::select(SampleID, GCA_019668465.1) %>%
  set_names(c("SampleID", "GCA_019668465.1_Sylph"))
compare_cm_sylph <- cm_top4_table %>%
  left_join(., sylph_brady_top, by = "SampleID") %>%
  left_join(., sylph_strepto_top, by = "SampleID") %>%
  left_join(., sylph_udaeo_top, by = "SampleID") %>%
  left_join(., sylph_myco_top, by = "SampleID") %>%
  left_join(., meta_AI, by = "SampleID") %>%
  arrange(AI) %>%
  dplyr::select(SampleID,
                GCA_013372845.1_coverM, GCA_013372845.1_Sylph,
                GCA_019668465.1_coverM, GCA_019668465.1_Sylph,
                GCA_016616885.1_coverM, GCA_016616885.1_Sylph,
                GCA_012922115.1_coverM, GCA_012922115.1_Sylph,
                AI)
compare_cm_sylph_mat <- compare_cm_sylph %>%
  dplyr::select(-AI) %>%
  column_to_rownames(var = "SampleID") %>%
  t()
compare_cm_sylph_ai <- compare_cm_sylph %>%
  dplyr::select(SampleID, AI)
ann_cols <- data.frame(row.names = colnames(compare_cm_sylph_mat), 
                       AI = compare_cm_sylph_ai$AI)
ann_rows <- data.frame(row.names = rownames(compare_cm_sylph_mat), 
                       Genus = c("Udaeobacter", "Udaeobacter",
                                 "Mycobacterium", "Mycobacterium",
                                 "Bradyrhizobium", "Bradyrhizobium",
                                 "Streptomyces", "Streptomyces")) %>%
  mutate(Genus = factor(Genus, levels = c("Udaeobacter", "Mycobacterium", 
                                          "Bradyrhizobium", "Streptomyces")))
pheatmap(compare_cm_sylph_mat,
         border_color = NA,
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 315,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         main = "Sylph vs. coverM",
         width = 8,
         height = 8,
         filename = "InitialFigs/compare_sylph_coverm.png")
dev.off()
dev.set(dev.next())

# Need to get list of samples to run StrainFinder
# Run on samples with the highest coverage
# Set mindepth parameter to half
cm_top4_table
hist(cm_top4_table$GCA_013372845.1_coverM)
hist(cm_top4_table$GCA_019668465.1_coverM)
hist(cm_top4_table$GCA_016616885.1_coverM)
hist(cm_top4_table$GCA_012922115.1_coverM)

# GCA_013372845.1_coverM has 15 samples > 4 # Udaeobacter
# GCA_019668465.1_coverM has 10 samples > 4 # Mycobacterium
# GCA_016616885.1_coverM has 57 samples > 4 # Bradyrhizobium
# GCA_012922115.1_coverM has 0 samples > 4 # Streptomyces

# Clearly, start with GCA_016616885.1 (Bradyrhizobium)
brady_cov4 <- cm_top4_table %>%
  filter(GCA_016616885.1_coverM > 4) %>%
  left_join(meta_AI, by = "SampleID") %>%
  mutate(HalfMean = round(GCA_016616885.1_coverM/2, digits = 0),
         UnderTwiceMean = round((GCA_016616885.1_coverM*2) - 1, digits = 0)) %>%
  rename(MeanCoverage = GCA_016616885.1_coverM) %>%
  dplyr::select(SampleID, AI, MeanCoverage, HalfMean, UnderTwiceMean)
myco_cov4 <- cm_top4_table %>%
  filter(GCA_019668465.1_coverM > 4) %>%
  left_join(meta_AI, by = "SampleID") %>%
  mutate(HalfMean = round(GCA_019668465.1_coverM/2, digits = 0),
         UnderTwiceMean = round((GCA_019668465.1_coverM*2) - 1, digits = 0)) %>%
  rename(MeanCoverage = GCA_019668465.1_coverM) %>%
  dplyr::select(SampleID, AI, MeanCoverage, HalfMean, UnderTwiceMean)



#### __Cov hist ####
# Note: Sample 401590 had no reads of any of the top 4 mapped so deleted file
# This is weird because for mean coverage and rel abund it did...
# Now only 103 samples
# Do for Brady first to keep file sizes down
# Then do for Myco so you can run StrainFinder to compare to MAG

# Brady
list.files("data/coverm_output_ch/")
cm_brady_ch_dfs <- list()
for (i in 1:length(list.files("data/coverm_output_ch/"))) {
  filename <- list.files("data/coverm_output_ch/")[i]
  cm_brady_ch_dfs[[i]] <- read.delim(paste("data/coverm_output_ch/", filename, sep = "")) %>%
    filter(Genome == "GCA_016616885.1_ASM1661688v1_genomic")
}

# Print number of rows by sample
brady_nrows <- as.data.frame(matrix(data = NA, 
                                    nrow = length(list.files("data/coverm_output_ch/")),
                                    ncol = 2)) %>%
  set_names(c("Sample", "BradyCoverages"))
for (i in 1:length(list.files("data/coverm_output_ch/"))) {
  brady_nrows$Sample[i] <- list.files("data/coverm_output_ch/")[i]
  brady_nrows$BradyCoverages[i] <- nrow(cm_brady_ch_dfs[[i]])
}
# Sample 12446 didn't have any Brady either but that one wasn't in brady_cov4

cm_brady_ch <- cm_brady_ch_dfs[[1]]
for (i in 2:length(list.files("data/coverm_output_ch/"))) {
  cm_brady_ch <- rbind(cm_brady_ch, cm_brady_ch_dfs[[i]])
}
cm_brady_ch <- cm_brady_ch %>%
  separate(Sample, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  filter(SampleID %in% brady_cov4$SampleID)
length(unique(cm_brady_ch$SampleID)) # 56, good

# Check the first 2 samples
# Check histogram
# Confirm that the total number of bases is the genome size (8,085,095)
# Missing 150 bp - that's because coverM trims the first and last 75 bp of the genome by default
cm_brady_ch_12816 <- cm_brady_ch %>%
  filter(SampleID == "12816")
plot(cm_brady_ch_12816$Coverage, cm_brady_ch_12816$Bases)
sum(cm_brady_ch_12816$Bases)

cm_brady_ch_138538 <- cm_brady_ch %>%
  filter(SampleID == "138538")
plot(cm_brady_ch_138538$Coverage, log10(cm_brady_ch_138538$Bases+1))
plot(cm_brady_ch_138538$Bases, cm_brady_ch_138538$Coverage)
sum(cm_brady_ch_138538$Bases)

# Check the number of bases with 0 coverage
# For example, see how many samples have < 50% with 0
# 34 samples have > 50% genome with at least coverage of 1
cm_brady_ch_zero <- cm_brady_ch %>%
  filter(Coverage == 0) %>%
  mutate(PercGenomeZero = round(Bases/8084945*100, digits = 2)) %>%
  dplyr::select(SampleID, PercGenomeZero)

# Check the max coverages per sample
cm_brady_ch_max <- cm_brady_ch %>%
  group_by(SampleID) %>%
  summarise(MaxCov = max(Coverage)) %>%
  mutate(HalfMax = MaxCov/2)
  
# How about setting the max at the first value that has only 10 or less bases
cm_brady_ch_cut10 <- cm_brady_ch %>%
  group_by(SampleID) %>%
  filter(Bases < 10) %>%
  slice_head(n = 1) %>%
  rename(Coverage10Bases = Coverage) %>%
  dplyr::select(SampleID, Coverage10Bases)

# Merge all of this info
brady_info <- brady_cov4 %>%
  filter(SampleID %in% unique(cm_brady_ch$SampleID)) %>%
  left_join(cm_brady_ch_zero, by = "SampleID") %>%
  left_join(cm_brady_ch_max, by = "SampleID") %>%
  left_join(cm_brady_ch_cut10, by = "SampleID") %>%
  filter(PercGenomeZero < 60) %>%
  arrange(AI)
#write_xlsx(brady_info, "data/Bradyrhizobium_StrainFinder_Info.xlsx", format_headers = F)

ggplot(brady_info, aes(MeanCoverage, PercGenomeZero)) + 
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  geom_smooth(method = "lm") +
  labs(x = "Mean Coverage (coverM)",
       y = "% bases with 0 coverage (coverM)") +
  theme_bw()

ggplot(brady_info, aes(AI, MeanCoverage)) + 
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  geom_smooth(method = "lm") +
  labs(x = "Aridity index",
       y = "Mean Coverage (coverM)") +
  theme_bw()



#### ___Myco ####
list.files("data/coverm_output_ch/")
cm_myco_ch_dfs <- list()
for (i in 1:length(list.files("data/coverm_output_ch/"))) {
  filename <- list.files("data/coverm_output_ch/")[i]
  cm_myco_ch_dfs[[i]] <- read.delim(paste("data/coverm_output_ch/", filename, sep = "")) %>%
    filter(Genome == "GCA_019668465.1_ASM1966846v1_genomic")
}

# Print number of rows by sample
myco_nrows <- as.data.frame(matrix(data = NA, 
                                    nrow = length(list.files("data/coverm_output_ch/")),
                                    ncol = 2)) %>%
  set_names(c("Sample", "mycoCoverages"))
for (i in 1:length(list.files("data/coverm_output_ch/"))) {
  myco_nrows$Sample[i] <- list.files("data/coverm_output_ch/")[i]
  myco_nrows$mycoCoverages[i] <- nrow(cm_myco_ch_dfs[[i]])
}

cm_myco_ch <- cm_myco_ch_dfs[[1]]
for (i in 2:length(list.files("data/coverm_output_ch/"))) {
  cm_myco_ch <- rbind(cm_myco_ch, cm_myco_ch_dfs[[i]])
}
cm_myco_ch <- cm_myco_ch %>%
  separate(Sample, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  filter(SampleID %in% myco_cov4$SampleID)
length(unique(cm_myco_ch$SampleID)) # 10, good

# Check histogram
# Confirm that the total number of bases is the genome size (5,831,001)
# Missing 150 bp - that's because coverM trims the first and last 75 bp of the genome by default
table(cm_myco_ch$SampleID)
cm_myco_ch_401606 <- cm_myco_ch %>%
  filter(SampleID == "401606")
plot(cm_myco_ch_401606$Coverage, cm_myco_ch_401606$Bases)
sum(cm_myco_ch_401606$Bases)

# Check the number of bases with 0 coverage
# For example, see how many samples have < 50% with 0
# 34 samples have > 50% genome with at least coverage of 1
cm_myco_ch_zero <- cm_myco_ch %>%
  filter(Coverage == 0) %>%
  mutate(PercGenomeZero = round(Bases/5831001*100, digits = 2)) %>%
  dplyr::select(SampleID, PercGenomeZero)

# Check the max coverages per sample
cm_myco_ch_max <- cm_myco_ch %>%
  group_by(SampleID) %>%
  summarise(MaxCov = max(Coverage)) %>%
  mutate(HalfMax = MaxCov/2)

# How about setting the max at the first value that has only 10 or less bases
# Use this! (Dylan Chivian confirmed it's a good idea)
cm_myco_ch_cut10 <- cm_myco_ch %>%
  group_by(SampleID) %>%
  filter(Bases < 10) %>%
  slice_head(n = 1) %>%
  rename(Coverage10Bases = Coverage) %>%
  dplyr::select(SampleID, Coverage10Bases)

# Merge all of this info
myco_info <- myco_cov4 %>%
  filter(SampleID %in% unique(cm_myco_ch$SampleID)) %>%
  left_join(cm_myco_ch_zero, by = "SampleID") %>%
  left_join(cm_myco_ch_max, by = "SampleID") %>%
  left_join(cm_myco_ch_cut10, by = "SampleID") %>%
  filter(PercGenomeZero < 60) %>%
  arrange(AI)
write_xlsx(myco_info, "data/Mycobacteria_StrainFinder_Info.xlsx", format_headers = F)

ggplot(myco_info, aes(MeanCoverage, PercGenomeZero)) + 
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  geom_smooth(method = "lm") +
  labs(x = "Mean Coverage (coverM)",
       y = "% bases with 0 coverage (coverM)") +
  theme_bw()

ggplot(myco_info, aes(AI, MeanCoverage)) + 
  geom_point(size = 2, pch = 16, alpha = 0.8) +
  geom_smooth(method = "lm") +
  labs(x = "Aridity index",
       y = "Mean Coverage (coverM)") +
  theme_bw()




#### __Rel abund ####
list.files("data/coverm_output_ra/")
cm_top4_ra_dfs <- list()
for (i in 1:length(list.files("data/coverm_output_ra/"))) {
  filename <- list.files("data/coverm_output_ra/")[i]
  cm_top4_ra_dfs[[i]] <- read.delim(paste("data/coverm_output_ra/", filename, sep = ""))
}
cm_top4_ra <- cm_top4_ra_dfs[[1]]
for (i in 2:length(list.files("data/coverm_output_ra/"))) {
  cm_top4_ra <- cbind(cm_top4_ra, cm_top4_ra_dfs[[i]][, 2, drop = FALSE])
}
cm_top4_ra_table <- cm_top4_ra %>%
  mutate(Genome = substr(Genome, start = 1, stop = 15)) %>%
  column_to_rownames(var = "Genome") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(SampleID = gsub("X", "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  set_names(paste(names(.), "_coverM", sep = "")) %>%
  rename(SampleID = SampleID_coverM)



#### 5. Strain Genomics ####
# Run StrainFinder on KBase
# Use the mean coverage and coverage histogram to inform the parameters
# Run for Bradyrhizobium on samples with mean coverage > 4
# First get metadata for those 53 samples
# Plot map and get distance matrix

# CheckM
checkM <- read.delim("data/CheckM_summary_table_Brady53.tsv")
range(checkM$Completeness)
range(checkM$Contamination)
checkM2to4 <- read.delim("data/CheckM_summary_table_Brady2to4.tsv")
range(checkM2to4$Completeness)
range(checkM2to4$Contamination)
checkM_all <- rbind(checkM, checkM2to4) %>%
  mutate(Bin.Name = gsub("-Reads_1-", "", Bin.Name)) %>%
  mutate(Bin.Name = gsub(".Genome", "", Bin.Name)) %>%
  mutate(Bin.Name = gsub("Strain", "_Strain", Bin.Name)) %>%
  mutate(Strain = substr(Bin.Name, 
                         start = nchar(Bin.Name) - 7, 
                         stop = nchar(Bin.Name))) %>%
  separate(Bin.Name, remove = F, sep = "_", into = c("Brady", "sampleID",
                                                     "Junk1", "Junk2")) %>%
  dplyr::select(-Brady, -Junk1, -Junk2) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  pivot_longer(cols = c("Completeness", "Contamination"))
figS2 <- ggplot(checkM_all, aes(Strain, value, group = sampleID)) +
  geom_point(size = 2, alpha = 0.75, pch = 16, 
             position = position_dodge(width = 0.1)) +
  geom_line(linewidth = 0.1,
            position = position_dodge(width = 0.1)) +
  labs(x = "StrainFinder Strain",
       y = "%") +
  facet_wrap(~ name, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
figS2
pdf("FinalFigs/FigureS2.pdf", width = 5, height = 7)
figS2
dev.off()
png("FinalFigs/FigureS2.png", width = 5, height = 7, units = "in", res = 300)
figS2
dev.off()

# Strain RelAbund from StrainFinder
strain_rel <- read_xlsx("data/Bradyrhizobium_StrainFinder_Info.xlsx") %>%
  dplyr::select(SampleID, Strain1, Strain2, Strain3, Strain4) %>%
  pivot_longer(cols = c("Strain1", "Strain2", "Strain3", "Strain4"))
pdf("InitialFigs/Strains_RelAbund.pdf", width = 7, height = 5)
ggplot(strain_rel, aes(name, value, group = SampleID)) +
  geom_point(size = 2, alpha = 0.75, pch = 16, 
             position = position_dodge(width = 0.2)) +
  geom_line(linewidth = 0.05,
            position = position_dodge(width = 0.2)) +
  labs(x = "StrainFinder Strain",
       y = "Relative abundance (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()
strain_rel_top2 <- strain_rel %>%
  filter(name %in% c("Strain1", "Strain2")) %>%
  group_by(SampleID) %>%
  summarise(Top2Rel = sum(value))
range(strain_rel_top2$Top2Rel) # 81.7 to 83.9%

# Number of missing bac120 markers
bac120 <- read.delim("data/brady_1081_markers_summary.tsv") %>%
  mutate(name2 = gsub(".fna", "", name)) %>%
  mutate(name2 = gsub("Brady_", "", name2)) %>%
  mutate(name2 = gsub("Strain_", "Strain", name2)) %>%
  separate(name2, remove = F, into = c("sampleID", "strainID")) %>%
  mutate(strainID = ifelse(grepl("Strain", strainID) == TRUE,
                           strainID, "GTDB")) %>%
  mutate(strainID = factor(strainID,
                           levels = c("Strain1", "Strain2", "Strain3", "Strain4",
                                       "GTDB")))
pdf("InitialFigs/Strains_MissingBac120.pdf", width = 7, height = 5)
ggplot(bac120, aes(strainID, number_missing_genes)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, pch = 16, width = 0.25) +
  labs(x = "Genome",
       y = "Number of missing bac120 genes") +
  ggtitle("Bradyrhizobium genomes, n = 1081") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()
bac120_strains <- bac120 %>%
  filter(grepl("Strain_1|Strain_2", name))
range(bac120_strains$number_missing_genes)

# Relationship between completeness and missings
checkM_all <- rbind(checkM, checkM2to4) %>%
  mutate(Bin.Name = gsub("-Reads_1-", "", Bin.Name)) %>%
  mutate(Bin.Name = gsub(".Genome", "", Bin.Name)) %>%
  mutate(Bin.Name = gsub("Strain", "_Strain", Bin.Name)) %>%
  mutate(Strain = substr(Bin.Name, 
                         start = nchar(Bin.Name) - 7, 
                         stop = nchar(Bin.Name))) %>%
  separate(Bin.Name, remove = F, sep = "_", into = c("Brady", "sampleID",
                                                     "Junk1", "Junk2")) %>%
  dplyr::select(-Brady, -Junk1, -Junk2) %>%
  mutate(sampleID = as.character(sampleID))
miscom <- bac120 %>%
  mutate(Bin.Name = gsub(".fna", "", name)) %>%
  left_join(., checkM_all, by = "Bin.Name") %>%
  drop_na()
plot(miscom$Completeness, miscom$number_missing_genes)

# Get genomes missing 4 or less bac120 genes, remove Strain 3 too
bac120_sub <- bac120 %>%
  filter(number_missing_genes <= 4) %>%
  filter(strainID != "Strain3")
ggplot(bac120_sub, aes(strainID, number_missing_genes)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, pch = 16, width = 0.25, height = 0) +
  labs(x = "Genome",
       y = "Number of missing bac120 genes") +
  ggtitle("Bradyrhizobium genomes, n = 915") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
brady915_alignment <- microseq::readFasta("data/Brady1081_msa.fasta") %>%
  filter(Header %in% bac120_sub$name)
#microseq::writeFasta(brady993_alignment, "data/Brady993_msa.fasta")
brady915 <- brady915_alignment %>%
  dplyr::select(Header) %>%
  mutate(Header = gsub(".fna", ".fna.gz", Header))
# write.table(brady915$Header,
#             "data/brady915.txt", sep = "\t", row.names = F, col.names = F)
#Don't forget to run dos2unix

# Check comp/cont for those GTDB
brady915_ncbi <- brady915 %>%
  mutate(Header = gsub(".fna.gz", "", Header)) %>%
  mutate(Header = gsub("Brady_Reference", "GCA_016616885.1", Header))
bacGT_915 <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% brady915_ncbi$Header) %>%
  mutate(GenomeSize = 100 * genome_size / checkm_completeness)
range(bacGT_915$checkm_completeness)
range(bacGT_915$checkm_contamination)

brady993_alignmentNG <- microseq::readFasta("data/brady993_99perc_msa.fasta")

# Continue
ref_ani <- read.delim("data/Brady_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                      "", Query)) %>%
  filter(grepl("Reference", Reference)) %>%
  filter(Query != "BradyReference.fna") %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID)) %>%
  dplyr::select(sampleID, ANI)
range(ref_ani$ANI)

d_brady <- d %>%
  filter(sampleID %in% ref_ani$sampleID) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID)
#write.table(d_brady$sampleID, "data/sampleIDs_n53.txt", sep = "\t")

# ANI for strains 1 and 2 and reference (report this in paper)
ref_ani_106 <- read.delim("data/strainfinder_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia_copy/StrainFinder/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia_copy/StrainFinder/",
                      "", Query)) %>%
  filter(grepl("Reference", Reference)) %>%
  filter(Query != "Brady_Reference.fna.gz") %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna.gz", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID)) %>%
  dplyr::select(sampleID, ANI) %>%
  drop_na()
range(ref_ani_106$ANI)
ani_106 <- read.delim("data/strainfinder_fastani.txt",
                          header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia_copy/StrainFinder/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia_copy/StrainFinder/",
                      "", Query)) %>%
  filter(!grepl("Reference", Reference)) %>%
  filter(Query != "Brady_Reference.fna.gz") %>%
  filter(Reference != Query) %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna.gz", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID))
range(ani_106$ANI)

# Map
coords <- d_brady %>%
  select(sampleID, latitude, longitude, AI)
coords_trans <- st_as_sf(coords, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
ozmap()
sf_oz <- ozmap("states")
pdf("InitialFigs/Mapn53.pdf", width = 7, height = 5)
ggplot(data = sf_oz) + 
  geom_sf(fill = "grey80", color = "white") +
  geom_sf(data = coords_trans,
          aes(color = AI)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal()
dev.off()
# Okay unfortunately these 53 with coverage > 4 are pretty clustered

# Geographic distance
dist.geog <- geosphere::distm(cbind(d_brady$longitude, d_brady$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- d_brady$sampleID
colnames(dist.geog) <- d_brady$sampleID
dist.geog[upper.tri(dist.geog, diag = TRUE)] <- NA
hist(dist.geog)

# Aridity distance
dist.ai <- as.matrix(dist(d_brady$AI, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.ai) <- d_brady$sampleID
colnames(dist.ai) <- d_brady$sampleID
dist.ai[upper.tri(dist.ai, diag = TRUE)] <- NA
hist(dist.ai)

# Environmental distance
d_brady_env <- d_brady %>%
  select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_brady_env)) {
  n_na[i] <- sum(is.na(d_brady_env[,i]))
}
n_na
d_brady_env <- d_brady_env %>%
  select(where(~ all(!is.na(.)))) %>%
  select(-sampleID, -latitude, -longitude)
# Conductivity, nitrate, carbon, temp, precip, AI
dist.env <- as.matrix(dist(d_brady_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- d_brady$sampleID
colnames(dist.env) <- d_brady$sampleID
dist.env[upper.tri(dist.env, diag = TRUE)] <- NA
hist(dist.env)



#### _Pairwise ANI/16S ####
# Need to check strain vs. reference ANI and full length 16S
# Need to check pairwise strain ANI and full length 16S
# Relate these to pairwise geographic and environmental distance

#### __ANI ####
ref_ani <- read.delim("data/Brady_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Query)) %>%
  filter(grepl("Reference", Reference)) %>%
  filter(Query != "BradyReference.fna") %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID)) %>%
  dplyr::select(sampleID, ANI) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID)
range(ref_ani$ANI) # Range is 95.4646 to 96.9770
hist(ref_ani$ANI)

d_pw <- d %>% 
  left_join(., ref_ani, by = "sampleID")
ggplot(d_pw, aes(AI, ANI)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  labs(x = "drier <= Aridity index => wetter",
       y = "ANI to reference (%)") +
  theme_bw()

strains_ani <- read.delim("data/Brady_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                      "", Query)) %>%
  filter(!grepl("Reference", Reference)) %>%
  filter(Query != "BradyReference.fna") %>%
  separate(Query, into = c("Brady", "Query"), sep = "_") %>%
  mutate(Query = gsub(".fna", "", Query)) %>%
  mutate(Query = as.character(Query)) %>%
  separate(Reference, into = c("Brady2", "Reference"), sep = "_") %>%
  mutate(Reference = gsub(".fna", "", Reference)) %>%
  mutate(Reference = as.character(Reference)) %>%
  dplyr::select(Reference, Query, ANI) %>%
  pivot_wider(names_from = Reference,
              values_from = ANI) %>%
  arrange(Query) %>%
  column_to_rownames(var = "Query") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Reference") %>%
  arrange(Reference) %>%
  column_to_rownames(var = "Reference")
dist.ani <- as.matrix(strains_ani)
dist.ani[upper.tri(dist.ani, diag = TRUE)] <- NA
hist(dist.ani)
range(dist.ani, na.rm = TRUE) # Range is 95.1578 97.0956



#### __16S ####
# 181 Sylph detected. 2 missing so 179.
pw_16S <- read.table("data/pairwise_identity_16S_181.txt",
                     header = TRUE) %>%
  column_to_rownames(var = "ID")
names(pw_16S) <- gsub(".fna", "", names(pw_16S))
names(pw_16S) <- gsub("Brady_", "", names(pw_16S))
rownames(pw_16S) <- gsub(".fna", "", rownames(pw_16S))
rownames(pw_16S) <- gsub("Brady_", "", rownames(pw_16S))
range(pw_16S) # 87.73% to 100%
strains_detected <- data.frame(GenomeID = brady181_tree$tip.label) %>%
  filter(grepl("Strain", GenomeID))
ref_16S <- pw_16S %>%
  dplyr::select(all_of(strains_detected$GenomeID), Reference) %>%
  filter(rownames(.) %in% c(strains_detected$GenomeID, "Reference")) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(strains_detected$GenomeID), Reference) %>%
  mutate_all(function(x) x * 100) %>%
  filter(rownames(.) != "Reference") %>%
  rownames_to_column(var = "sampleID") %>%
  rename("rRNA16S" = "Reference") %>%
  dplyr::select(sampleID, rRNA16S)
sum(ref_ani$sampleID != ref_16S$sampleID)
range(ref_16S$rRNA16S) # Range is 90.54 to 97.12
hist(ref_16S$rRNA16S)

pw_16S_df <- as.data.frame(pw_16S)
pw_16S_df$Row <- rownames(pw_16S_df)
long_16S <- pivot_longer(pw_16S_df, cols = -Row, names_to = "Column", values_to = "Value") %>%
  filter(Row != Column) %>%
  mutate(Pair1 = pmin(Row, Column), Pair2 = pmax(Row, Column)) %>%
  distinct(Pair1, Pair2, .keep_all = TRUE) %>%
  dplyr::select(-Pair1, -Pair2)
long_16S_strains <- long_16S %>%
  filter(str_detect(Row, "Strain") | str_detect(Column, "Strain"))
range(long_16S_strains$Value) # 87.73 to 99.73
long_16S_strain2 <- long_16S %>%
  filter(str_detect(Row, "Strain_2") | str_detect(Column, "Strain_2"))
range(long_16S_strain2$Value) # 87.73 to 97.85
long_16S_rest <- long_16S %>%
  filter(!str_detect(Row, "Strain_2") & !str_detect(Column, "Strain_2"))
range(long_16S_rest$Value) # 92.76 to 100
long_16S_GTDB <- long_16S %>%
  filter(!str_detect(Row, "Strain") & !str_detect(Column, "Strain"))
range(long_16S_GTDB$Value) # 95.77 to 100



#### ___OTUs ####
# 1. hclust
set.seed(100)
identity_matrix <- as.matrix(pw_16S)
dissimilarity_matrix <- 1 - identity_matrix
summary(dissimilarity_matrix)
diag(dissimilarity_matrix) <- 0
dist_matrix <- as.dist(dissimilarity_matrix)
hclust_result <- hclust(dist_matrix, method = "complete") # Use complete!
hclust_result$height <- sort(hclust_result$height)
plot(hclust_result)
otu_16S <- data.frame("Cutoff" = c(seq(87, 100, 1))) %>%
  mutate(Metric = "16S") %>%
  mutate("Count" = c(
    length(unique(cutree(hclust_result, h = 0.13))),
    length(unique(cutree(hclust_result, h = 0.12))),
    length(unique(cutree(hclust_result, h = 0.11))),
    length(unique(cutree(hclust_result, h = 0.10))),
    length(unique(cutree(hclust_result, h = 0.09))),
    length(unique(cutree(hclust_result, h = 0.08))),
    length(unique(cutree(hclust_result, h = 0.07))),
    length(unique(cutree(hclust_result, h = 0.06))),
    length(unique(cutree(hclust_result, h = 0.05))),
    length(unique(cutree(hclust_result, h = 0.04))),
    length(unique(cutree(hclust_result, h = 0.03))),
    length(unique(cutree(hclust_result, h = 0.02))),
    length(unique(cutree(hclust_result, h = 0.01))),
    length(unique(cutree(hclust_result, h = 0.00)))))

# 2. Igraph, adjacency matrix
adj_87 <- (identity_matrix >= 0.87) * 1
adj_88 <- (identity_matrix >= 0.88) * 1
adj_97 <- (identity_matrix >= 0.97) * 1
adj_99 <- (identity_matrix >= 0.99) * 1
adj_100 <- (identity_matrix == 1.00) * 1
diag(adj_87) <- 0
diag(adj_88) <- 0
diag(adj_97) <- 0
diag(adj_99) <- 0
diag(adj_99) <- 0
graph_87 <- graph_from_adjacency_matrix(adj_87, mode = "undirected", diag = FALSE)
graph_88 <- graph_from_adjacency_matrix(adj_88, mode = "undirected", diag = FALSE)
graph_97 <- graph_from_adjacency_matrix(adj_97, mode = "undirected", diag = FALSE)
graph_99 <- graph_from_adjacency_matrix(adj_99, mode = "undirected", diag = FALSE)
graph_100 <- graph_from_adjacency_matrix(adj_100, mode = "undirected", diag = FALSE)
OTUs_87 <- components(graph_87)$no
OTUs_88 <- components(graph_88)$no
OTUs_97 <- components(graph_97)$no
OTUs_99 <- components(graph_99)$no
OTUs_100 <- components(graph_100)$no

# 3. Count unique taxa below cutoff. Add 1 to capture the rest.
otu_16S <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 2)) %>%
  set_names(c("Cutoff", "Count"))
for (i in 87:100) {
  df <- long_16S %>%
    filter(Value*100 < i) %>%
    dplyr::select(Row, Column) %>%
    unlist() %>%
    unique()
  otu_16S$Cutoff[i] <- i
  otu_16S$Count[i] <- length(df) + 1
}

# 4. DECIPHER clustering
library(DECIPHER)
dna <- readDNAStringSet("data/combined_16S.fasta")
clusters <- Clusterize(dna, 
                       cutoff = c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                                  0.08, 0.09, 0.10, 0.11, 0.12, 0.13), 
                       processors = 4,
                       method = "longest")
otu_16S_decipher <- data.frame("Cutoff" = c(seq(87, 100, 1))) %>%
  mutate("Count" = c(
    length(unique(clusters$cluster_0_13)),
    length(unique(clusters$cluster_0_12)),
    length(unique(clusters$cluster_0_11)),
    length(unique(clusters$cluster_0_1)),
    length(unique(clusters$cluster_0_09)),
    length(unique(clusters$cluster_0_08)),
    length(unique(clusters$cluster_0_07)),
    length(unique(clusters$cluster_0_06)),
    length(unique(clusters$cluster_0_05)),
    length(unique(clusters$cluster_0_04)),
    length(unique(clusters$cluster_0_03)),
    length(unique(clusters$cluster_0_02)),
    length(unique(clusters$cluster_0_01)),
    length(unique(clusters$cluster_0))))

# Now also do ANI
ani_181 <- read.delim("data/ANI_ani_181.txt") %>%
  column_to_rownames(var = "layers") %>%
  dplyr::select(all_of(rownames(.)))
sum(ani_181[ani_181 == 1])
diag <- diag(as.matrix(ani_181))
sum(diag < 1)
lt <- as.matrix(ani_181)[lower.tri(as.matrix(ani_181))]
range(lt)
set.seed(100)
identity_matrix <- as.matrix(ani_181)
dissimilarity_matrix <- 1 - identity_matrix
diag(dissimilarity_matrix) <- 0
dist_matrix <- as.dist(dissimilarity_matrix)
hclust_result <- hclust(dist_matrix, method = "complete")
hclust_result$height <- sort(hclust_result$height)
plot(hclust_result)
otu_ANI <- data.frame("Cutoff" = c(seq(79, 100, 1))) %>%
  mutate("Metric" = "ANI") %>%
  mutate("Count" = c(
    length(unique(cutree(hclust_result, h = 0.21))),
    length(unique(cutree(hclust_result, h = 0.20))),
    length(unique(cutree(hclust_result, h = 0.19))),
    length(unique(cutree(hclust_result, h = 0.18))),
    length(unique(cutree(hclust_result, h = 0.17))),
    length(unique(cutree(hclust_result, h = 0.16))),
    length(unique(cutree(hclust_result, h = 0.15))),
    length(unique(cutree(hclust_result, h = 0.14))),
    length(unique(cutree(hclust_result, h = 0.13))),
    length(unique(cutree(hclust_result, h = 0.12))),
    length(unique(cutree(hclust_result, h = 0.11))),
    length(unique(cutree(hclust_result, h = 0.10))),
    length(unique(cutree(hclust_result, h = 0.09))),
    length(unique(cutree(hclust_result, h = 0.08))),
    length(unique(cutree(hclust_result, h = 0.07))),
    length(unique(cutree(hclust_result, h = 0.06))),
    length(unique(cutree(hclust_result, h = 0.05))),
    length(unique(cutree(hclust_result, h = 0.04))),
    length(unique(cutree(hclust_result, h = 0.03))),
    length(unique(cutree(hclust_result, h = 0.02))),
    length(unique(cutree(hclust_result, h = 0.01))),
    length(unique(cutree(hclust_result, h = 0.00)))))
otu_comb <- rbind(otu_16S, otu_ANI)

# Now just EMP 515-806 part of the 16S
combined_16S <- microseq::readFasta("data/combined_16S.fasta")
EMP <- combined_16S %>%
  mutate(Sequence = substr(Sequence, start = 515, stop = 806))
microseq::writeFasta(EMP, "data/EMP_16S.fasta")
pw_16S_EMP <- read.table("data/pairwise_identity_16S_181_EMP.txt",
                     header = TRUE) %>%
  column_to_rownames(var = "ID")
names(pw_16S_EMP) <- gsub(".fna", "", names(pw_16S_EMP))
names(pw_16S_EMP) <- gsub("Brady_", "", names(pw_16S_EMP))
rownames(pw_16S_EMP) <- gsub(".fna", "", rownames(pw_16S_EMP))
rownames(pw_16S_EMP) <- gsub("Brady_", "", rownames(pw_16S_EMP))
range(pw_16S_EMP)
set.seed(100)
identity_matrix <- as.matrix(pw_16S_EMP)
dissimilarity_matrix <- 1 - identity_matrix
summary(dissimilarity_matrix)
diag(dissimilarity_matrix) <- 0
dist_matrix <- as.dist(dissimilarity_matrix)
hclust_result <- hclust(dist_matrix, method = "complete") # Use complete!
hclust_result$height <- sort(hclust_result$height)
plot(hclust_result)
otu_16S_EMP <- data.frame("Cutoff" = c(seq(84, 100, 1))) %>%
  mutate(Metric = "16S EMP") %>%
  mutate("Count" = c(
    length(unique(cutree(hclust_result, h = 0.16))),
    length(unique(cutree(hclust_result, h = 0.15))),
    length(unique(cutree(hclust_result, h = 0.14))),
    length(unique(cutree(hclust_result, h = 0.13))),
    length(unique(cutree(hclust_result, h = 0.12))),
    length(unique(cutree(hclust_result, h = 0.11))),
    length(unique(cutree(hclust_result, h = 0.10))),
    length(unique(cutree(hclust_result, h = 0.09))),
    length(unique(cutree(hclust_result, h = 0.08))),
    length(unique(cutree(hclust_result, h = 0.07))),
    length(unique(cutree(hclust_result, h = 0.06))),
    length(unique(cutree(hclust_result, h = 0.05))),
    length(unique(cutree(hclust_result, h = 0.04))),
    length(unique(cutree(hclust_result, h = 0.03))),
    length(unique(cutree(hclust_result, h = 0.02))),
    length(unique(cutree(hclust_result, h = 0.01))),
    length(unique(cutree(hclust_result, h = 0.00)))))

otu_comb <- rbind(otu_16S_EMP, otu_16S, otu_ANI) %>%
  mutate(Metric = factor(Metric, levels = c("16S EMP", "16S", "ANI")))

# Figure S5
# Need to plot % identity on X and # comparisons on y
# Draft 1
label <- data.frame(x = c(97, 99, 100),
                    y = c(2400, 2400, 2200),
                    label = c("33 OTUs", "69 OTUs", "149 OTUs"))
figS5 <- ggplot(long_16S, aes(x = Value*100)) +
  geom_histogram(breaks = c(87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
                            97, 98, 99, 100, 101),
                 closed = "left",
                 fill = "grey80", colour = "black") +
  stat_bin(aes(y= ..count.., label = ..count..), geom = "text", vjust = -0.5,
           breaks = c(87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
                      97, 98, 99, 100, 101),
           closed = "left") +
  geom_vline(xintercept = 97, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 99, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_label(data = label, aes(x, y, label = label), size = 3) +
  labs(x = "% 16S rRNA gene identity",
       y = "# pairwise comparisons") +
  scale_x_continuous(breaks = c(87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
                                97, 98, 99, 100, 101)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
figS5

# Draft 2
facet_names <- c("16S EMP" = "515-806 16S rRNA gene (n OTUs)",
                 "16S" = "Full length 16S rRNA gene (n OTUs)",
                 "ANI" = "ANI (n genomes)")
figS5 <- ggplot(otu_comb, aes(Cutoff, Count)) +
  geom_bar(stat = "identity") +
  labs(x = "Similarity cutoff (%)",
       y = NULL) +
  facet_wrap(~ Metric, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 8))
figS5
pdf("FinalFigs/FigureS5.pdf", width = 7, height = 5)
figS5
dev.off()
png("FinalFigs/FigureS5.png", width = 7, height = 5, units = "in", res = 300)
figS5
dev.off()

# Now check Angela's Australian 16S data
ang <- read.delim("~/Desktop/Fierer/Strains/australia_taxa_table.txt") %>%
  column_to_rownames(var = "X.ASV_ID") %>%
  mutate(Sum = rowSums(.[, 1:ncol(.) - 1])) %>%
  arrange(desc(Sum))
ang_prev <- read.delim("~/Desktop/Fierer/Strains/australia_taxa_table.txt") %>%
  column_to_rownames(var = "X.ASV_ID") %>%
  mutate(Prev = rowSums(.[, sapply(., is.numeric)] > 0)) %>%
  arrange(desc(Prev))
ang_brady <- ang %>%
  filter(grepl("Bradyrhizobium", taxonomy))

# Now check LUCAS soil (Labouryie et al. 2023)
lucas_taxonomy <- read.csv("~/Desktop/Fierer/Strains/LUCAS_Taxonomy.csv")
lucas_otu <- read.csv("~/Desktop/Fierer/Strains/LUCAS_OTU_16S.csv")
names(lucas_otu)
head(lucas_otu$X)
lucas_abund <- read.csv("~/Desktop/Fierer/Strains/LUCAS_OTU_16S.csv") %>%
  mutate(Sum = rowSums(.[, 2:ncol(.)])) %>%
  mutate(Prev = rowSums(.[, sapply(., is.numeric)] > 0)) %>%
  arrange(desc(Sum)) %>%
  left_join(., lucas_taxonomy, by = c("X" = "zOTU"))
# 1 and 2 are Bradyrhizobiaceae, 3 and 4 are Bradyrhizobium
lucas_prev <- read.csv("~/Desktop/Fierer/Strains/LUCAS_OTU_16S.csv") %>%
  mutate(Sum = rowSums(.[, 2:ncol(.)])) %>%
  mutate(Prev = rowSums(.[, sapply(., is.numeric)] > 0)) %>%
  arrange(desc(Prev)) %>%
  left_join(., lucas_taxonomy, by = c("X" = "zOTU"))
# 2 of top 5 most prevalent zOTU are Brady


# Old
pw_16S <- read.table("data/pairwise_identity_matrix.txt",
                     header = TRUE) %>%
  column_to_rownames(var = "ID")
names(pw_16S) <- gsub(".fna", "", names(pw_16S))
names(pw_16S) <- gsub("Brady_", "", names(pw_16S))
rownames(pw_16S) <- gsub(".fna", "", rownames(pw_16S))
rownames(pw_16S) <- gsub("Brady_", "", rownames(pw_16S))

# Strains to Ref
ref_16S <- pw_16S %>%
  dplyr::select(all_of(names(strains_ani)), Reference) %>%
  filter(rownames(.) %in% c(names(strains_ani), "Reference")) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(names(strains_ani)), Reference) %>%
  mutate_all(function(x) x * 100) %>%
  filter(rownames(.) != "Reference") %>%
  rownames_to_column(var = "sampleID") %>%
  rename("rRNA16S" = "Reference") %>%
  dplyr::select(sampleID, rRNA16S)
sum(ref_ani$sampleID != ref_16S$sampleID)
range(ref_16S$rRNA16S) # Range is 95.61231 to 97.18402
hist(ref_16S$rRNA16S)

d_pw <- d_pw %>% 
  left_join(., ref_16S, by = "sampleID")
ggplot(d_pw, aes(AI, rRNA16S)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  labs(x = "drier <= Aridity index => wetter",
       y = "16S % ID to reference") +
  theme_bw()

# Strains
strains_16S <- pw_16S %>%
  dplyr::select(all_of(names(strains_ani))) %>%
  filter(rownames(.) %in% names(strains_ani)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(names(strains_ani))) %>%
  mutate_all(function(x) x * 100)
dist.16S <- as.matrix(strains_16S)
dist.16S[upper.tri(dist.16S, diag = TRUE)] <- NA
hist(dist.16S)
range(dist.16S, na.rm = TRUE) # Range is 96.59463 99.73805



#### __Test ####
# Mantel
mantel(dist.ani, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.env, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.geog, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.ani, dist.env, method = "spearman", permutations = 2000)
mantel(dist.ani, dist.geog, method = "spearman", permutations = 2000)

mantel(dist.16S, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.env, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.geog, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.16S, dist.env, method = "spearman", permutations = 2000)
mantel(dist.16S, dist.geog, method = "spearman", permutations = 2000)

mantel(dist.16S, dist.ani, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.ani, method = "spearman", permutations = 2000)

# DBRDA
dist.16S.d <- as.dist(dist.16S)
mod0 <- dbrda(dist.16S.d ~ 1, d_brady_env)  # Model with intercept only
mod1 <- dbrda(dist.16S.d ~ ., d_brady_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # None

dist.ani.d <- as.dist(dist.ani)
mod0 <- dbrda(dist.ani.d ~ 1, d_brady_env)  # Model with intercept only
mod1 <- dbrda(dist.ani.d ~ ., d_brady_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # None

# Lengthen so we have all 1378 points
vec.ani <- as.vector(dist.ani)
vec.16S <- as.vector(dist.16S)
vec.geog <- as.vector(dist.geog)
vec.ai <- as.vector(dist.ai)
vec.env <- as.vector(dist.env)
mats <- data.frame(ANI = vec.ani,
                   rRNA16S = vec.16S,
                   Geog = vec.geog,
                   AI = vec.ai,
                   Env = vec.env) %>%
  drop_na()

# Test correlations, linear regressions, polynomial regressions (not valid though!)
cor.test(mats$ANI, mats$Geog, method = "pearson")
cor.test(mats$ANI, mats$AI, method = "pearson")
cor.test(mats$ANI, mats$Env, method = "pearson")
cor.test(mats$ANI, mats$Geog, method = "spearman")
cor.test(mats$ANI, mats$AI, method = "spearman")
cor.test(mats$ANI, mats$Env, method = "spearman")
summary(lm(ANI ~ Geog, data = mats))
summary(lm(ANI ~ AI, data = mats))
summary(lm(ANI ~ Env, data = mats))
summary(lm(ANI ~ poly(Geog, 2, raw = TRUE), data = mats))
summary(lm(ANI ~ poly(AI, 2, raw = TRUE), data = mats))
summary(lm(ANI ~ poly(Env, 2, raw = TRUE), data = mats))

cor.test(mats$rRNA16S, mats$Geog, method = "pearson")
cor.test(mats$rRNA16S, mats$AI, method = "pearson")
cor.test(mats$rRNA16S, mats$Env, method = "pearson")
cor.test(mats$rRNA16S, mats$Geog, method = "spearman")
cor.test(mats$rRNA16S, mats$AI, method = "spearman")
cor.test(mats$rRNA16S, mats$Env, method = "spearman")
summary(lm(rRNA16S ~ Geog, data = mats))
summary(lm(rRNA16S ~ AI, data = mats))
summary(lm(rRNA16S ~ Env, data = mats))
summary(lm(rRNA16S ~ poly(Geog, 2, raw = TRUE), data = mats))
summary(lm(rRNA16S ~ poly(AI, 2, raw = TRUE), data = mats))
summary(lm(rRNA16S ~ poly(Env, 2, raw = TRUE), data = mats))



#### __Plot ####
ggplot(mats, aes(Geog, ANI)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Geographic distance (m)",
       y = "ANI (%)") +
  theme_bw()
ggplot(mats, aes(AI, ANI)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Aridity distance (Euclidean)",
       y = "ANI (%)") +
  theme_bw()
ggplot(mats, aes(Env, ANI)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Environmental distance (Euclidean)",
       y = "ANI (%)") +
  theme_bw()

ggplot(mats, aes(Geog, rRNA16S)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Geographic distance (m)",
       y = "16S (% ID)") +
  theme_bw()
ggplot(mats, aes(AI, rRNA16S)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Aridity distance (Euclidean)",
       y = "16S (% ID)") +
  theme_bw()
ggplot(mats, aes(Env, rRNA16S)) +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  labs(x = "Environmental distance (Euclidean)",
       y = "16S (% ID)") +
  theme_bw()

# Multipanel Plot
facet_names <- c("AI" = "Aridity Index",
                 "Env" = "Env. (Euclidean)",
                 "Geog" = "Geog. (m)")
mats_long <- mats %>%
  pivot_longer(cols = c("Geog", "AI", "Env"))
pdf("InitialFigs/ANI_dist.pdf", width = 7, height = 4)
ggplot(mats_long, aes(value, ANI)) +
  geom_point(size = 1.5, alpha = 0.2, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  facet_wrap(~ name, ncol = 3, scales = "free_x", labeller = as_labeller(facet_names)) +
  labs(x = "Distance",
       y = "ANI (%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(c(1, 1), "lines"))
dev.off()

pdf("InitialFigs/16S_dist.pdf", width = 7, height = 4)
ggplot(mats_long, aes(value, rRNA16S)) +
  geom_point(size = 1.5, alpha = 0.2, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  facet_wrap(~ name, ncol = 3, scales = "free_x", labeller = as_labeller(facet_names)) +
  labs(x = "Distance",
       y = "16S (% ID)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(c(1, 1), "lines"))
dev.off()

facet_names <- c("AI" = "Aridity Index",
                 "Env" = "Env. (Euclidean)",
                 "Geog" = "Geog. (m)",
                 "ANI" = "ANI",
                 "rRNA16S" = "16S")
mats_long1 <- mats_long %>%
  mutate(Dataset = "ANI") %>%
  dplyr::select(-rRNA16S) %>%
  rename("y" = "ANI")
mats_long2 <- mats_long %>%
  mutate(Dataset = "rRNA16S") %>%
  dplyr::select(-ANI) %>%
  rename("y" = "rRNA16S")
mats_long_both <- rbind(mats_long1, mats_long2)
pdf("InitialFigs/Dist_ANI_16S.pdf", width = 7, height = 6)
ggplot(mats_long_both, aes(value, y)) +
  geom_point(size = 1.5, alpha = 0.2, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  facet_grid(Dataset ~ name, scales = "free_x", labeller = as_labeller(facet_names)) +
  labs(x = "Distance",
       y = "% ID") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(c(1, 1), "lines"))
dev.off()

# Also ANI vs 16S
mantel(dist.ani, dist.16S, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.16S, method = "spearman", permutations = 2000)
cor.test(mats$ANI, mats$rRNA16S, method = "pearson")
cor.test(mats$ANI, mats$rRNA16S, method = "spearman")
summary(lm(ANI ~ rRNA16S, data = mats))
pdf("InitialFigs/ANI_16S.pdf", width = 7, height = 5)
ggplot(mats, aes(rRNA16S, ANI)) +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.25, pch = 16) +
  geom_smooth(method = "lm", se = FALSE) +
  # xlim(95, 100) +
  # ylim(95, 100) +
  labs(x = "16S % ID",
       y = "ANI (%)") +
  theme_bw()
dev.off()

# PCoA with vectors
pcoa_16S <- cmdscale(dist.16S.d, k = nrow(d_brady) - 1, eig = T)
set.seed(100)
ef_16S <- envfit(pcoa_16S, d_brady_env, permutations = 999, na.rm = TRUE)
ef_16S # None
ordiplot(pcoa_16S)
plot(ef_16S, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_16S)
vec.df_16S <- as.data.frame(ef_16S$vectors$arrows*sqrt(ef_16S$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_16S,
         Dim2 = Dim2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  #filter(ef_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Conductivity", "NO3", "Org. C", "Temp", "Precip", "AI"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
d_brady$Axis01 <- vegan::scores(pcoa_16S)[,1]
d_brady$Axis02 <- vegan::scores(pcoa_16S)[,2]
pdf("InitialFigs/PCoA_16S.pdf", width = 7, height = 5)
ggplot(d_brady, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_16S,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "red", fontface = "bold") +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S,
       colour = "Site") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

pcoa_ani <- cmdscale(dist.ani.d, k = nrow(d_brady) - 1, eig = T)
set.seed(100)
ef_ani <- envfit(pcoa_ani, d_brady_env, permutations = 999, na.rm = TRUE)
ef_ani # None
ordiplot(pcoa_ani)
plot(ef_ani, cex = 0.5)
multiplier_ani <- ordiArrowMul(ef_ani)
vec.df_ani <- as.data.frame(ef_ani$vectors$arrows*sqrt(ef_ani$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_ani,
         Dim2 = Dim2 * multiplier_ani) %>%
  mutate(variables = rownames(.)) %>%
  #filter(ef_ani$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Conductivity", "NO3", "Org. C", "Temp", "Precip", "AI"))
pcoaA1_ani <- paste("PC1: ", round((eigenvals(pcoa_ani)/sum(eigenvals(pcoa_ani)))[1]*100, 1), "%")
pcoaA2_ani <- paste("PC2: ", round((eigenvals(pcoa_ani)/sum(eigenvals(pcoa_ani)))[2]*100, 1), "%")
d_brady$Axis01 <- vegan::scores(pcoa_ani)[,1]
d_brady$Axis02 <- vegan::scores(pcoa_ani)[,2]
pdf("InitialFigs/PCoA_ani.pdf", width = 7, height = 5)
ggplot(d_brady, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df_ani,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_ani,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "red", fontface = "bold") +
  labs(x = pcoaA1_ani, 
       y = pcoaA2_ani,
       colour = "Site") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



#### _Genes ####
# Import DRAM annotation for the strain genomes
# Get number of KO per strain genome
# Compute shared and unique KOs
# Make complex upset plot
a <- read.delim("~/Desktop/Fierer/Strains/Australia/Brady_annotations.tsv") %>%
  mutate(fasta = gsub("-Reads_1-Strain_1.Genome", "", fasta)) %>%
  mutate(fasta = gsub("Brady_", "", fasta)) %>%
  dplyr::select(fasta, ko_id, kegg_hit) %>%
  rename(sampleID = fasta) %>%
  filter(ko_id != "")
length(unique(a$ko_id))
ko_func <- a %>%
  dplyr::select(ko_id, kegg_hit) %>%
  group_by(ko_id) %>%
  slice_head(n = 1) %>%
  ungroup()
ko_count <- a %>%
  group_by(sampleID) %>%
  summarize(count = n()) %>%
  left_join(., d_AI, by = "sampleID")
hist(ko_count$count)
min(ko_count$count) # 3423, including copies!
max(ko_count$count) # 3479, including copies!
ggplot(ko_count, aes(AI, count)) +
  geom_point(size = 2, alpha = 0.5, pch = 16) +
  geom_smooth() +
  theme_bw()

brady_mat <- a %>%
  select(sampleID, ko_id) %>%
  rename(name = ko_id) %>%
  mutate(value = 1) %>%
  group_by(sampleID, name) %>%
  slice_head(n = 1) %>%
  pivot_wider() %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "sampleID") %>%
  t() %>%
  as.data.frame()
ko_count_uniq <- data.frame("count" = colSums(brady_mat))
# 2187 total unique KOs
# 1544 in all, 309 all but 1, 17 1
# 17 with 1, that's 13 strains, 13 with 1 unique and 4 with 2 unique
2187-1544-309-17 # 317 in all other combinations
# Can also just plot this info
int_info <- data.frame("Strains" = c("53", "52", "2-51", "1"),
                       "KOs" = c(1544, 309, 317, 17))
ggplot(int_info, aes(Strains, KOs)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of strains",
       y = "Number of KOs") +
  theme_bw()
# Or even get info for all combinations 1-53
int_info_all <- data.frame(Strains = seq(1:53)) %>%
  mutate(Strains = factor(Strains,
                          levels = seq(1:53))) %>%
  mutate(KOs = "NA")
for (i in 1:53) {
  k <- brady_mat %>%
    mutate(sum = rowSums(.)) %>%
    filter(sum == i) %>%
    select(-sum)
  int_info_all$KOs[i] <- nrow(k)
}
int_info_all$KOs <- as.integer(int_info_all$KOs)
pdf("InitialFigs/Genes_SharedKOs.pdf", width = 7, height = 4)
ggplot(int_info_all, aes(Strains, KOs)) +
  geom_bar(stat = "identity") +
  geom_text(aes(Strains, KOs + 20, label = KOs), size = 2) +
  labs(x = "Number of strains",
       y = "Number of shared KOs") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6))
dev.off()

# Get KOs in all 53, in 52/53, in 1
k53 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 53) %>%
  select(-sum)
k52 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 52) %>%
  select(-sum)
k1 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 1) %>%
  select(-sum)
kotab <- rbind(k1, 
               #k52, 
               k53)
# Upset plot (order by AI, remember high AI is wet)
d_AI_ordered <- d_AI %>%
  arrange(AI) %>%
  filter(sampleID %in% names(kotab))
kotab_reorder <- kotab %>%
  dplyr::select(all_of(d_AI_ordered$sampleID))
up1 <- upset(kotab,
             intersect = names(kotab_reorder),
             sort_sets = FALSE,
             set_sizes = FALSE,
             height_ratio = 4,
             name = NULL,
             base_annotations = list('KOs' = intersection_size(counts=T)),
             themes = upset_modify_themes(
               list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40),
                                                                              face = "italic"),
                                                   axis.title.x = element_text(size = 16, hjust = 0.5),
                                                   plot.margin = margin(l = 30)))))
pdf("InitialFigs/Genes_Upset.pdf", width = 6, height = 10)
up1
dev.off()

# Check actual functions in the 17 KOs in just 1 strain
k1_func <- a %>%
  subset(ko_id %in% rownames(k1))

# Need to do some logistic regression
# Which KOs are significantly more likely to present based on AI?
# Can't do if only 1's
brady_mat2 <- brady_mat %>%
  mutate(tot = rowSums(.)) %>%
  filter(tot != 53) %>%
  dplyr::select(-tot)
brady_mat_reorder <- brady_mat2 %>%
  dplyr::select(all_of(d_AI_ordered$sampleID))
sum(names(brady_mat_reorder) != d_AI_ordered$sampleID)
# Loop through each KO, perform logistic regression with AI
brady_mat_reorder_t <- as.data.frame(t(brady_mat_reorder))
sigKOs <- data.frame(KO = names(brady_mat_reorder_t),
                     p = NA,
                     PseudoR2 = NA)
for (i in 1:ncol(brady_mat_reorder_t)) {
  m <- glm(brady_mat_reorder_t[,i] ~ d_AI_ordered$AI, family = "binomial")
  r <- nagelkerke(m, null = NULL, restrictNobs = FALSE)
  sigKOs$p[i] <- coef(summary(m))[2,4]
  sigKOs$PseudoR2[i] <- r$Pseudo.R.squared.for.model.vs.null[3]
}
sigKOs$Pfdr <- p.adjust(sigKOs$p, method = "fdr")
sigKOs$Pbon <- p.adjust(sigKOs$p, method = "bonferroni")
brady_mat3 <- brady_mat %>%
  mutate(tot = rowSums(.)) %>%
  filter(rownames(.) %in% sigKOs$KO) %>%
  rownames_to_column(var = "KO") %>%
  dplyr::select(KO, tot)
sigKOs <- sigKOs %>%
  filter(p < 0.05) %>%
  dplyr::select(-Pfdr, -Pbon) %>%
  left_join(., ko_func, by = c("KO" = "ko_id")) %>%
  left_join(., brady_mat3, by = "KO")
sigKOs_good <- sigKOs %>%
  filter(tot >= 10) %>%
  filter(tot <= 43)
sigKO_df <- brady_mat_reorder_t %>%
  dplyr::select(all_of(sigKOs$KO)) %>%
  mutate(AI = d_AI_ordered$AI)

pdf("InitialFigs/LogRegK00982.pdf", width = 7, height = 5)
ggplot(sigKO_df, aes(x = AI, y = K00982)) + 
  geom_point() +
  stat_smooth(method = "glm", color = "blue", se = F, 
              method.args = list(family = binomial)) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Predicted probability") +
  ggtitle("K00982: adenylyltransferase [EC:2.7.7.42]") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()
pdf("InitialFigs/LogRegK01577.pdf", width = 7, height = 5)
ggplot(sigKO_df, aes(x = AI, y = K01577)) + 
  geom_point() +
  stat_smooth(method = "glm", color = "blue", se = F, 
              method.args = list(family = binomial)) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Predicted probability") +
  ggtitle("K01577: oxalyl-CoA decarboxylase [EC:4.1.1.8]") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()
pdf("InitialFigs/LogRegK02068.pdf", width = 7, height = 5)
ggplot(sigKO_df, aes(x = AI, y = K02068)) + 
  geom_point() +
  stat_smooth(method = "glm", color = "blue", se = F, 
              method.args = list(family = binomial)) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Predicted probability") +
  ggtitle("K02068: UDP-glucose/iron transport system ATP-binding protein") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()
pdf("InitialFigs/LogRegK16074.pdf", width = 7, height = 5)
ggplot(sigKO_df, aes(x = AI, y = K16074)) + 
  geom_point() +
  stat_smooth(method = "glm", color = "blue", se = F, 
              method.args = list(family = binomial)) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Predicted probability") +
  ggtitle("K16074: zinc transporter") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()
pdf("InitialFigs/LogRegK03046.pdf", width = 7, height = 5)
ggplot(sigKO_df, aes(x = AI, y = K03046)) + 
  geom_point() +
  stat_smooth(method = "glm", color = "blue", se = F, 
              method.args = list(family = binomial)) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Predicted probability") +
  ggtitle("K03046: DNA-directed RNA polymerase subunit beta' [EC:2.7.7.6]") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()


#### _Gene Clusters ####
# From anvio pangenomics analysis
# Pangenome summary file (anvi-summarize)
# Note: 7,842 gene clusters with 399,027 genes
gc_all <- read.delim("~/Desktop/Fierer/Strains/Australia/BRADY-SUMMARY/Bradyrhizobium_Pan_gene_clusters_summary.txt")
names(gc_all)

# Check if all genes in the gc have the same function
gc_check_fun <- gc_all %>%
  group_by(gene_cluster_id, COG20_FUNCTION_ACC) %>%
  summarise(n = n())
nrow(gc_check_fun) # 8205 is more than 7842, so no.



#### __Patterns ####
# Collapse to 1 row per gene cluster, with 1 COG category (most frequent) each
# Wrangle COG Categories. Note 2392 unknown function
gc <- read.delim("~/Desktop/Fierer/Strains/Australia/BRADY-SUMMARY/Bradyrhizobium_Pan_gene_clusters_summary.txt") %>%
  arrange(gene_cluster_id, COG20_CATEGORY) %>%
  group_by(gene_cluster_id) %>%
  add_count(COG20_CATEGORY, sort = TRUE) %>% # Count occurrences of 'COG20_CATEGORY' within each group
  slice_max(n, n = 1, with_ties = FALSE) %>% # Select the most frequent
  #slice_head(n = 1) %>% # Old version, just getting first COG
  mutate(COG20_CATEGORY = ifelse(COG20_CATEGORY == "", "Unknown", COG20_CATEGORY)) %>%
  separate(COG20_CATEGORY, remove = F, sep = "\\|",
           into = c("COG1", "COG2", "COG3", "COG4", "COG5", "COG6", "COG7", "COG8"))
dup_fun <- gc %>%
  rowwise() %>%
  transmute(all_equal = n_distinct(c_across(c("COG1", "COG2", "COG3", "COG4", 
                                              "COG5", "COG6", "COG7", "COG8")), 
                                   na.rm = TRUE) == 1) %>%
  ungroup()
gc$all_equal <- dup_fun$all_equal
gc <- gc %>%
  mutate(COG20_Uniq = ifelse(all_equal == TRUE, COG1, "Multiple"))

# Prep data and plot % of gene clusters by COG
cog_df_all <- as.data.frame(table(gc$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate(All = Freq/nrow(gc)*100) %>%
  dplyr::select(COG, All)
gc_53 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 53)
cog_df_53 <- as.data.frame(table(gc_53$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_53" = Freq/nrow(gc_53)*100) %>%
  dplyr::select(COG, Genomes_53)
gc_252 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 52 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_252 <- as.data.frame(table(gc_252$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_2_to_52" = Freq/nrow(gc_252)*100) %>%
  dplyr::select(COG, Genomes_2_to_52)
gc_1 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 1)
cog_df_1 <- as.data.frame(table(gc_1$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_1" = Freq/nrow(gc_1)*100) %>%
  dplyr::select(COG, Genomes_1)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_53, by = c("COG")) %>%
  left_join(., cog_df_252, by = c("COG")) %>%
  left_join(., cog_df_1, by = c("COG")) %>%
  #replace(is.na(.), 0) %>%
  column_to_rownames(var = "COG")
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
pheatmap(cog_comb,
         legend = T,
         legend_breaks = c(0.0255037, 20, 40, 60, 80, 88.99676),
         legend_labels = c("0", "20", "40", "60", "80", ""),
         main = "            % Gene Clusters by COG",
         cluster_rows = F,
         cluster_cols = F,
         labels_col = c("All (n = 7842)", 
                        "In 53 genomes (n = 6560)",
                        "In 2 to 52 genomes (n = 973)",
                        "In 1 genome (n = 309)"),
         angle_col = 315,
         display_numbers = T,
         number_color = "black",
         fontsize_number = 10,
         border_color = "white",
         filename = "InitialFigs/GeneClusters_COGs.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Prep data and plot homogeneity by COG
cog_df_all <- gc %>%
  group_by(COG20_Uniq) %>%
  summarise(All_FH = mean(functional_homogeneity_index),
            All_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
cog_df_53 <- gc_53 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes53_FH = mean(functional_homogeneity_index),
            Genomes53_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_2752 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 52 & num_genomes_gene_cluster_has_hits >= 27)
cog_df_2752 <- gc_2752 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes2752_FH = mean(functional_homogeneity_index),
            Genomes2752_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_226 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 26 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_226 <- gc_226 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes226_FH = mean(functional_homogeneity_index),
            Genomes226_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_53, by = c("COG")) %>%
  left_join(., cog_df_2752, by = c("COG")) %>%
  left_join(., cog_df_226, by = c("COG")) %>%
  column_to_rownames(var = "COG") %>%
  dplyr::select(All_FH, Genomes53_FH, Genomes2752_FH, Genomes226_FH,
                All_GH, Genomes53_GH, Genomes2752_GH, Genomes226_GH)
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
ann_cols <- data.frame(row.names = colnames(cog_comb),
                       "Index" = c(rep("Functional", 4),
                                   rep("Geometric", 4)))
ann_colors <- list(Index = c(Functional = "#F8766D",
                             Geometric = "#619CFF"))
pheatmap(cog_comb,
         legend = T,
         main = "               Gene Clusters Homogeneity by COG",
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = 4,
         labels_col = c("All (n = 7842)",
                        "In 53 genomes (n = 6560)",
                        "In 27 to 52 genomes (n = 508)",
                        "In 2 to 26 genomes (n = 465)",
                        "All (n = 7842)",
                        "In 53 genomes (n = 6560)",
                        "In 27 to 52 genomes (n = 508)",
                        "In 2 to 26 genomes (n = 465)"),
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         display_numbers = T,
         number_format = "%.3f",
         number_color = "black",
         fontsize_number = 6,
         border_color = "white",
         filename = "InitialFigs/GeneClusters_Homogeneity.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### __Heterogeneity ####
# Selected gene clusters in all 53 genomes but most heterogeneous
# Functional homogeneity and geometric homogeneity, below 0.75
hetF_names <- read.delim("data/HetF_GC_names.txt", header = F) %>%
  rename(gene_cluster_id = V1) %>%
  mutate(gene_cluster_id = gsub("BRADY_", "", gene_cluster_id)) %>%
  mutate(gene_cluster_id = gsub(".fa", "", gene_cluster_id))
hetG_names <- read.delim("data/HetG_GC_names.txt", header = F) %>%
  rename(gene_cluster_id = V1) %>%
  mutate(gene_cluster_id = gsub("BRADY_", "", gene_cluster_id)) %>%
  mutate(gene_cluster_id = gsub(".fa", "", gene_cluster_id))

# Functional Homogeneity
gcF <- gc_all %>%
  filter(gene_cluster_id %in% hetF_names$gene_cluster_id)
length(unique(gcF$gene_cluster_id))
# Get one function per gene cluster
gcF <- gc_all %>%
  filter(gene_cluster_id %in% hetF_names$gene_cluster_id) %>%
  group_by(gene_cluster_id, COG20_FUNCTION) %>%
  summarise(n = n(), FH = mean(functional_homogeneity_index)) %>%
  slice_tail(n = 1) %>%
  arrange(FH)
  
# Geometric Homogeneity
gcG <- gc_all %>%
  filter(gene_cluster_id %in% hetG_names$gene_cluster_id)
length(unique(gcG$gene_cluster_id))
# Get one function per gene cluster
gcG <- gc_all %>%
  filter(gene_cluster_id %in% hetG_names$gene_cluster_id) %>%
  group_by(gene_cluster_id, COG20_FUNCTION) %>%
  summarise(n = n(), GH = mean(geometric_homogeneity_index)) %>%
  slice_tail(n = 1) %>%
  arrange(GH)
# Paste those into Excel spreadsheet



#### _Phylogeny ####
# Need to know if the strains are the same species or not
# Full
gtdb_brady <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus == "Bradyrhizobium")
# 869 genomes
# Make sure the reference used in StrainFinder is there
# GCA_016616885.1
"GCA_016616885.1" %in% gtdb_brady$ncbi_genbank_assembly_accession # TRUE
# Save these 869 accessions as a list
# write.table(gtdb_brady$ncbi_genbank_assembly_accession,
#             "data/brady869.txt", sep = "\t", row.names = F, col.names = F)
# Don't forget to run dos2unix

# dRep
gtdb_brady_dRep <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") %>%
  filter(gtdb_representative == "t") %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";") %>%
  filter(Genus == "Bradyrhizobium")
# 278 genomes
# Make sure the reference used in StrainFinder is there
# GCA_016616885.1
"GCA_016616885.1" %in% gtdb_brady_dRep$ncbi_genbank_assembly_accession # TRUE
# Save these 278 accessions as a list
# write.table(gtdb_brady_dRep$ncbi_genbank_assembly_accession,
#             "data/brady278.txt", sep = "\t", row.names = F, col.names = F)
# Don't forget to run dos2unix



#### __16S ####
# Extracted 16S with barrnap
counts_16S <- read.table("data/16S_counts.txt") %>%
  dplyr::select(V1, V2) %>%
  set_names(c("Genome", "Count16S"))
table(counts_16S$Count16S)
# 36 with 0, 276 with 1, 14 with 2, 1 with 3, 5 with 4
# Of the 20 with multiple, 8 are identical, 12 are different
# Just use the 276 with 1 16S

# Best tree from RAxML with bootstrap support
tree_16S <- ape::read.tree("data/RAxML_bipartitions.16S")
ggtree(tree_16S)
# Okay all of the Brady are in the same branch
# There are also 2 GTDB completely rescaling the tree - prune those
# Check those, maybe they didn't have full length 16S...
# GCA_023435725.1.fna
# GCA_028290305.1.fna
# BLAST says these are both uncultured bacterium clones
tree_16S_pruned <- drop.tip(tree_16S,
                            tip = c("GCA_023435725.1.fna",
                                    "GCA_028290305.1.fna"))
tree_16S_pruned$tip.label <- gsub(".fna", "", tree_16S_pruned$tip.label)
tree_16S_pruned$tip.label <- gsub("Brady_", "", tree_16S_pruned$tip.label)
tree_16S_pruned$tip.label <- gsub("_", "", tree_16S_pruned$tip.label)
tree_16S_pruned$tip.label
tipcols <- c(rep("black", 7),
             rep("red", 53),
             rep("black", 181),
             "blue",
             rep("black", 31),
             "purple")
pdf("InitialFigs/BradyTree_16S.pdf", width = 7, height = 9)
ggtree(tree_16S_pruned, linewidth = 0.1) +
  geom_tiplab(size = 1, vjust = 0.5,
              color = tipcols) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.75)
dev.off()
#write.tree(tree_16S_pruned, "data/tree_16S_pruned.nwk")
# Note the closest branch is Bradyrhizobium vignae genome assembly ASM411442v1
# That one was isolated from Arachis nodules



#### __bac120 ####
# Best tree without bootstrap values (bootstrap takes awhile!)
# Remember, reference was Bradyrhizobium diazoefficiens
# But closest 16S was Bradyrhizobium vignae
# In bac120, equally close are the reference and Bradyrhizobium sp. WSM3983
tree_bac120 <- read.tree("data/RAxML_bipartitions.bac120")
tree_bac120$tip.label <- gsub(".fna", "", tree_bac120$tip.label)
tree_bac120$tip.label <- gsub("Brady_", "", tree_bac120$tip.label)
tree_bac120$tip.label <- gsub("_", "", tree_bac120$tip.label)
tree_bac120$tip.label
tipcols <- c(rep("black", 138),
             rep("red", 53),
             "black",
             "blue",
             rep("black", 138),
             "purple")
pdf("InitialFigs/BradyTree_Bac120.pdf", width = 7, height = 9)
ggtree(tree_bac120, linewidth = 0.1) +
  geom_tiplab(size = 0.55, vjust = 0.5,
              color = tipcols) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.75)
dev.off()

# Brady only
tree_data <- d_brady %>%
  dplyr::select(sampleID, conductivity, nitrate_nitrogen, organic_carbon, bio1,
                bio12, AI) %>%
  column_to_rownames(var = "sampleID") %>%
  scale() %>%
  as.data.frame()
tree_bac120$tip.label
tree_bac120_pruned <- drop.tip(tree_bac120,
                               tip = c(tree_bac120$tip.label[1:138],
                                       tree_bac120$tip.label[192:332]))
p <- ggtree(tree_bac120_pruned, linewidth = 0.2) +
  geom_tiplab(size = 2, vjust = 0.5, align = TRUE) +
  geom_treescale(y = -2) +
  geom_nodelab(aes(label = label, fill = "white"), size = 2, hjust = 0) +
  new_scale_fill()
pdf("InitialFigs/BradyTree_Strains_Heatmap.pdf", width = 8, height = 6)
gheatmap(p, tree_data, offset = 0.01, width = 0.8, font.size = 3, 
         colnames_angle = -45, hjust = 0) +
  scale_fill_viridis_c(option = "B", name = "z-score") +
  scale_y_continuous(expand = c(0.15, 0)) +
  theme(legend.box.margin = margin(0,0,0,-20, unit = "pt"),
        plot.margin = margin(-55,0,-15,0, unit = "pt"))
dev.off()



#### 6. Pop. Genomics ####
# Assess strain population genomics from .vcf files
# Input_POGENOM and POGENOM
# Ran with fasta and with gff for additional gene by gene metrics
allele_freqs <- read.delim("data/pogenom_fasta_results/Brady.allele-freqs.txt")
pi_locus <- read.delim("data/pogenom_fasta_results/Brady.intradiv-per-locus.txt")
pi <- read.delim("data/pogenom_fasta_results/Brady.intradiv.txt") %>%
  filter(Sample != "All_samples_combined") %>%
  arrange(Sample) %>%
  select(Sample, Intra_pi) %>%
  rename(sampleID = Sample)



#### _Pairwise FST ####
fst <- read.delim("data/pogenom_fasta_results/Brady.fst.txt") %>%
  mutate(X = as.character(X)) %>%
  arrange(X) %>%
  column_to_rownames(var = "X") %>%
  set_names(gsub("X", "", names(.))) %>%
  dplyr::select(all_of(rownames(.)))
dist.fst <- as.matrix(fst)
dist.fst[upper.tri(dist.fst, diag = TRUE)] <- NA
hist(dist.fst)
range(fst, na.rm = T) # Range is -0.0919 to 0.3732

# Since fst only has 52 samples, need to reselect the other distance matrices
# Also add in pi!
d_brady <- d %>%
  filter(sampleID %in% rownames(fst)) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID) %>%
  left_join(., pi, by = "sampleID")
d_brady_env <- d_brady %>%
  select(sampleID, all_of(env_vars)) %>%
  select(where(~ all(!is.na(.)))) %>%
  select(-sampleID, -latitude, -longitude)
dist.16S <- dist.16S[rownames(dist.fst), rownames(dist.fst)]
dist.ai <- dist.ai[rownames(dist.fst), rownames(dist.fst)]
dist.ani <- dist.ani[rownames(dist.fst), rownames(dist.fst)]
dist.env <- dist.env[rownames(dist.fst), rownames(dist.fst)]
dist.geog <- dist.geog[rownames(dist.fst), rownames(dist.fst)]


# Map with 52 samples
coords <- d_brady %>%
  select(sampleID, latitude, longitude, AI)
coords_trans <- st_as_sf(coords, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
ozmap()
sf_oz <- ozmap("states")
pdf("InitialFigs/Mapn52.pdf", width = 7, height = 5)
ggplot(data = sf_oz) + 
  geom_sf(fill = "grey80", color = "white") +
  geom_sf(data = coords_trans,
          aes(color = AI)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal()
dev.off()



#### __Test ####
mantel(dist.fst, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.env, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.geog, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.fst, dist.env, method = "spearman", permutations = 2000)
mantel(dist.fst, dist.geog, method = "spearman", permutations = 2000)

dist.fst.d <- as.dist(dist.fst)
set.seed(1100)
mod0 <- dbrda(dist.fst.d ~ 1, d_brady_env)  # Model with intercept only
mod1 <- dbrda(dist.fst.d ~ ., d_brady_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # AI and bio1



#### __Plot ####
vec.ani <- as.vector(dist.ani)
vec.16S <- as.vector(dist.16S)
vec.fst <- as.vector(dist.fst)
vec.geog <- as.vector(dist.geog)
vec.ai <- as.vector(dist.ai)
vec.env <- as.vector(dist.env)
mats <- data.frame(ANI = vec.ani,
                   rRNA16S = vec.16S,
                   Fst = vec.fst,
                   Geog = vec.geog,
                   AI = vec.ai,
                   Env = vec.env) %>%
  drop_na()
mats_long <- mats %>%
  pivot_longer(cols = c("Geog", "AI", "Env"))
facet_names <- c("AI" = "'Aridity Index'",
                 "Env" = "'Env. (Euclidean)'",
                 "Geog" = "'Geog. (m)'",
                 "ANI" = "'ANI (% ID)'",
                 "rRNA16S" = "'16S (% ID)'",
                 "Fst" = "F[ST]")
mats_long1 <- mats_long %>%
  mutate(Dataset = "ANI") %>%
  dplyr::select(-rRNA16S, -Fst) %>%
  rename("y" = "ANI")
mats_long2 <- mats_long %>%
  mutate(Dataset = "rRNA16S") %>%
  dplyr::select(-ANI, -Fst) %>%
  rename("y" = "rRNA16S")
mats_long3 <- mats_long %>%
  mutate(Dataset = "Fst") %>%
  dplyr::select(-ANI, -rRNA16S) %>%
  rename("y" = "Fst")
mats_long_three <- rbind(mats_long1, mats_long2, mats_long3) %>%
  mutate(Dataset = factor(Dataset, levels = c("ANI", "rRNA16S", "Fst")))
# Pearson
mantel(dist.ani, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.env, method = "pearson", permutations = 2000)
mantel(dist.ani, dist.geog, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.env, method = "pearson", permutations = 2000)
mantel(dist.16S, dist.geog, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.ai, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.env, method = "pearson", permutations = 2000)
mantel(dist.fst, dist.geog, method = "pearson", permutations = 2000)
stats_label <- data.frame(Dataset = c("ANI", "ANI", "ANI", 
                                      "rRNA16S", "rRNA16S", "rRNA16S",
                                      "Fst", "Fst", "Fst"),
                          name = c("AI", "Env", "Geog", "AI", "Env", "Geog",
                                   "AI", "Env", "Geog"),
                          x = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                          y = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                          label = c("r = -0.17\np = 0.96",
                                    "r = -0.12\np = 0.88",
                                    "r = -0.25\np = 0.99",
                                    "r = -0.45\np = 1.00",
                                    "r = -0.41\np = 1.00",
                                    "r = -0.18\np = 0.97",
                                    "r = 0.02\np = 0.38",
                                    "r = -0.02\np = 0.48",
                                    "r = 0.10\np = 0.14"))  %>%
  mutate(Dataset = factor(Dataset, levels = c("ANI", "rRNA16S", "Fst")))
pdf("InitialFigs/Dist_ANI_16S_Fst_Pearson.pdf", width = 7, height = 8)
ggplot(mats_long_three, aes(value, y)) +
  geom_point(size = 1.5, alpha = 0.2, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  geom_text(data = stats_label, aes(x, y, label = label), 
            size = 2, hjust = 1.25, vjust = 1.4) +
  facet_grid(Dataset ~ name, scales = "free", 
             labeller = as_labeller(facet_names, default = label_parsed)) +
  labs(x = "Distance",
       y = NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(c(1, 1), "lines"))
dev.off()

# Spearman
mantel(dist.ani, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.ani, dist.env, method = "spearman", permutations = 2000)
mantel(dist.ani, dist.geog, method = "spearman", permutations = 2000)
mantel(dist.16S, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.16S, dist.env, method = "spearman", permutations = 2000)
mantel(dist.16S, dist.geog, method = "spearman", permutations = 2000)
mantel(dist.fst, dist.ai, method = "spearman", permutations = 2000)
mantel(dist.fst, dist.env, method = "spearman", permutations = 2000)
mantel(dist.fst, dist.geog, method = "spearman", permutations = 2000)
stats_label <- data.frame(Dataset = c("ANI", "ANI", "ANI", 
                                      "rRNA16S", "rRNA16S", "rRNA16S",
                                      "Fst", "Fst", "Fst"),
                          name = c("AI", "Env", "Geog", "AI", "Env", "Geog",
                                   "AI", "Env", "Geog"),
                          x = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                          y = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                          label = c("r = -0.14\np = 0.97",
                                    "r = -0.08\np = 0.85",
                                    "r = -0.23\np = 0.99",
                                    "r = -0.43\np = 1.00",
                                    "r = -0.40\np = 1.00",
                                    "r = -0.31\np = 1.00",
                                    "r = 0.05\np = 0.19",
                                    "r = -0.02\np = 0.48",
                                    "r = 0.14\np = 0.05"))  %>%
  mutate(Dataset = factor(Dataset, levels = c("ANI", "rRNA16S", "Fst")))
pdf("InitialFigs/Dist_ANI_16S_Fst.pdf", width = 7, height = 8)
ggplot(mats_long_three, aes(value, y)) +
  geom_point(size = 1.5, alpha = 0.2, pch = 16) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  geom_smooth(data = subset(mats_long_three, name == "Geog" & Dataset == "Fst"),
                            method = "lm", se = FALSE) +
  geom_text(data = stats_label, aes(x, y, label = label), 
            size = 2, hjust = 1.25, vjust = 1.4) +
  facet_grid(Dataset ~ name, scales = "free", 
             labeller = as_labeller(facet_names, default = label_parsed)) +
  labs(x = "Distance",
       y = NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(c(1, 1), "lines"))
dev.off()



#### _IntraPi ####
# Test for best model of nucleotide diversity
X <- d_brady_env %>%
  scale()
y <- scale(d_brady$Intra_pi)
Xy <- as.data.frame(cbind(X, y))
bestM <- bestglm(Xy,
                 family = gaussian,
                 IC = "AIC",
                 t = "default",
                 CVArgs = "default", 
                 qLevel = 0.99, 
                 TopModels = 10, 
                 method = "exhaustive", 
                 intercept = TRUE, 
                 weights = NULL, 
                 nvmax = "default", 
                 RequireFullEnumerationQ = FALSE)
bestM # Conductivity, Nitrate, Precip
bestM$BestModels

# Univariate
summary(lm(Intra_pi ~ conductivity, data = d_brady))
summary(lm(y ~ X[,1]))
summary(lm(Intra_pi ~ nitrate_nitrogen, data = d_brady))
summary(lm(y ~ X[,2]))

# Plot intraspecies nucleotide diversity by environment
ggplot(d_brady, aes(conductivity, Intra_pi)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  theme_bw()
pdf("InitialFigs/Pi_Nitrate.pdf", width = 7, height = 5)
ggplot(d_brady, aes(nitrate_nitrogen, Intra_pi)) +
  geom_point(size = 3, alpha = 0.75, pch = 16) +
  geom_smooth(method = "lm") +
  geom_text(x = Inf, y = Inf, 
            label = "R^2 == 0.11~~~p == 0.015", 
            check_overlap = T, hjust = 1.1, vjust = 1.4, parse = T) +
  labs(x = "Soil nitrate (mg/kg soil)",
       y = expression(Nucleotide~diversity~"("*pi*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12))
dev.off()
ggplot(d_brady, aes(organic_carbon, Intra_pi)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  theme_bw()
ggplot(d_brady, aes(bio1, Intra_pi)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  theme_bw()
ggplot(d_brady, aes(bio12, Intra_pi)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  theme_bw()
ggplot(d_brady, aes(AI, Intra_pi)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  theme_bw()



#### _Gene-wise ####
pNpS <- read.delim("data/pogenom_gff_results/Brady.pNpS-per-gene.txt")
pNpS_rows <- pNpS %>%
  na.omit() # 0
pNpS_cols <- pNpS %>%
  select(where(~ all(!is.na(.))))
pNpS_loci <- pNpS %>%
  filter(Num_loci > 0) %>%
  select(where(~ !all(is.na(.))))



#### 7. Strain Distributions ####
# Reran sylph 1) 53 strains, 2) 53 Strains + All Brady GTDB Genomes
# Assess distribution in 1) 53 sampkes, 2) 104 samples
# Thus, need to make 4 heatmaps
# Sort by AI
ref_ani <- read.delim("data/Brady_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                      "", Query)) %>%
  filter(grepl("Reference", Reference)) %>%
  filter(Query != "BradyReference.fna") %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID)) %>%
  dplyr::select(sampleID, ANI)
d_ai_sort <- d %>%
  arrange(AI)
d_brady <- d %>%
  filter(sampleID %in% ref_ani$sampleID) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID)
d_brady_ai_sort <- d_brady %>%
  arrange(AI)



#### _Strains ####
# 53 Samples
strains_sylph <- read.delim("data/sylph_profile_bradystrains_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrains/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID) %>%
  filter(SampleID %in% d_brady$sampleID)
length(unique(strains_sylph$Sample_file)) # 53
length(unique(strains_sylph$Genome_file)) # 53
# How many strains in their StrainFinder sample?
sum(strains_sylph$Same) # 53!
# How many strains in another sample?
sum(strains_sylph$Same == FALSE) # 4

# Plot (will mostly be 1 to 1 but also 4 other squares)
sylph_strains_53 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange(d_brady_ai_sort$sampleID) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(rownames(.)))
sum(rownames(sylph_strains_53) != colnames(sylph_strains_53))
sylph_strains_53 <- move_to_lower_triangle(sylph_strains_53)
sylph_strains_53 <- replace_lower_triangle_na(sylph_strains_53)
ann_cols <- data.frame(row.names = colnames(sylph_strains_53), 
                       AI = d_brady_ai_sort$AI)
pheatmap(sylph_strains_53,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Strains % Abundance, n = 53",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/Strains_Abund_53.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# 104 Samples
strains_sylph <- read.delim("data/sylph_profile_bradystrains_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrains/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID)
length(unique(strains_sylph$Sample_file)) # 75
length(unique(strains_sylph$Genome_file)) # 53
# How many strains in their StrainFinder sample?
sum(strains_sylph$Same) # 53!
# How many strains in another sample?
sum(strains_sylph$Same == FALSE) # 27

# Plot (will mostly be 1 to 1 but also 27 other squares)
missings <- d %>%
  filter(sampleID %notin% strains_sylph$SampleID)
sylph_strains_104 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange(match(GenomeID,  d_brady_ai_sort$sampleID)) %>%
  column_to_rownames(var = "GenomeID") %>%
  mutate(!!!set_names(rep(list(NA), length(missings$sampleID)), missings$sampleID)) %>% # Add missings
  dplyr::select(all_of(d_ai_sort$sampleID))
# sum(rownames(sylph_strains_104) != colnames(sylph_strains_104))
# sylph_strains_104 <- move_to_lower_triangle(sylph_strains_104)
# sylph_strains_104 <- replace_lower_triangle_na(sylph_strains_104)
ann_cols <- data.frame(row.names = colnames(sylph_strains_104), 
                       AI = d_ai_sort$AI) %>%
  mutate(StrainFinder = ifelse(rownames(.) %in% d_brady$sampleID,
                               "Yes",
                               "No"))
ann_colors <- list(StrainFinder = c(No = "#F8766D",
                                    Yes = "#619CFF"))
pheatmap(sylph_strains_104,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Strains % Abundance, n = 104",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 5,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/Strains_Abund_104.png",
         width = 9,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())




#### _Strains + GTDB ####
# Use the Strains + GTDB Sylph run but still only plot Strains
# 53 Samples
strains_sylph <- read.delim("data/sylph_profile_bradystrainsGTDB_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID) %>%
  filter(SampleID %in% d_brady$sampleID)
length(unique(strains_sylph$Sample_file)) # 53
length(unique(strains_sylph$Genome_file)) # 73, so 20 extra detected
# How many strains detected at all?
sum(d_brady$sampleID %in% strains_sylph$GenomeID) # 53
# How many strains in their StrainFinder sample?
sum(strains_sylph$Same) # 53
# How many strains in another sample?
sum(strains_sylph$Same == FALSE) # 62, including GTDB though
# How many samples with Reference
sum(strains_sylph$GenomeID == "Reference") # 6

# Plot (will mostly be 1 to 1 but also 4 other squares)
sylph_strains_53 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  filter(GenomeID %in% d_brady$sampleID) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange(d_brady_ai_sort$sampleID) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(rownames(.)))
sum(rownames(sylph_strains_53) != colnames(sylph_strains_53))
sylph_strains_53 <- move_to_lower_triangle(sylph_strains_53)
sylph_strains_53 <- replace_lower_triangle_na(sylph_strains_53)
ann_cols <- data.frame(row.names = colnames(sylph_strains_53), 
                       AI = d_brady_ai_sort$AI)
pheatmap(sylph_strains_53,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Strains % Abundance, n = 53",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/StrainsGTDB_Abund_53.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



# 104 Samples
strains_sylph <- read.delim("data/sylph_profile_bradystrainsGTDB_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID)
length(unique(strains_sylph$Sample_file)) # 92
length(unique(strains_sylph$Genome_file)) # 92
# How many strains detected at all?
sum(d_brady$sampleID %in% strains_sylph$GenomeID) # 53
# How many strains in their StrainFinder sample?
sum(strains_sylph$Same) # 53
# How many strains in another sample?
sum(strains_sylph$Same == FALSE) # 137 but includes other GTDB
# How many samples with Reference
sum(strains_sylph$GenomeID == "Reference") # 8

# Plot (will mostly be 1 to 1 but also other squares)
sylph_strains_104 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  filter(GenomeID %in% d_brady$sampleID)
missings <- d %>%
  filter(sampleID %notin% colnames(sylph_strains_104))
sylph_strains_104 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  filter(GenomeID %in% d_brady$sampleID) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange(match(GenomeID,  d_brady_ai_sort$sampleID)) %>%
  column_to_rownames(var = "GenomeID") %>%
  mutate(!!!set_names(rep(list(NA), length(missings$sampleID)), missings$sampleID)) %>% # Add missings
  dplyr::select(all_of(d_ai_sort$sampleID))
# sum(rownames(sylph_strains_104) != colnames(sylph_strains_104))
# sylph_strains_104 <- move_to_lower_triangle(sylph_strains_104)
# sylph_strains_104 <- replace_lower_triangle_na(sylph_strains_104)
ann_cols <- data.frame(row.names = colnames(sylph_strains_104), 
                       AI = d_ai_sort$AI) %>%
  mutate(StrainFinder = ifelse(rownames(.) %in% d_brady$sampleID,
                               "Yes",
                               "No"))
ann_colors <- list(StrainFinder = c(No = "#F8766D",
                                    Yes = "#619CFF"))
pheatmap(sylph_strains_104,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Strains % Abundance, n = 104",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 5,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/StrainsGTDB_Abund_104.png",
         width = 9,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _Strains w GTDB ####
# Use the Strains + GTDB Sylph run but plot all genomes hit at least once
# So instead of 53 genome IDs there will be more
# Matrix will not be symmetric

# 53
strains_sylph <- read.delim("data/sylph_profile_bradystrainsGTDB_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID) %>%
  filter(SampleID %in% d_brady$sampleID)

# Plot (will mostly be 1 to 1 but also other squares)
sylph_strains_53 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange("GenomeID") %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(d_brady_ai_sort$sampleID)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(d_brady_ai_sort$sampleID), Reference, everything()) %>%
  t() %>%
  as.data.frame()
ann_cols <- data.frame(row.names = colnames(sylph_strains_53), 
                       AI = d_brady_ai_sort$AI)
ann_rows <- data.frame(row.names = rownames(sylph_strains_53)) %>%
  mutate(Genome = ifelse(rownames(.) %in% d_brady$sampleID,
                         "Strain",
                         ifelse(rownames(.) == "Reference",
                                "Reference", "Other GTDB")))
ann_colors <- list(Genome = c("Strain" = "#00BFC4",
                              "Reference" = "#F8766D",
                              "Other GTDB" = "#C77CFF"))
pheatmap(sylph_strains_53,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Brady Genomes % Abundance, n = 53",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/Strains_wGTDB_Abund_53.png",
         width = 8,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



# 104
strains_sylph <- read.delim("data/sylph_profile_bradystrainsGTDB_ani95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia/QC_reads/", "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/",
                         "", GenomeID)) %>%
  mutate(Same = SampleID == GenomeID)
length(unique(strains_sylph$SampleID))

# Plot (will mostly be 1 to 1 but also other squares)
missings <- d %>%
  filter(sampleID %notin% strains_sylph$SampleID)
sylph_strains_104 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  arrange("GenomeID") %>%
  column_to_rownames(var = "GenomeID") %>%
  mutate(!!!set_names(rep(list(NA), length(missings$sampleID)), missings$sampleID)) %>%
  dplyr::select(all_of(d_ai_sort$sampleID)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(d_brady_ai_sort$sampleID), Reference, everything()) %>%
  t() %>%
  as.data.frame()
ann_cols <- data.frame(row.names = colnames(sylph_strains_104), 
                       AI = d_ai_sort$AI) %>%
  mutate(StrainFinder = ifelse(rownames(.) %in% d_brady$sampleID,
                               "Yes",
                               "No"))
ann_rows <- data.frame(row.names = rownames(sylph_strains_104)) %>%
  mutate(Genome = ifelse(rownames(.) %in% d_brady$sampleID,
                         "Strain",
                         ifelse(rownames(.) == "Reference",
                                "Reference", "Other GTDB")))
ann_colors <- list(Genome = c("Strain" = "#00BFC4",
                              "Reference" = "#F8766D",
                              "Other GTDB" = "#C77CFF"),
                   StrainFinder = c(No = "#F8766D",
                                    Yes = "#619CFF"))
pheatmap(sylph_strains_104,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Brady Genomes % Abundance, n = 104",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 5,
         na_col = "white",
         border_color = "white",
         filename = "InitialFigs/Strains_wGTDB_Abund_104.png",
         width = 8,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 8. Sylph 1081 n 53 ####
# Analyze Bradyrhizobium community composition across Australia
# Use taxonomic profile from Sylph with 1081 Bradyrhizobium genomes
# 212 Strains from StrainFinder plus 869 from full GTDB
# Note on Sylph output:
# Taxonomic_abundance: normalized taxonomic abundance as a percentage. Coverage-normalized - same as MetaPhlAn abundance
# Sequence_abundance: normalized sequence abundance as a percentage. The "percentage of reads" assigned to each genome - same as Kraken abundance



#### _Setup ####
strains_sylph_raw <- read.delim("data/sylph_profile_brady1081.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/Brady1081/",
                         "", Genome_file))
strains_sylph <- read.delim("data/sylph_profile_brady1081.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/POGENOM/Input_POGENOM/RAW_DATA/Reads/Aus/", 
                         "", Sample_file)) %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/", 
                         "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1"), sep = "_") %>%
  dplyr::select(-Junk1) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/Brady1081/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna.gz",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia/Brady1081/",
                         "", GenomeID)) %>%
  #mutate(Same = SampleID == GenomeID) %>% # Need grep not ==
  filter(SampleID %in% d_brady$sampleID)
length(unique(strains_sylph$SampleID)) # 53 samples
length(unique(strains_sylph$GenomeID)) # 118 genomes - need to get these and build tree!
brady118 <- as.data.frame(unique(strains_sylph_raw$GenomeID))
# write.table(brady118,
#             "data/brady118.txt", sep = "\t", row.names = F, col.names = F)
hist(strains_sylph$Eff_cov)
plot(strains_sylph$Sequence_abundance, strains_sylph$Taxonomic_abundance)

# Get info on these genomes (the 26 GTDB genomes)
gtdb_26_names <- strains_sylph %>%
  filter(grepl("GCA_|Reference", Genome_file)) %>%
  filter(!duplicated(Genome_file)) %>%
  mutate(GenomeID = gsub("Reference", "GCA_016616885.1", GenomeID))
gtdb_26 <- bacGT %>%
  subset(grepl(paste(gtdb_26_names$GenomeID, collapse="|"), ncbi_genbank_assembly_accession))
table(gtdb_26$ncbi_isolation_source) # 12 soil, 3 nodule, 11 unknown
table(gtdb_26$ncbi_genome_category)
gtdb_26$ncbi_taxonomy
gtdb_26$gtdb_taxonomy

# Plot - use sequence abundance
sylph_strains_53 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(grepl(paste(d_brady$sampleID, collapse = "|"), GenomeID),
                         1,
                         ifelse(GenomeID == "Reference",
                                2, 3))) %>%
  separate(GenomeID, remove = F, into = c("Col1", "Col2", "Col3")) %>%
  mutate(Col1 = gsub("Reference", "aReference", Col1)) %>%
  arrange(Type, match(Col1, d_brady_ai_sort$sampleID)) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(d_brady_ai_sort$sampleID))
table(d_brady$vegetation_type)
ann_cols <- data.frame(row.names = colnames(sylph_strains_53), 
                       Aridity = d_brady_ai_sort$AI,
                       Vegetation = d_brady_ai_sort$vegetation_type)
ann_rows <- data.frame(row.names = rownames(sylph_strains_53)) %>%
  mutate(Genome = ifelse(grepl(paste(d_brady$sampleID, collapse = "|"), rownames(.)),
                         "Strain",
                         ifelse(rownames(.) == "Reference",
                                "Reference", "Other GTDB")))
ann_colors <- list(Genome = c("Strain" = "#00BFC4",
                              "Reference" = "#F8766D",
                              "Other GTDB" = "#C77CFF"),
                   Vegetation = c("Forest" = "darkgreen",
                                  "Grassland" = "gold",
                                  "Woodland" = "brown"))
pheatmap(sylph_strains_53,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Bradyrhizobium Genomes % Abundance, n = 53",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         labels_row = gsub("_", " ", rownames(sylph_strains_53)),
         filename = "InitialFigs/Brady_1081_abund.png",
         width = 8,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Import Tree
# Fasttree
brady118_tree <- read.tree("data/brady118_fasttree.tree")
plot(brady118_tree)
ggtree(brady118_tree)
brady118_tree$tip.label
brady118_tree$tip.label <- gsub(".fna", "", brady118_tree$tip.label)
brady118_tree$tip.label <- gsub("Brady_", "", brady118_tree$tip.label)
brady118_tree$tip.label <- gsub("_", " ", brady118_tree$tip.label)
brady118_tree$tip.label
tipcols <- c(rep("red", 34),
             rep("black", 15),
             rep("red", 4),
             "blue",
             rep("black", 10),
             rep("red", 54))
pdf("InitialFigs/Brady118_Fasttree_Bac120.pdf", width = 9, height = 7)
ggtree(brady118_tree, linewidth = 0.1) +
  geom_tiplab(size = 1.5, vjust = 0.5,
              color = tipcols) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 1)
dev.off()

# IQ-Tree
brady118_tree <- read.tree("data/brady118_iqtree_ml.tree")
plot(brady118_tree)
ggtree(brady118_tree)
brady118_tree$tip.label
brady118_tree$tip.label <- gsub(".fna", "", brady118_tree$tip.label)
brady118_tree$tip.label <- gsub("Brady_", "", brady118_tree$tip.label)
brady118_tree$tip.label <- gsub("_", " ", brady118_tree$tip.label)
brady118_tree$tip.label
tipcols <- c(rep("red", 44),
             "blue",
             rep("black", 10),
             rep("red", 25),
             rep("black", 15),
             rep("red", 23))
pdf("InitialFigs/Brady118_IQtree_Bac120.pdf", width = 9, height = 7)
ggtree(brady118_tree, linewidth = 0.1) +
  geom_tiplab(size = 1.5, vjust = 0.5,
              color = tipcols) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 1)
dev.off()



#### _Comp ####
# Composition and drivers, Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
sp_comp <- sylph_strains_53 %>%
  replace(is.na(.), 0)
input <- list()
input$map_loaded <- d_brady_ai_sort
rownames(input$map_loaded) <- d_brady_ai_sort$sampleID
sum(d_brady_ai_sort$sampleID != names(sp_comp))
input$data_loaded <- sp_comp
input$taxonomy_loaded <- sp_comp %>%
  mutate(taxonomy1 = "Bacteria",
         taxonomy2 = "Proteobacteria",
         taxonomy3 = "Alphaproteobacteria",
         taxonomy4 = "Hyphomicrobiales",
         taxonomy5 = "Nitrobacteraceae",
         taxonomy6 = "Bradyrhizobium",
         taxonomy7 = rownames(.)) %>%
  dplyr::select(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                taxonomy5, taxonomy6, taxonomy7)
sum(rownames(input$data_loaded) != rownames(input$taxonomy_loaded))

# Alpha
input$map_loaded$rich <- specnumber(input$data_loaded, 
                                    MARGIN = 2)
input$map_loaded$shannon <- vegan::diversity(input$data_loaded, 
                                             index = "shannon", 
                                             MARGIN = 2)
range(input$map_loaded$rich)
range(input$map_loaded$shannon)
pdf("InitialFigs/BradyComp_RichAI.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(se = F) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Dissimilarity matrices
bc <- calc_dm(input$data_loaded)
ja <- calc_dm(input$data_loaded, "jaccard")

# Import phylogenetic tree. Use Phyloseq to calculate weighted UniFrac distance
otu <- phyloseq::otu_table(input$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input$taxonomy_loaded))
map <- phyloseq::sample_data(input$map_loaded)
tree <- read_tree("data/brady118_fasttree.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree)
Wun <- distance(input.phy, 
                method = "wunifrac", 
                type = "samples")
un <- distance(input.phy, 
                method = "unifrac", 
                type = "samples")



# PCoAs (for each distance matrix)
pcoa <- cmdscale(bc, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = "pH")
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_BC.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(ja, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = "pH")
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Jac.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, -Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = -Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(Wun, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.1) %>%
  mutate(shortnames = c("Precip", "AI", "pH"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Wun.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Weighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(un, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.1) %>%
  mutate(shortnames = "pH")
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_un.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Unweighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



#### _dbRDA ####
# Full dataset, no pH
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, latitude, longitude)

mod0 <- dbrda(bc ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # AI

mod0 <- dbrda(ja ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(ja ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # Nothing

mod0 <- dbrda(Wun ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(Wun ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # precip, carbon, longitude

mod0 <- dbrda(un ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(un ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # 

# Subset dataset, pH included
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, latitude, longitude, ph)
sum(is.na(d_env))
env_nona <- d_env %>%
  drop_na()
input_sub <- filter_data(input,
                         filter_cat = "sampleID",
                         keep_vals = rownames(env_nona))
bc_sub <- calc_dm(input_sub$data_loaded)
ja_sub <- calc_dm(input_sub$data_loaded, "jaccard")
otu <- phyloseq::otu_table(input_sub$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input_sub$taxonomy_loaded))
map <- phyloseq::sample_data(input_sub$map_loaded)
tree <- read_tree("data/brady118_fasttree.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
tl <- as.data.frame(tree$tip.label) %>%
  set_names("tip.label") %>%
  filter(tip.label %notin% rownames(input_sub$data_loaded))
tree_pruned <- drop.tip(tree,
                        tip = tl$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree_pruned)
Wun_sub <- distance(input.phy, 
                    method = "wunifrac", 
                    type = "samples")
un_sub <- distance(input.phy, 
                   method = "unifrac", 
                   type = "samples")

mod0 <- dbrda(bc_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(bc_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # pH

mod0 <- dbrda(ja_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(ja_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # pH

mod0 <- dbrda(Wun_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(Wun_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # 

mod0 <- dbrda(un_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(un_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # 



#### _Varpart ####
# varpart can take 2, 3, or 4 explanatory matrices
# Partition variation into geography, climate, soil
mod <- varpart(bc, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))

mod <- varpart(ja, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))


mod <- varpart(Wun, 
                   ~ latitude + longitude, 
                   ~ AI + bio1 + bio12,
                   ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
pdf("InitialFigs/BradyComp_varpart.pdf", width = 7, height = 5)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))
dev.off()

mod <- varpart(un, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))



#### _Mantel ####
# Need to remake geog.dist and env.dist, recycle code from above with input$map_loaded
# Run Partial Mantels to control for one and test the other
# Geographic distance
dist.geog <- geosphere::distm(cbind(input$map_loaded$longitude, input$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input$map_loaded$sampleID
colnames(dist.geog) <- input$map_loaded$sampleID
dist.geog[upper.tri(dist.geog, diag = TRUE)] <- NA
hist(dist.geog)

# Environmental distance
d_brady_env <- input$map_loaded %>%
  dplyr::select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_brady_env)) {
  n_na[i] <- sum(is.na(d_brady_env[,i]))
}
n_na
d_brady_env <- d_brady_env %>%
  dplyr::select(where(~ all(!is.na(.)))) %>%
  dplyr::select(-sampleID, -latitude, -longitude)
# Conductivity, nitrate, carbon, temp, precip, AI
dist.env <- as.matrix(dist(d_brady_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- d_brady$sampleID
colnames(dist.env) <- d_brady$sampleID
dist.env[upper.tri(dist.env, diag = TRUE)] <- NA
hist(dist.env)

set.seed(100)
mantel(bc, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(bc, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(bc, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(bc, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(ja, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(ja, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(ja, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(ja, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(Wun, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(Wun, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(Wun, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(Wun, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(un, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(un, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(un, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(un, dist.geog, dist.env, permutations = 2000)

pdf("InitialFigs/BradyComp_BC_geog.pdf", width = 7, height = 5)
qplot(as.dist(dist.geog), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

qplot(as.dist(dist.geog), Wun, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Weighted UniFrac") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))



#### _GDM ####
# Again, test geog.dist and env.dist

# Remake as full, not just lower triangle
dist.env <- as.matrix(dist(d_brady_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- d_brady$sampleID
colnames(dist.env) <- d_brady$sampleID

# GDM - BC
sum(rownames(as.matrix(bc)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(bc))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(bc))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "longitude",
                           YColumn = "latitude",
                           predData = gdm.data,
                           distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE) # The algorithm was unable to fit a model to your data.

# GDM - Jaccard
sum(rownames(as.matrix(ja)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(ja))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(ja))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.jac <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "longitude",
                           YColumn = "latitude",
                           predData = gdm.data,
                           distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE) # The algorithm was unable to fit a model to your data.

# GDM - Weighted UniFrac
sum(rownames(as.matrix(Wun)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(Wun))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(Wun))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.Wun <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.Wun, geo = TRUE) # Worked!
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.Wun, geo = TRUE)

# GDM - UniFrac
sum(rownames(as.matrix(un)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(un))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(un))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.un <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                         bioFormat = 3,
                         siteColumn = "gdm.sampID",
                         XColumn = "longitude",
                         YColumn = "latitude",
                         predData = gdm.data,
                         distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.un, geo = TRUE) # The algorithm was unable to fit a model to your data.



#### 9. Sylph 993 n 104 ####
# Reran Sylph but removed Strain 4 and GTDB with > 13 missing bac120
# Run on all 104 samples not just 53
# Make tree with those detected by Sylph
# Rerun code from section 8 with the new input data, new UniFrac

#### _Setup ####
strains_sylph_raw <- read.delim("data/sylph_profile_brady993_n104.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/",
                         "", Genome_file))
strains_sylph <- read.delim("data/sylph_profile_brady993_n104.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/", 
                         "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna.gz",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/",
                         "", GenomeID))
length(unique(strains_sylph$SampleID)) # 92/104 samples
length(unique(strains_sylph$GenomeID)) # 134/993 genomes - need to get these and build tree!
sort(unique(strains_sylph$GenomeID)) # 89 strains, 
brady134 <- as.data.frame(unique(strains_sylph_raw$GenomeID))
# write.table(brady134,
#             "data/brady134.txt", sep = "\t", row.names = F, col.names = F)
hist(strains_sylph$Eff_cov)
plot(strains_sylph$Sequence_abundance, strains_sylph$Taxonomic_abundance)

# Get info on these genomes (the 45 GTDB genomes) and Brady genomes in general
gtdb_brady <- bacGT %>%
  filter(grepl("Bradyrhizobium", gtdb_taxonomy))
table(gtdb_brady$ncbi_genome_category) # 85 mags, 784 isolates
table(gtdb_brady$ncbi_isolation_source)
table(gtdb_brady$gtdb_representative) # 278 reps
hist(gtdb_brady$gc_percentage)
hist(gtdb_brady$genome_size)
gtdb_45_names <- strains_sylph %>%
  filter(grepl("GCA_|Reference", Genome_file)) %>%
  filter(!duplicated(Genome_file)) %>%
  mutate(GenomeID = gsub("Reference", "GCA_016616885.1", GenomeID))
gtdb_45 <- bacGT %>%
  subset(grepl(paste(gtdb_45_names$GenomeID, collapse="|"), ncbi_genbank_assembly_accession))
table(gtdb_45$ncbi_isolation_source) # 13 soil, 17 nodule, 1 feces, 14 unknown
table(gtdb_45$ncbi_genome_category)
gtdb_45$ncbi_taxonomy
gtdb_45$gtdb_taxonomy
bac120_sub$name3 <- gsub("Reference", "GCA_016616885.1", bac120_sub$name2)
tax <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% bac120_sub$name3) %>%
  dplyr::select(ncbi_genbank_assembly_accession, gtdb_taxonomy) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species")) %>%
  dplyr::select(ncbi_genbank_assembly_accession, Species)

# Plot - use sequence abundance
d_sylph <- d %>%
  filter(sampleID %in% strains_sylph$SampleID)
d_sylph_ai_sort <- d_sylph %>%
  arrange(AI)
sylph_strains_104 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(grepl(paste(d_sylph$sampleID, collapse = "|"), GenomeID),
                       1,
                       ifelse(GenomeID == "Reference",
                              2, 3))) %>%
  separate(GenomeID, remove = F, into = c("Col1", "Col2", "Col3")) %>%
  mutate(Col1 = gsub("Reference", "aReference", Col1)) %>%
  arrange(Type, match(Col1, d_sylph_ai_sort$sampleID)) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(d_sylph_ai_sort$sampleID))
table(d_sylph_ai_sort$vegetation_type)
ann_cols <- data.frame(row.names = colnames(sylph_strains_104), 
                       Aridity = d_sylph_ai_sort$AI,
                       Vegetation = d_sylph_ai_sort$vegetation_type)
ann_rows <- data.frame(row.names = rownames(sylph_strains_104)) %>%
  mutate(Genome = ifelse(grepl(paste(d_sylph$sampleID, collapse = "|"), rownames(.)),
                         "Strain",
                         ifelse(grepl("Reference", rownames(.)),
                                "Reference", "Other GTDB")))
table(ann_rows$Genome) # 92 Strain, 89 GTDB
ann_colors <- list(Genome = c("Strain" = "#00BFC4",
                              "Reference" = "#F8766D",
                              "Other GTDB" = "#C77CFF"),
                   Vegetation = c("Forest" = "darkgreen",
                                  "Shrubland" = "chartreuse3",
                                  "Woodland" = "brown",
                                  "Grassland" = "gold",
                                  "Savannah" = "darkgoldenrod",
                                  "Heathland" = "purple"))
pheatmap(sylph_strains_104,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Bradyrhizobium Genomes % Abundance, n = 92/104",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         labels_row = gsub("_", " ", rownames(sylph_strains_104)),
         filename = "InitialFigs/Brady_993_abund104.png",
         width = 8,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Big Tree
sylph_detected <- brady134 %>%
  set_names("GenomeID") %>%
  mutate(GenomeID = gsub(".fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = gsub("Brady_", "", GenomeID))
#brady993_tree <- read.tree("data/brady993_fasttree.tree")
brady993_tree <- read.tree("data/brady993_fasttree_nogap.tree")
ggtree(brady993_tree)
brady993_tree$tip.label
brady993_tree$tip.label <- gsub(".fna", "", brady993_tree$tip.label)
brady993_tree$tip.label <- gsub("Brady_", "", brady993_tree$tip.label)
sum(sylph_detected$GenomeID %in% brady993_tree$tip.label) # 113
tipcols <- ifelse(brady993_tree$tip.label %in% sylph_detected$GenomeID,
                  "Detected",
                  "Not detected")
brady993_tree$tip.cols <- tipcols
tc <- data.frame(GenomeID = brady993_tree$tip.label,
                 Sylph = brady993_tree$tip.cols) %>%
  left_join(., tax, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  mutate(Species = coalesce(Species, GenomeID)) %>%
  mutate(Species = gsub("Reference", "Bradyrhizobium diazoefficiens_F (Reference)", Species))
tcd <- tc %>%
  filter(Sylph == "Detected")
pdf("InitialFigs/Brady993_Fasttree_Bac120_2409AA.pdf", width = 8.5, height = 17)
ggtree(brady993_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 0.5, vjust = 0.5, aes(color = Sylph, label = Species)) +
  #geom_tippoint(aes(shape = Sylph, color = Sylph)) +
  scale_color_manual(values = c("red", "grey50")) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()

# Small Tree
brady134_tree <- read.tree("data/brady134_fasttree_nogap.tree")
ggtree(brady134_tree)
brady134_tree$tip.label
brady134_tree$tip.label <- gsub(".fna", "", brady134_tree$tip.label)
brady134_tree$tip.label <- gsub("Brady_", "", brady134_tree$tip.label)
sum(sylph_detected$GenomeID %in% brady134_tree$tip.label) # 134
tc <- data.frame(GenomeID = brady134_tree$tip.label) %>%
  left_join(., tax, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  mutate(Species = coalesce(Species, GenomeID)) %>%
  mutate(Species = gsub("Reference", "Bradyrhizobium diazoefficiens_F (Reference)", Species))
pdf("InitialFigs/Brady134_Fasttree_Bac120_2971AA.pdf", width = 8.5, height = 17)
ggtree(brady134_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 1, vjust = 0.5, aes(label = Species)) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()



#### _Comp ####
# Composition and drivers, Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
sp_comp <- sylph_strains_104 %>%
  replace(is.na(.), 0)
input <- list()
input$map_loaded <- d_sylph_ai_sort
rownames(input$map_loaded) <- d_sylph_ai_sort$sampleID
sum(d_sylph_ai_sort$sampleID != names(sp_comp))
input$data_loaded <- sp_comp
input$taxonomy_loaded <- sp_comp %>%
  mutate(taxonomy1 = "Bacteria",
         taxonomy2 = "Proteobacteria",
         taxonomy3 = "Alphaproteobacteria",
         taxonomy4 = "Hyphomicrobiales",
         taxonomy5 = "Nitrobacteraceae",
         taxonomy6 = "Bradyrhizobium",
         taxonomy7 = rownames(.)) %>%
  dplyr::select(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                taxonomy5, taxonomy6, taxonomy7)
sum(rownames(input$data_loaded) != rownames(input$taxonomy_loaded))

# Alpha
input$map_loaded$rich <- specnumber(input$data_loaded, 
                                    MARGIN = 2)
input$map_loaded$shannon <- vegan::diversity(input$data_loaded, 
                                             index = "shannon", 
                                             MARGIN = 2)
range(input$map_loaded$rich)
range(input$map_loaded$shannon)
pdf("InitialFigs/BradyComp_RichAI.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(se = F) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Dissimilarity matrices
bc <- calc_dm(input$data_loaded)
ja <- calc_dm(input$data_loaded, "jaccard")

# Import phylogenetic tree. Use Phyloseq to calculate weighted UniFrac distance
otu <- phyloseq::otu_table(input$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input$taxonomy_loaded))
map <- phyloseq::sample_data(input$map_loaded)
tree <- read_tree("data/brady134_fasttree_nogap.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree)
Wun <- distance(input.phy, 
                method = "wunifrac", 
                type = "samples")
un <- distance(input.phy, 
               method = "unifrac", 
               type = "samples")



# PCoAs (for each distance matrix)
pcoa <- cmdscale(bc, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("NO3", "C", "Temp", "Precip", "AI", "pH", "Long."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_BC.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(ja, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("NO3", "C", "Temp", "Precip", "AI", "pH", "Long."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Jac.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, -Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = -Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(Wun, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.1) %>%
  mutate(shortnames = c("C", "Temp", "Precip", "AI", "pH", "Long."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Wun.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Weighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(un, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.06) %>%
  mutate(shortnames = c("C", "Temp", "Precip", "AI", "pH", "Long."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_un.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Unweighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



#### _dbRDA ####
# Full dataset, no pH
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, latitude, longitude)

mod0 <- dbrda(bc ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, C, NO3

mod0 <- dbrda(ja ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(ja ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, C, NO3

mod0 <- dbrda(Wun ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(Wun ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, C

mod0 <- dbrda(un ~ 1, d_env)  # Model with intercept only
mod1 <- dbrda(un ~ ., d_env)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, C

# Subset dataset, pH included
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, latitude, longitude, ph)
sum(is.na(d_env))
env_nona <- d_env %>%
  drop_na()
input_sub <- filter_data(input,
                         filter_cat = "sampleID",
                         keep_vals = rownames(env_nona))
bc_sub <- calc_dm(input_sub$data_loaded)
ja_sub <- calc_dm(input_sub$data_loaded, "jaccard")
otu <- phyloseq::otu_table(input_sub$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input_sub$taxonomy_loaded))
map <- phyloseq::sample_data(input_sub$map_loaded)
tree <- read_tree("data/brady134_fasttree_nogap.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
tl <- as.data.frame(tree$tip.label) %>%
  set_names("tip.label") %>%
  filter(tip.label %notin% rownames(input_sub$data_loaded))
tree_pruned <- drop.tip(tree,
                        tip = tl$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree_pruned)
Wun_sub <- distance(input.phy, 
                    method = "wunifrac", 
                    type = "samples")
un_sub <- distance(input.phy, 
                   method = "unifrac", 
                   type = "samples")

mod0 <- dbrda(bc_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(bc_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, NO3, pH, cond.

mod0 <- dbrda(ja_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(ja_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, precip, AI, pH, NO3

mod0 <- dbrda(Wun_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(Wun_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # pH, temp, AI

mod0 <- dbrda(un_sub ~ 1, env_nona)  # Model with intercept only
mod1 <- dbrda(un_sub ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # pH, temp, AI



#### _Varpart ####
# varpart can take 2, 3, or 4 explanatory matrices
# Partition variation into geography, climate, soil
mod <- varpart(bc, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))

mod <- varpart(ja, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))


mod <- varpart(Wun, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
pdf("InitialFigs/BradyComp_varpart.pdf", width = 7, height = 5)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))
dev.off()

mod <- varpart(un, 
               ~ latitude + longitude, 
               ~ AI + bio1 + bio12,
               ~ conductivity + nitrate_nitrogen + organic_carbon,
               data = d_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil'))



#### _Mantel ####
# Need to remake geog.dist and env.dist, recycle code from above with input$map_loaded
# Run Partial Mantels to control for one and test the other
# Geographic distance
dist.geog <- geosphere::distm(cbind(input$map_loaded$longitude, input$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input$map_loaded$sampleID
colnames(dist.geog) <- input$map_loaded$sampleID
dist.geog[upper.tri(dist.geog, diag = TRUE)] <- NA
hist(dist.geog)

# Environmental distance
d_sylph_env <- input_sub$map_loaded %>%
  dplyr::select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_sylph_env)) {
  n_na[i] <- sum(is.na(d_sylph_env[,i]))
}
n_na
d_sylph_env <- d_sylph_env %>%
  dplyr::select(where(~ all(!is.na(.)))) %>%
  dplyr::select(-sampleID, -latitude, -longitude)
# Conductivity, nitrate, carbon, temp, precip, AI
dist.env <- as.matrix(dist(d_sylph_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- d_sylph$sampleID
colnames(dist.env) <- d_sylph$sampleID
dist.env[upper.tri(dist.env, diag = TRUE)] <- NA
hist(dist.env)

set.seed(100)
mantel(bc, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(bc, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(bc, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(bc, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(ja, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(ja, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(ja, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(ja, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(Wun, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(Wun, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(Wun, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(Wun, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(un, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(un, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(un, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(un, dist.geog, dist.env, permutations = 2000)

pdf("InitialFigs/BradyComp_BC_geog.pdf", width = 7, height = 5)
qplot(as.dist(dist.geog), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

qplot(as.dist(dist.geog), Wun, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Weighted UniFrac") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))



#### _GDM ####
# Again, test geog.dist and env.dist

# Remake as full, not just lower triangle
dist.env <- as.matrix(dist(d_sylph_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- d_sylph$sampleID
colnames(dist.env) <- d_sylph$sampleID

# GDM - BC
sum(rownames(as.matrix(bc)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(bc))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(bc))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "longitude",
                           YColumn = "latitude",
                           predData = gdm.data,
                           distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE) # Worked
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.bray, geo = TRUE)
# Geographic         22.126
# matrix_1           25.634

# GDM - Jaccard
sum(rownames(as.matrix(ja)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(ja))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(ja))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.jac <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE) # Worked
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.jac, geo = TRUE)
# Geographic         18.690
# matrix_1           32.055

# GDM - Weighted UniFrac
sum(rownames(as.matrix(Wun)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(Wun))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(Wun))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.Wun <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.Wun, geo = TRUE) # Worked!
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.Wun, geo = TRUE)
# Geographic         19.825
# matrix_1           21.734

# GDM - UniFrac
sum(rownames(as.matrix(un)) != input$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(un))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(un))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input$map_loaded$gdm.sampID <- input$map_loaded$sampleID
gdm.data <- dplyr::select(input$map_loaded, gdm.sampID, latitude, longitude)
gdm.un <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                         bioFormat = 3,
                         siteColumn = "gdm.sampID",
                         XColumn = "longitude",
                         YColumn = "latitude",
                         predData = gdm.data,
                         distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.un, geo = TRUE) # Worked!
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.un, geo = TRUE)
# Geographic         15.846
# matrix_1           38.547



#### 10. Sylph 993 n 331 ####
# Reran Sylph but removed Strain 4 and GTDB with > 13 missing bac120
# Run on all 331 BASE metagenomes that met criteria!
# Make tree with those detected by Sylph
# Rerun code from section 8/9 with the new input data, new UniFrac
d_331 <- read.csv("data/metadata_331.csv") %>%
  mutate("Climate Class" = ifelse(AI < 0.5, "Arid to semi-arid",
                                  ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                         "Humid")))
table(d_331$`Climate Class`)
coords_trans <- st_as_sf(d_331, 
                         coords = c('longitude', 'latitude'), 
                         crs=4326)
sf_oz <- ozmap("states")
map <- ggplot(data = sf_oz) + 
  geom_sf(fill = "grey90", color = "white") +
  geom_sf(data = coords_trans,
          aes(fill = AI, shape = `Climate Class`), 
          size = 3, alpha = 1, color = "black", stroke = 0.3) +
  scale_shape_manual(values = c(21, 24, 22)) +
  annotation_scale() +
  annotation_north_arrow(pad_x = unit(1, "cm"), pad_y = unit(1, "cm"),
                         height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  guides(shape = guide_legend(order = 2),
         fill = guide_colorbar(order = 1)) +
  xlim(111, 155) +
  ylim(43, 12) +
  labs(fill = "Aridity\nindex") +
  theme_minimal()
pdf("FinalFigs/Figure1.pdf", width = 7, height = 5)
map
dev.off()
png("FinalFigs/Figure1.png", width = 7, height = 5, units = "in", res = 300)
map
dev.off()

bacGT <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") # takes a while to load
bac120 <- read.delim("data/brady_1081_markers_summary.tsv") %>%
  mutate(name2 = gsub(".fna", "", name)) %>%
  mutate(name2 = gsub("Brady_", "", name2)) %>%
  mutate(name2 = gsub("Strain_", "Strain", name2)) %>%
  separate(name2, remove = F, into = c("sampleID", "strainID")) %>%
  mutate(strainID = ifelse(grepl("Strain", strainID) == TRUE,
                           strainID, "GTDB")) %>%
  mutate(strainID = factor(strainID,
                           levels = c("Strain1", "Strain2", "Strain3", "Strain4",
                                      "GTDB")))
bac120_sub <- bac120 %>%
  filter(number_missing_genes <= 13) %>%
  mutate(name3 = gsub("Reference", "GCA_016616885.1", name2))

ref_ani <- read.delim("data/Brady_fastani.txt",
                      header = F) %>%
  set_names("Reference", "Query", "ANI", "Aligned1", "Aligned2") %>%
  mutate(Reference = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                          "", Reference)) %>%
  mutate(Query = gsub("/scratch/alpine/clbd1748/Australia/BradyStrainsRef/",
                      "", Query)) %>%
  filter(grepl("Reference", Reference)) %>%
  filter(Query != "BradyReference.fna") %>%
  separate(Query, into = c("Brady", "sampleID"), sep = "_", remove = FALSE) %>%
  mutate(sampleID = gsub(".fna", "", sampleID)) %>%
  mutate(sampleID = as.integer(sampleID)) %>%
  dplyr::select(sampleID, ANI)
range(ref_ani$ANI)

d_brady <- d %>%
  filter(sampleID %in% ref_ani$sampleID) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID)

# Check aridity
# <0.03 Hyper Arid 
# 0.03–0.2 Arid 
# 0.2–0.5 Semi-Arid 
# 0.5–0.65 Dry sub-humid 
# >0.65 Humid
sum(d_331$AI < 0.03) # 0
sum(d_331$AI >= 0.03 & d_331$AI < 0.2) # 48
sum(d_331$AI >= 0.2 & d_331$AI < 0.5) # 55
sum(d_331$AI >= 0.5 & d_331$AI < 0.65) # 110
sum(d_331$AI > 0.65) # 118

# So total:
# 103 arid to semi-arid
# 110 dry sub-humid
# 118 humid



#### _Setup ####
strains_sylph_raw <- read.delim("data/sylph_profile_brady993_n331.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/",
                         "", Genome_file))
strains_sylph <- read.delim("data/sylph_profile_brady993_n331.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/", 
                         "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna.gz",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady993/",
                         "", GenomeID))
length(unique(strains_sylph$SampleID)) # 271/331 samples
length(unique(strains_sylph$GenomeID)) # 182/993 genomes - need to get these and build tree!
sort(unique(strains_sylph$GenomeID)) # 89 strains
brady182 <- as.data.frame(unique(strains_sylph_raw$GenomeID))
# write.table(brady182,
#             "data/brady182.txt", sep = "\t", row.names = F, col.names = F)
hist(strains_sylph$Eff_cov)
plot(strains_sylph$Sequence_abundance, strains_sylph$Taxonomic_abundance)

# Get info on these genomes (the 93 GTDB genomes) and Brady genomes in general
gtdb_brady <- bacGT %>%
  filter(grepl("Bradyrhizobium", gtdb_taxonomy))
table(gtdb_brady$ncbi_genome_category) # 85 mags, 784 isolates
table(gtdb_brady$ncbi_isolation_source)
table(gtdb_brady$gtdb_representative) # 278 reps
hist(gtdb_brady$gc_percentage)
hist(gtdb_brady$genome_size)
gtdb_93_names <- strains_sylph %>%
  filter(grepl("GCA_|Reference", Genome_file)) %>%
  filter(!duplicated(Genome_file)) %>%
  mutate(GenomeID = gsub("Reference", "GCA_016616885.1", GenomeID))
gtdb_93 <- bacGT %>%
  subset(grepl(paste(gtdb_93_names$GenomeID, collapse="|"), ncbi_genbank_assembly_accession))
table(gtdb_93$ncbi_isolation_source) # 20 soil, 43 nodule, 30 other
table(gtdb_93$ncbi_genome_category) # 9 MAG, 84 not
gtdb_93$ncbi_taxonomy
gtdb_93$gtdb_taxonomy
tax <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% bac120_sub$name3) %>%
  dplyr::select(ncbi_genbank_assembly_accession, gtdb_taxonomy) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species")) %>%
  dplyr::select(ncbi_genbank_assembly_accession, Species)

# CheckM for all 993 input. Already checked strains. Check GTDB
bacGT_993 <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% strains_sylph$GenomeID)
min(bacGT_993$checkm2_completeness) # 80.94
max(bacGT_993$checkm2_contamination) # 5.15

# Plot - use sequence abundance, filter d to samples with >= 1 genome detected
d_sylph <- d_331 %>%
  filter(sampleID %in% strains_sylph$SampleID)

# Check numbers
sum(d_sylph$AI < 0.03) # 0
sum(d_sylph$AI >= 0.03 & d_sylph$AI < 0.2) # 40
sum(d_sylph$AI >= 0.2 & d_sylph$AI < 0.5) # 45
sum(d_sylph$AI >= 0.5 & d_sylph$AI < 0.65) # 95
sum(d_sylph$AI > 0.65) # 91
# For d_sylph, 85, 95, 91.

d_sylph_ai_sort <- d_sylph %>%
  arrange(AI)
sylph_strains_331 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(grepl(paste(d_sylph$sampleID, collapse = "|"), GenomeID),
                       1,
                       ifelse(GenomeID == "Reference",
                              2, 3))) %>%
  separate(GenomeID, remove = F, into = c("Col1", "Col2", "Col3")) %>%
  mutate(Col1 = gsub("Reference", "aReference", Col1)) %>%
  arrange(Type, match(Col1, d_sylph_ai_sort$sampleID)) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(as.character(d_sylph_ai_sort$sampleID)))
table(d_sylph_ai_sort$vegetation_type)
ann_cols <- data.frame(row.names = colnames(sylph_strains_331), 
                       Aridity = d_sylph_ai_sort$AI,
                       Vegetation = d_sylph_ai_sort$vegetation_type) %>%
  mutate(StrainFinder = ifelse(rownames(.) %in% d_brady$sampleID,
                               "Yes",
                               "No"))
ann_rows <- data.frame(row.names = rownames(sylph_strains_331)) %>%
  mutate(Genome = ifelse(grepl(paste(paste0(d_sylph$sampleID, "_"), collapse = "|"), 
                               rownames(.)),
                         "Strain",
                         ifelse(grepl("Reference", rownames(.)),
                                "Reference", "Other GTDB")))
table(ann_rows$Genome) # 89 Strain, 93 GTDB
ann_colors <- list(Genome = c("Strain" = "#00BFC4",
                              "Reference" = "#F8766D",
                              "Other GTDB" = "#C77CFF"),
                   Vegetation = c("Forest" = "darkgreen",
                                  "Shrubland" = "chartreuse3",
                                  "Woodland" = "brown",
                                  "Grassland" = "gold",
                                  "Savannah" = "darkgoldenrod",
                                  "Heathland" = "purple",
                                  "Dune" = "antiquewhite"),
                   StrainFinder = c(No = "#F8766D",
                                    Yes = "#619CFF"))
pheatmap(sylph_strains_331,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "Bradyrhizobium % Abundance (n genomes=182/993, n samples=271/331)",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_row = 4,
         fontsize_col = 4,
         na_col = "white",
         border_color = "white",
         labels_row = gsub("_", " ", rownames(sylph_strains_331)),
         filename = "InitialFigs/Brady_993_abund331.png",
         width = 12,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

sylph_strains_331_nona <- sylph_strains_331 %>%
  replace(is.na(.), 0)
pheatmap(sylph_strains_331_nona,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         main = "            Bradyrhizobium Genomes % Abundance, n = 271/331",
         cluster_rows = T,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_number = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "white",
         border_color = "white",
         labels_row = gsub("_", " ", rownames(sylph_strains_331)),
         filename = "InitialFigs/Brady_993_abund331_clust.png",
         width = 8,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Tree (for big context tree, see section 8)
sylph_detected <- brady182 %>%
  set_names("GenomeID") %>%
  mutate(GenomeID = gsub(".fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = gsub("Brady_", "", GenomeID))
brady182_tree <- read.tree("data/brady182_fasttree_nogap.tree")
ggtree(brady182_tree)
brady182_tree$tip.label
brady182_tree$tip.label <- gsub(".fna", "", brady182_tree$tip.label)
brady182_tree$tip.label <- gsub("Brady_", "", brady182_tree$tip.label)
sum(sylph_detected$GenomeID %in% brady182_tree$tip.label) # 182
tc <- data.frame(GenomeID = brady182_tree$tip.label) %>%
  left_join(., tax, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  mutate(Species = coalesce(Species, GenomeID)) %>%
  mutate(Species = gsub("Reference", "Bradyrhizobium diazoefficiens_F (Reference)", Species))
pdf("InitialFigs/Brady182_Fasttree_Bac120_2215AA.pdf", width = 8.5, height = 11)
ggtree(brady182_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 1, vjust = 0.5, aes(label = Species)) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.75) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()

# Tree with Aridity info
# Need genomeID rownames, then column for min AI and max AI
d_sylph_ai <- d_sylph %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(SampleID, AI)

tree_data <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID) %>%
  left_join(., d_sylph_ai, by = "SampleID") %>%
  group_by(GenomeID) %>%
  summarise(minAI = min(AI),
            maxAI = max(AI)) %>%
  as.data.frame() %>%
  left_join(., tc, by = "GenomeID")
rownames(tree_data) <- tree_data$GenomeID
sum(tree_data$GenomeID %in% brady182_tree$tip.label)
p <- ggtree(brady182_tree)
get_taxa_name(p)
tree_data <- tree_data[order(match(tree_data$GenomeID, rev(get_taxa_name(p)))), ]
tree_data$GenomeID <- factor(tree_data$GenomeID,
                             levels = rev(get_taxa_name(p)))

p <- ggtree(brady182_tree, linewidth = 0.2) +
  theme(plot.margin = margin(0,-20,0,-15))
p

s <- ggplot(tree_data) +
  geom_point(aes(x = minAI, y = GenomeID), size = 1, colour = "red") +
  geom_point(aes(x = maxAI, y = GenomeID), size = 1, colour = "blue") +
  geom_segment(aes(x = minAI, xend = maxAI,
                   y = GenomeID, yend = GenomeID)) +
  labs(x = "Aridity index",
       y = NULL) +
  #scale_y_discrete(expand = c(0.01, 0.01)) +
  scale_y_discrete(labels = tree_data$Species) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 3),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
s

pdf("InitialFigs/BradyTree_Aridity.pdf", width = 8, height = 8)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()



#### _Comp ####
# Composition and drivers, Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
sp_comp <- sylph_strains_331 %>%
  replace(is.na(.), 0)
input <- list()
input$map_loaded <- d_sylph_ai_sort
rownames(input$map_loaded) <- d_sylph_ai_sort$sampleID
sum(d_sylph_ai_sort$sampleID != names(sp_comp))
input$data_loaded <- sp_comp
input$taxonomy_loaded <- sp_comp %>%
  mutate(taxonomy1 = "Bacteria",
         taxonomy2 = "Proteobacteria",
         taxonomy3 = "Alphaproteobacteria",
         taxonomy4 = "Hyphomicrobiales",
         taxonomy5 = "Nitrobacteraceae",
         taxonomy6 = "Bradyrhizobium",
         taxonomy7 = rownames(.)) %>%
  dplyr::select(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                taxonomy5, taxonomy6, taxonomy7)
sum(rownames(input$data_loaded) != rownames(input$taxonomy_loaded))

# Alpha
input$map_loaded$rich <- specnumber(input$data_loaded, 
                                    MARGIN = 2)
input$map_loaded$shannon <- vegan::diversity(input$data_loaded, 
                                             index = "shannon", 
                                             MARGIN = 2)
range(input$map_loaded$rich)
range(input$map_loaded$shannon)
pdf("InitialFigs/BradyComp_RichAI.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(se = F) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Dissimilarity matrices
bc <- calc_dm(input$data_loaded)
ja <- calc_dm(input$data_loaded, "jaccard")

# Import phylogenetic tree. Use Phyloseq to calculate weighted UniFrac distance
otu <- phyloseq::otu_table(input$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input$taxonomy_loaded))
map <- phyloseq::sample_data(input$map_loaded)
tree <- read_tree("data/brady182_fasttree_nogap.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree)
Wun <- distance(input.phy, 
                method = "wunifrac", 
                type = "samples")
un <- distance(input.phy, 
               method = "unifrac", 
               type = "samples")



# PCoAs (for each distance matrix)
pcoa <- cmdscale(bc, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("NO3", "C", "Temp", "Precip", "AI", "pH", "Lat", "Long"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_BC.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.4) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(ja, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("NO3", "C", "Temp", "Precip", "AI", "pH", "Lat", "Long"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Jac.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(Wun, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("C", "Temp", "AI", "pH", "Long."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Wun.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Weighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Need to remove 4 outliers
out <- input$map_loaded %>%
  filter(Axis01 > 0.5)
input_filt <- filter_data(input,
                          filter_cat = "sampleID",
                          filter_vals = out$sampleID)


pcoa <- cmdscale(un, k = nrow(input$map_loaded) - 1, eig = T)
d_env <- input$map_loaded %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, 
                bio1, bio12, AI, ph, latitude, longitude)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Cond", "NO3", "C", "Temp", "Precip", "AI", "pH", "Lat", "Long"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_un.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "blue") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Unweighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()


#### _n219 ####
# Subset to 219 samples with no NA in 16 soil variables
d_sylph_env <- input$map_loaded %>%
  filter(is.na(ph) == FALSE) %>%
  dplyr::select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_sylph_env)) {
  n_na[i] <- sum(is.na(d_sylph_env[,i]))
}
n_na

# Remove variables with correlations > 0.7
d_sylph_env <- d_sylph_env %>%
  column_to_rownames(var = "sampleID") %>%
  dplyr::select(where(~ all(!is.na(.)))) %>%
  dplyr::select(-latitude, -longitude) %>% # Can add later
  dplyr::select(-boron_hot_cacl2, -exc_sodium, -sulphur) %>% # Correlated w cond.
  dplyr::select(-bio12) # Correlated with AI
names(d_sylph_env)

#plot(d_sylph_env$conductivity, d_sylph_env$boron_hot_cacl2)
#plot(d_sylph_env$conductivity, d_sylph_env$exc_sodium)
#plot(d_sylph_env$conductivity, d_sylph_env$sulphur)
# plot(d_sylph_env$bio12, d_sylph_env$AI)
# plot(d_sylph_env$bio1, d_sylph_env$AI)
m <- cor(d_sylph_env)
pdf("InitialFigs/Env_Corrplot219.pdf", width = 8, height = 6)
corrplot(m, 
         method = "number",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5,
         number.cex = 0.5)
dev.off()

# Remake inputs and dissimilarity matrices
input_sub <- filter_data(input,
                         filter_cat = "sampleID",
                         keep_vals = rownames(d_sylph_env))
bc_sub <- calc_dm(input_sub$data_loaded)
ja_sub <- calc_dm(input_sub$data_loaded, "jaccard")
otu <- phyloseq::otu_table(input_sub$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input_sub$taxonomy_loaded))
map <- phyloseq::sample_data(input_sub$map_loaded)
tree <- read_tree("data/brady182_fasttree_nogap.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
tl <- as.data.frame(tree$tip.label) %>%
  set_names("tip.label") %>%
  filter(tip.label %notin% rownames(input_sub$data_loaded))
tree_pruned <- drop.tip(tree,
                        tip = tl$tip.label)
tree_pruned$tip.label
input.phy <- phyloseq::phyloseq(otu, tax, map, tree_pruned)
Wun_sub <- distance(input.phy, 
                    method = "wunifrac", 
                    type = "samples")
un_sub <- distance(input.phy, 
                   method = "unifrac", 
                   type = "samples")

# Check distribution of values
qplot(bc_sub, Wun_sub) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25))
# Wun has values above 1! Need to fix
range(Wun_sub) # 0 to 1.123804
# Diagnose
any(phy_tree(input.phy)$edge.length < 0) # FALSE, good, no negative branch lengths
# Normalize
Wun_sub <- Wun_sub / max(Wun_sub)
range(Wun_sub) # 0 to 1
qplot(bc_sub, Wun_sub)
range(un_sub) # Good
qplot(ja_sub, Wun_sub)
qplot(un_sub, Wun_sub)

# bc, ja, un all suffer from lots of 1s. 
# Use Weighted UniFrac! It was developed for this purpose!



#### __dbRDA ####
mod0 <- dbrda(bc_sub ~ 1, d_sylph_env)  # Model with intercept only
mod1 <- dbrda(bc_sub ~ ., d_sylph_env)  # Model with all explanatory variables
set.seed(100)
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, Mg, Mn, Zn, Al, pH, NO3, AI, C, K, Cond

mod0 <- dbrda(ja_sub ~ 1, d_sylph_env)  # Model with intercept only
mod1 <- dbrda(ja_sub ~ ., d_sylph_env)  # Model with all explanatory variables
set.seed(100)
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, Mg, Mn, Al, Zn, pH, NO3, AI, C, K, Cond 

mod0 <- dbrda(Wun_sub ~ 1, d_sylph_env)  # Model with intercept only
mod1 <- dbrda(Wun_sub ~ ., d_sylph_env)  # Model with all explanatory variables
set.seed(100)
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # P, C
# Wasn't working, updated vegan and did outside test, then worked.
#saveRDS(Wun_sub, "~/Desktop/dist.rds")
#saveRDS(d_sylph_env, "~/Desktop/env.rds")

mod0 <- dbrda(un_sub ~ 1, d_sylph_env)  # Model with intercept only
mod1 <- dbrda(un_sub ~ ., d_sylph_env)  # Model with all explanatory variables
set.seed(100)
mod <- ordistep(mod0, scope = formula(mod1))
mod$anova # temp, Mg, pH, Mn, Zn, Al, NO3, AI, Fe, P, C



#### __PCoA ####
# PCoAs (for each distance matrix)
pcoa <- cmdscale(bc_sub, k = nrow(input_sub$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_sylph_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Cu", "Fe", "Mn", "Zn", "Al", "Ca", "Mg", "K", "NO3",
                        "C", "pH", "P", "Temp", "AI"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_sub$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_sub$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_bc_sub.pdf", width = 7, height = 5)
ggplot(input_sub$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.4) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(ja_sub, k = nrow(input_sub$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_sylph_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Cu", "Fe", "Mn", "Zn", "Al", "Ca", "Mg", "K", "NO3",
                        "C", "pH", "P", "Temp", "AI"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_sub$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_sub$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_jac_sub.pdf", width = 7, height = 5)
ggplot(input_sub$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



pcoa <- cmdscale(Wun_sub, k = nrow(input_sub$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_sylph_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Cu", "Mn", "Zn", "Al", "Ca", "K", "NO3",
                        "C", "P", "Temp", "Aridity"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_sub$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_sub$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_Wun_sub.pdf", width = 7, height = 5)
ggplot(input_sub$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df, shortnames %notin% c("Zn", "Aridity")),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  geom_text(data = subset(vec.df, shortnames %in% c("Zn")),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", nudge_y = -0.02) +
  geom_text(data = subset(vec.df, shortnames %in% c("Aridity")),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red", nudge_y = 0.02) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

pcoa <- cmdscale(un_sub, k = nrow(input_sub$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_sylph_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("Cond", "Cu", "Fe", "Mn", "Zn", "Al", "Ca", "Mg", "K", "NO3",
                        "C", "pH", "P", "Temp", "AI"))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_sub$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_sub$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
pdf("InitialFigs/BradyComp_un_sub.pdf", width = 7, height = 5)
ggplot(input_sub$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  ggtitle("Unweighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()



#### __Varpart ####
# varpart can take 2, 3, or 4 explanatory matrices
# Partition variation into geography, climate, soil
var_env <- d_sylph_env %>%
  mutate(latitude = input_sub$map_loaded$latitude,
         longitude = input_sub$map_loaded$longitude)
mod <- varpart(bc_sub, 
               ~ latitude + longitude, 
               ~ bio1 + AI,
               ~ conductivity + dtpa_copper + dtpa_iron + dtpa_manganese + dtpa_zinc + exc_aluminium + exc_calcium + exc_magnesium + exc_potassium + nitrate_nitrogen + organic_carbon + ph + phosphorus_colwell,
               data = var_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil')) # 18%

mod <- varpart(ja_sub, 
               ~ latitude + longitude, 
               ~ bio1 + AI,
               ~ conductivity + dtpa_copper + dtpa_iron + dtpa_manganese + dtpa_zinc + exc_aluminium + exc_calcium + exc_magnesium + exc_potassium + nitrate_nitrogen + organic_carbon + ph + phosphorus_colwell,
               data = var_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil')) # 15%


mod <- varpart(Wun_sub, 
               ~ latitude + longitude, 
               ~ bio1 + AI,
               ~ conductivity + dtpa_copper + dtpa_iron + dtpa_manganese + dtpa_zinc + exc_aluminium + exc_calcium + exc_magnesium + exc_potassium + nitrate_nitrogen + organic_carbon + ph + phosphorus_colwell,
               data = var_env)
mod
summary(mod)
pdf("InitialFigs/BradyComp_varpart.pdf", width = 7, height = 5)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil')) # 7%
dev.off()

mod <- varpart(un_sub, 
               ~ latitude + longitude, 
               ~ bio1 + AI,
               ~ conductivity + dtpa_copper + dtpa_iron + dtpa_manganese + dtpa_zinc + exc_aluminium + exc_calcium + exc_magnesium + exc_potassium + nitrate_nitrogen + organic_carbon + ph + phosphorus_colwell,
               data = var_env)
mod
summary(mod)
plot(mod, bg = 2:4, Xnames = c('Geog.', 'Clim.', 'Soil')) # 20

# Weighted UniFrac has the least % variation explained
# This might be because the others have artifcats with lots of 1s so less variation



#### __Mantel ####
# Need to remake geog.dist and env.dist, recycle code from above with input$map_loaded
# Run Partial Mantels to control for one and test the other
# Geographic distance
dist.geog <- geosphere::distm(cbind(input_sub$map_loaded$longitude, 
                                    input_sub$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input_sub$map_loaded$sampleID
colnames(dist.geog) <- input_sub$map_loaded$sampleID
#dist.geog[upper.tri(dist.geog, diag = TRUE)] <- NA
#hist(dist.geog)

# Environmental distance
dist.env <- as.matrix(dist(d_sylph_env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID
#dist.env[upper.tri(dist.env, diag = TRUE)] <- NA
#hist(dist.env)

set.seed(100)
mantel(bc_sub, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(bc_sub, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(bc_sub, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(bc_sub, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(ja_sub, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(ja_sub, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(ja_sub, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(ja_sub, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(Wun_sub, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(Wun_sub, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(Wun_sub, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(Wun_sub, dist.geog, dist.env, permutations = 2000)

set.seed(100)
mantel(un_sub, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(un_sub, dist.env, dist.geog, permutations = 2000)
set.seed(100)
mantel(un_sub, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(un_sub, dist.geog, dist.env, permutations = 2000)

pdf("InitialFigs/BradyComp_BC_geog.pdf", width = 7, height = 5)
qplot(as.dist(dist.geog), bc_sub, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

qplot(as.dist(dist.geog), Wun_sub, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Weighted UniFrac") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

qplot(as.dist(dist.geog), un_sub, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Weighted UniFrac") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
# Says missing or non-finite value...
Wun_sub_mat <- as.matrix(Wun_sub)
un_sub_mat <- as.matrix(un_sub)
sum(is.na(un_sub_mat))
sum(is.infinite(un_sub_mat))
sum(is.finite(un_sub_mat))
sum(is.finite(Wun_sub_mat))



#### __GDM ####
# Again, test geog.dist and env.dist

# Remake as full, not just lower triangle
dist.env <- as.matrix(dist(d_sylph_env, method = "euclidean", diag = T, upper = T))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID

# GDM - BC
sum(rownames(as.matrix(bc_sub)) != input_sub$map_loaded$sampleID)
gdm.sampID <- as.numeric(rownames(as.matrix(bc_sub)))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(bc_sub))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input_sub$map_loaded$gdm.sampID <- input_sub$map_loaded$sampleID
gdm.data <- dplyr::select(input_sub$map_loaded, gdm.sampID, latitude, longitude)
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "longitude",
                           YColumn = "latitude",
                           predData = gdm.data,
                           distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE) # Worked
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.bray, geo = TRUE)
#Geographic         80.498
#matrix_1           15.055

# GDM - Jaccard
sum(rownames(as.matrix(ja_sub)) != input_sub$map_loaded$sampleID)
gdm.sampID <- as.numeric(rownames(as.matrix(ja_sub)))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(ja_sub))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input_sub$map_loaded$gdm.sampID <- input_sub$map_loaded$sampleID
gdm.data <- dplyr::select(input_sub$map_loaded, gdm.sampID, latitude, longitude)
gdm.jac <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE) # Worked
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.jac, geo = TRUE)
# Geographic         78.328
# matrix_1           17.992

# GDM - Weighted UniFrac
# Try with fewer env variables since didn't converge, or has 0 sum coeff.
d_sylph_env2 <- d_sylph_env %>%
  dplyr::select(conductivity, nitrate_nitrogen, organic_carbon, ph,
                phosphorus_colwell, bio1, AI)
dist.env2 <- as.matrix(dist(d_sylph_env2, method = "euclidean", diag = T, upper = T))
rownames(dist.env2) <- input_sub$map_loaded$sampleID
colnames(dist.env2) <- input_sub$map_loaded$sampleID

sum(rownames(as.matrix(Wun_sub)) != input_sub$map_loaded$sampleID)
gdm.sampID <- as.numeric(rownames(as.matrix(Wun_sub)))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(Wun_sub))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env2)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input_sub$map_loaded$gdm.sampID <- input_sub$map_loaded$sampleID
gdm.data <- dplyr::select(input_sub$map_loaded, gdm.sampID, latitude, longitude)
gdm.Wun <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.Wun, geo = TRUE) # Worked!
summary(gdm.1)
plot(gdm.1)
# Variable Importance
gdm.varImp(gdm.Wun, geo = TRUE)
#Geographic         89.069 (p = 0.00)
#matrix_1           12.409 (p = 0.16)

# GDM - UniFrac
sum(rownames(as.matrix(un_sub)) != input_sub$map_loaded$sampleID)
gdm.sampID <- as.numeric(rownames(as.matrix(un_sub)))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(un_sub))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env2)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input_sub$map_loaded$gdm.sampID <- input_sub$map_loaded$sampleID
gdm.data <- dplyr::select(input_sub$map_loaded, gdm.sampID, latitude, longitude)
gdm.un <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "longitude",
                          YColumn = "latitude",
                          predData = gdm.data,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.un, geo = TRUE) # Response data have values greater than 1. !? (Not true)
# Weird error, can't trouble shoot, but not needed anyway
summary(gdm.1)
plot(gdm.1)
# Variable Importance
gdm.varImp(gdm.1, geo = TRUE)
#Geographic         23.522
#matrix_1           73.438



#### 11. Sylph 915 n 331 ####
# Reran Sylph but removed Strain 3 and 4 and GTDB with > 4 missing bac120
# Run on all 331 BASE metagenomes that met criteria!
# Make tree with those detected by Sylph
# Rerun code from section 8/9 with the new input data, new UniFrac
d_331 <- read.csv("data/metadata_331.csv") %>%
  mutate("Climate Class" = ifelse(AI < 0.5, "Arid to semi-arid",
                                  ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                         "Humid")))
bacGT <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") # takes a while to load


#### Commercial
# 19 genomes from Kohlmeier et al. 2025
refID <- "GCA_016616885.1"
com <- c("GCA_021052265.1", "GCA_021052365.1", "GCA_029714225.1", "GCA_029714765.1",
         "GCA_029714345.1", "GCA_021052285.1", "GCA_025200885.1", "GCA_029714545.1",
         "GCA_025200925.1", "GCA_029761915.1", "GCA_029714325.1", "GCA_029714425.1",
         "GCA_029714305.1", "GCA_029714405.1", "GCA_029714945.1", "GCA_025200865.1",
         "GCA_021052245.1", "GCA_021052305.1", "GCA_025200905.1")
comGTDB <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% com)
com_ubiq <- ubiq %>%
  filter(GenomeID %in% comGTDB$ncbi_genbank_assembly_accession)



#### _Setup ####
contig_check <- read.delim("data/sylph_profile_brady915_n331.tsv") %>%
  group_by(Sample_file, Genome_file) %>%
  summarise(count = n()) # 1 contig per genome
strains_sylph_raw <- read.delim("data/sylph_profile_brady915_n331.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady915/",
                         "", Genome_file))
strains_sylph <- read.delim("data/sylph_profile_brady915_n331.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/", 
                         "", SampleID)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady915/Brady_",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna.gz",
                         "", GenomeID)) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/Brady915/",
                         "", GenomeID))
length(unique(strains_sylph$SampleID)) # 268/331 samples
length(unique(strains_sylph$GenomeID)) # 181/915 genomes - need to get these and build tree!
sort(unique(strains_sylph$GenomeID)) # 92 strains
brady181 <- as.data.frame(unique(strains_sylph_raw$GenomeID))
sylph_detected <- brady181 %>%
  set_names("GenomeID") %>%
  mutate(GenomeID = gsub(".fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = gsub("Brady_", "", GenomeID))
sum(comGTDB$ncbi_genbank_assembly_accession %in% sylph_detected$GenomeID) # 6 commercial
# write.table(brady181,
#             "data/brady181.txt", sep = "\t", row.names = F, col.names = F)
#sum(brady181$`unique(strains_sylph_raw$GenomeID)` %notin% brady182$V1) # 8
#sum(brady182$V1 %notin% brady181$`unique(strains_sylph_raw$GenomeID)`) # 9
hist(strains_sylph$Eff_cov)
plot(strains_sylph$Sequence_abundance, strains_sylph$Taxonomic_abundance)

# Get info on these genomes (the 89 GTDB genomes)
gtdb_89_names <- strains_sylph %>%
  filter(grepl("GCA_|Reference", Genome_file)) %>%
  filter(!duplicated(Genome_file)) %>%
  mutate(GenomeID = gsub("Reference", "GCA_016616885.1", GenomeID))
gtdb_89 <- bacGT %>%
  subset(grepl(paste(gtdb_89_names$GenomeID, collapse="|"), ncbi_genbank_assembly_accession))
table(gtdb_89$ncbi_isolation_source) # 18 soil, 41 nodule, 30 other
table(gtdb_89$ncbi_genome_category) # 5 MAG, 84 not
par(mar = c(3,3,3,3))
hist(gtdb_89$genome_size)
hist(gtdb_89$gc_percentage)
range(gtdb_89$gc_percentage)
gtdb_89$ncbi_taxonomy
gtdb_89$gtdb_taxonomy
tax <- bacGT %>%
  subset(grepl(paste(gtdb_89_names$GenomeID, collapse="|"), ncbi_genbank_assembly_accession)) %>%
  dplyr::select(ncbi_genbank_assembly_accession, gtdb_taxonomy) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species")) %>%
  dplyr::select(ncbi_genbank_assembly_accession, Species)

# CheckM for all 915 input. Already checked strains. Check GTDB
# Also how many Brady?
brady <- bacGT %>%
  filter(grepl("g__Bradyrhizobium", gtdb_taxonomy))
table(brady$ncbi_isolation_source)
table(brady$ncbi_genome_category) # 85 MAG, 784 none
bacGT_915 <- bacGT %>%
  filter(ncbi_genbank_assembly_accession %in% strains_sylph$GenomeID)
range(bacGT_915$genome_size)
range(bacGT_915$checkm2_completeness) # 99.53
range(bacGT_915$checkm2_contamination) # 4.78
range(bacGT_915$checkm_completeness) # 96.15
range(bacGT_915$checkm_contamination) # 4.77

# Plot - use sequence abundance, filter d to samples with >= 1 genome detected
d_sylph <- d_331 %>%
  filter(sampleID %in% strains_sylph$SampleID)
#write.csv(d_sylph$sampleID, "data/sampleID_n268.csv")

# Check numbers
sum(d_sylph$AI < 0.03) # 0
sum(d_sylph$AI >= 0.03 & d_sylph$AI < 0.2) # 40
sum(d_sylph$AI >= 0.2 & d_sylph$AI < 0.5) # 45
sum(d_sylph$AI >= 0.5 & d_sylph$AI < 0.65) # 94
sum(d_sylph$AI > 0.65) # 89
# For d_sylph, 85, 94, 89.

d_sylph_ai_sort <- d_sylph %>%
  arrange(AI)
sylph_strains_331 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  column_to_rownames(var = "GenomeID") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(rev(tree_data_maxmin$GenomeID))) %>% # Run tree section first
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(as.character(d_sylph_ai_sort$sampleID)))
gtdb89_use <- gtdb_89 %>%
  dplyr::select(ncbi_genbank_assembly_accession, genome_size, ncbi_genome_category, 
                ncbi_isolate, ncbi_isolation_source, ncbi_assembly_level,
                ncbi_refseq_category, ncbi_genome_representation,
                gtdb_genome_representative) %>%
  mutate(Type = ifelse(ncbi_assembly_level %in% c("Chromosome", "Complete Genome"),
                       "GTDB Isolate", "GTDB MAG")) %>%
  mutate(Source = ifelse(grepl("root|nodule|plant", ncbi_isolation_source),
                         "Plant", "Other/NA")) %>%
  mutate(Source = ifelse(grepl("soil|rhizo", ncbi_isolation_source),
                         "Soil", Source)) %>%
  dplyr::select(ncbi_genbank_assembly_accession, genome_size, Type, Source)
gtdb_checkM <- gtdb_89 %>%
  dplyr::select(ncbi_genbank_assembly_accession, checkm_completeness) %>%
  rename(GenomeID = ncbi_genbank_assembly_accession)
checkm_181 <- checkM_all %>%
  filter(name == "Completeness") %>%
  filter(Strain %in% c("Strain_1", "Strain_2")) %>%
  group_by(Bin.Name) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(GenomeID = paste(sampleID, Strain, sep = "_")) %>%
  rename(checkm_completeness = value) %>%
  dplyr::select(GenomeID, checkm_completeness) %>%
  filter(GenomeID %in% sylph_detected$GenomeID) %>%
  rbind(., gtdb_checkM) %>%
  mutate(GenomeID = gsub("GCA_016616885.1", "Reference", GenomeID))
row_dat <- data.frame(GenomeID = rownames(sylph_strains_331)) %>%
  left_join(., gtdb89_use, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  replace_na(list(genome_size = 8085095, Type = "StrainFinder", Source = "Soil")) %>%
  mutate(Type = ifelse(GenomeID == "Reference", "StrainFinder Ref", Type)) %>%
  mutate(Type = ifelse(GenomeID %in% com_ubiq$GenomeID, "Commercial", Type)) %>% # Add later
  left_join(., checkm_181, by = "GenomeID") %>%
  mutate(GenomeSize = 100 * genome_size / checkm_completeness) %>%
  left_join(., brady_func, by = "GenomeID") %>%
  mutate(Function = ifelse(Nfix >= 5 & Nod >= 4, "N fix. Sym.",
                           ifelse(Nfix >= 5 & Nod < 4, "N fix. Free",
                                  ifelse(Photo > 0, "Photosyn.", " "))))
#saveRDS(row_dat, "data/row_data.rds")
row_dat <- readRDS("data/row_data.rds")
table(row_dat$Function)
ann_cols <- data.frame(row.names = colnames(sylph_strains_331), 
                       Aridity = d_sylph_ai_sort$AI) %>%
  mutate(Vegetation = d_sylph_ai_sort$vegetation_type) %>%
  #mutate(`Climate Class` = d_sylph_ai_sort$`Climate Class`) %>%
  mutate(StrainFinder = ifelse(rownames(.) %in% d_brady$sampleID,
                               "Yes",
                               "No"))
ann_rows <- data.frame(row.names = rownames(sylph_strains_331)) %>%
  mutate(`Genome Type` = row_dat$Type,
         `Genome Source` = row_dat$Source,
         #`Genome Size` = row_dat$GenomeSize
         `Function` = row_dat$Function)
table(ann_rows$Function)
ann_colors <- list(`Genome Type` = c("Commercial" = "#EE3377",
                                     "GTDB Isolate" = "#66CCEE",
                                     "GTDB MAG" = "#332288",
                                     "StrainFinder" = "#EE7733",
                                     "StrainFinder Ref" = "yellow"),
                   # Genome = c("Strain" = "#00BFC4",
                   #            "Reference" = "red",
                   #            "Other GTDB" = "#C77CFF"),
                   `Genome Source` = c("Plant" = "#44AA99",
                                       "Soil" = "#DDCC77",
                                       "Other/NA" = "#DDDDDD"),
                   `Function` = c("N fix. Sym." = "#CCEEFF",
                                  "N fix. Free" = "red",
                                  "Photosyn." = "#225522",
                                  " " = "white"),
                   # `Climate Class` = c("Arid to semi-arid" = "tan",
                   #           "Dry sub-humid" = "yellowgreen",
                   #           "Humid" = "deepskyblue"),
                   Vegetation = c("Forest" = "darkgreen",
                                  "Shrubland" = "chartreuse3",
                                  "Woodland" = "brown",
                                  "Grassland" = "gold",
                                  "Savannah" = "darkgoldenrod",
                                  "Heathland" = "purple",
                                  "Dune" = "antiquewhite"),
                   StrainFinder = c(No = "#F8766D",
                                    Yes = "#619CFF"),
                   Aridity = colorRampPalette(brewer.pal(n = 3, name = "RdYlBu"))(100))
tree_data_maxmin <- tree_data_maxmin %>%
  mutate(SpeciesShort = gsub("Bradyrhizobium ", "", Species))
pheatmap(sylph_strains_331,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 1,
         cellheight = 3,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = F,
         fontsize = 6,
         fontsize_row = 3,
         na_col = "white",
         border_color = "white",
         labels_row = rev(tree_data_maxmin$SpeciesShort),
         show_colnames = F,
         filename = "FinalFigs/Figure3.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



### _Tree ####
# Tree (for big context tree, see section 9)
brady181_tree <- read.tree("data/brady181_fasttree_nogap.tree")
brady181_tree$tip.label
brady181_tree$tip.label <- gsub(".fna", "", brady181_tree$tip.label)
brady181_tree$tip.label <- gsub("Brady_", "", brady181_tree$tip.label)
sum(sylph_detected$GenomeID %in% brady181_tree$tip.label) # 181
is.rooted(brady181_tree)
brady181_tree <- ape::root(brady181_tree, outgroup = "Nitrobacter", resolve.root = TRUE)
is.rooted(brady181_tree)
ggtree(brady181_tree)
tc <- data.frame(GenomeID = brady181_tree$tip.label) %>%
  left_join(., tax, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  mutate(Species = coalesce(Species, GenomeID)) %>%
  mutate(Species = gsub("Reference", "Bradyrhizobium diazoefficiens_F (Reference)", Species))
pdf("InitialFigs/brady181_Fasttree_Bac120_2026AA.pdf", width = 8.5, height = 11)
ggtree(brady181_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 1, vjust = 0.5, aes(label = Species)) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.75) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()

# Tree with Aridity info
# Need genomeID rownames, then column for min AI and max AI
d_sylph_ai <- d_sylph %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(SampleID, AI)
tree_data_maxmin <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID) %>%
  left_join(., d_sylph_ai, by = "SampleID") %>%
  group_by(GenomeID) %>%
  summarise(minAI = min(AI),
            maxAI = max(AI)) %>%
  as.data.frame() %>%
  left_join(., tc, by = "GenomeID") %>%
  mutate(SpeciesShort = gsub("Bradyrhizobium ", "", Species)) %>%
  #left_join(., row_dat, by = "GenomeID") %>%
  mutate(SpeciesShort = gsub(" \\(Reference\\)", "", SpeciesShort))
# Or all AI data
tree_data <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID) %>%
  left_join(., d_sylph_ai, by = "SampleID") %>%
  as.data.frame() %>%
  left_join(., tc, by = "GenomeID") %>%
  left_join(., row_dat, by = "GenomeID")
sum(tree_data$GenomeID %in% brady181_tree$tip.label)
brady181_tree <- ape::drop.tip(brady181_tree, "Nitrobacter")
p <- ggtree(brady181_tree)
get_taxa_name(p)
tree_data_maxmin <- tree_data_maxmin[order(match(tree_data_maxmin$GenomeID, rev(get_taxa_name(p)))), ]
tree_data$GenomeID <- factor(tree_data$GenomeID,
                             levels = rev(get_taxa_name(p)))

p <- ggtree(brady181_tree, linewidth = 0.2) +
  theme(plot.margin = margin(0,-20,0,-15))
p

# s <- ggplot(tree_data) +
#   geom_vline(xintercept = 0.5, linetype = "dotted") +
#   geom_vline(xintercept = 0.65, linetype = "dotted") +
#   geom_point(aes(x = AI, y = GenomeID), size = 0.5, colour = "black") +
#   geom_point(data = tree_data_maxmin, 
#              aes(x = minAI, y = GenomeID), size = 1, colour = "red") +
#   geom_point(data = tree_data_maxmin, 
#              aes(x = maxAI, y = GenomeID), size = 1, colour = "blue") +
#   geom_segment(data = tree_data_maxmin, 
#                aes(x = minAI, xend = maxAI, y = GenomeID, yend = GenomeID)) +
#   labs(x = "Aridity index",
#        y = NULL) +
#   #scale_y_discrete(expand = c(0.01, 0.01)) +
#   scale_y_discrete(labels = tree_data_maxmin$Species) +
#   theme_bw() +
#   theme(axis.text.y = element_text(size = 3),
#         axis.ticks.y = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank())
s <- ggplot(tree_data, aes(x = AI, y = GenomeID, colour = Type)) +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  geom_vline(xintercept = 0.65, linetype = "dotted") +
  geom_segment(data = tree_data_maxmin, 
               aes(x = minAI, xend = maxAI, y = GenomeID, yend = GenomeID,
                   colour = Type)) +
  geom_point(size = 0.5) +
  labs(x = "Aridity index",
       y = NULL) +
  scale_y_discrete(labels = tree_data_maxmin$SpeciesShort) +
  scale_colour_manual(values = c("#EE3377", "#66CCEE", "#332288", "#EE7733", "#FDE725FF")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 3),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
s

pdf("FinalFigs/Figure6.pdf", width = 8, height = 8)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()
png("FinalFigs/Figure6.png", width = 8, height = 8, units = "in", res = 300)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()

# Try to make 3 groupings
# Genomes detected in at least 2 samples
# All samples only in a specific climate class
one <- ubiq %>%
  filter(Ubiquity < 2) # 109 genomes
two <- ubiq %>%
  filter(Ubiquity >= 2) # 72 genomes
arid <- tree_data_maxmin %>%
  mutate(Group = "Arid to semi-arid") %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(minAI < 0.5) %>%
  filter(maxAI < 0.5) # 12 genomes
dry <- tree_data_maxmin %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(minAI >= 0.5) %>%
  filter(maxAI < 0.65) # 10 genomes
humid <- tree_data_maxmin %>%
  mutate(Group = "Humid") %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(minAI >= 0.65) %>%
  filter(maxAI >= 0.65) # 11 genomes
compare_ai <- rbind(arid, humid)
hist(compare_ai$GenomeSize)
car::leveneTest(GenomeSize ~ Group, data = compare_ai)
t.test(GenomeSize ~ Group, data = compare_ai)

general <- tree_data_maxmin %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(GenomeID %notin% arid$GenomeID) %>%
  filter(GenomeID %notin% dry$GenomeID) %>%
  filter(GenomeID %notin% humid$GenomeID) # 39
tree_data_maxmin$range <- tree_data_maxmin$maxAI - tree_data_maxmin$minAI
hist(tree_data_maxmin$range)
# Add to ubiq and plot versus size
ubiq <- ubiq %>%
  dplyr::select(GenomeID, Ubiquity) %>%
  left_join(., tree_data_maxmin, by = "GenomeID")
ggplot(ubiq, aes(GenomeSize, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  labs(x = "Estimated genome size (bp)",
       y = "Aridity range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ubiq_two <- ubiq %>%
  filter(Ubiquity >= 2)
ggplot(ubiq_two, aes(GenomeSize, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Estimated genome size (bp)",
       y = "Aridity range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggplot(ubiq_two, aes(Ubiquity, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Ubiquity",
       y = "Aridity range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
summary(lm(range ~ Ubiquity, data = ubiq_two)) # R2 = 0.34, p < 0.001
summary(lm(range ~ GenomeSize, data = ubiq_two)) # R2 = 0.01, p = 0.52



# Tree with Temperature info (instead of aridity)
# First check temps across the climate classes
ggplot(d_331, aes(`Climate Class`, bio1)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 2, pch = 16, alpha = 0.75, width = 0.2) +
  theme_bw()
hist(d_331$bio1)
# Need genomeID rownames, then column for min temp and max temp
d_sylph_temp <- d_sylph %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(SampleID, bio1)
tree_data_maxmin <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID) %>%
  left_join(., d_sylph_temp, by = "SampleID") %>%
  group_by(GenomeID) %>%
  summarise(minAI = min(bio1),
            maxAI = max(bio1)) %>%
  as.data.frame() %>%
  left_join(., tc, by = "GenomeID") %>%
  mutate(SpeciesShort = gsub("Bradyrhizobium ", "", Species)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(SpeciesShort = gsub(" \\(Reference\\)", "", SpeciesShort))
# Or all temp data
tree_data <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID) %>%
  left_join(., d_sylph_temp, by = "SampleID") %>%
  as.data.frame() %>%
  left_join(., tc, by = "GenomeID") %>%
  left_join(., row_dat, by = "GenomeID")
sum(tree_data$GenomeID %in% brady181_tree$tip.label)
brady181_tree <- ape::drop.tip(brady181_tree, "Nitrobacter")
p <- ggtree(brady181_tree)
get_taxa_name(p)
tree_data_maxmin <- tree_data_maxmin[order(match(tree_data_maxmin$GenomeID, rev(get_taxa_name(p)))), ]
tree_data$GenomeID <- factor(tree_data$GenomeID,
                             levels = rev(get_taxa_name(p)))

p <- ggtree(brady181_tree, linewidth = 0.2) +
  theme(plot.margin = margin(0,-20,0,-15))
p
# s <- ggplot(tree_data) +
#   geom_point(aes(x = bio1, y = GenomeID), size = 0.5, colour = "black") +
#   geom_point(data = tree_data_maxmin, 
#              aes(x = minT, y = GenomeID), size = 1, colour = "red") +
#   geom_point(data = tree_data_maxmin, 
#              aes(x = maxT, y = GenomeID), size = 1, colour = "blue") +
#   geom_segment(data = tree_data_maxmin, 
#                aes(x = minT, xend = maxT, y = GenomeID, yend = GenomeID)) +
#   labs(x = "MAT (\u00b0C)",
#        y = NULL) +
#   scale_y_discrete(labels = tree_data_maxmin$Species) +
#   theme_bw() +
#   theme(axis.text.y = element_text(size = 3),
#         axis.ticks.y = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank())
# s
s <- ggplot(tree_data, aes(x = bio1, y = GenomeID, colour = Type)) +
  geom_segment(data = tree_data_maxmin, 
               aes(x = minAI, xend = maxAI, y = GenomeID, yend = GenomeID,
                   colour = Type)) +
  geom_point(size = 0.5) +
  labs(x = "MAT (\u00b0C)",
       y = NULL) +
  scale_y_discrete(labels = tree_data_maxmin$SpeciesShort) +
  scale_colour_manual(values = c("#EE3377", "#66CCEE", "#332288", "#EE7733", "#FDE725FF")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0.65),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 3),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
s

pdf("FinalFigs/FigureS6.pdf", width = 8, height = 8)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()
png("FinalFigs/FigureS6.png", width = 8, height = 8, units = "in", res = 300)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()

# Try to make 2 groupings
# Genomes detected in at least 2 samples
# Narrow range vs. wide range
one <- ubiq %>%
  filter(Ubiquity < 2) # 109 genomes
two <- ubiq %>%
  filter(Ubiquity >= 2) # 72 genomes
narrow <- tree_data_maxmin %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(maxT - minT < 5) # 40 genomes
wide <- tree_data_maxmin %>%
  filter(GenomeID %in% two$GenomeID) %>%
  filter(maxT - minT > 10) # 19 genomes
tree_data_maxmin$range <- tree_data_maxmin$maxT - tree_data_maxmin$minT
hist(tree_data_maxmin$range)
# Add to ubiq and plot versus size
ubiq <- ubiq %>%
  left_join(., tree_data_maxmin, by = "GenomeID") %>%
  mutate(Group = ifelse(range < 5, "Narrow",
                        ifelse(range >=5 & range <= 10, "Intermediate",
                               "Wide")))
ggplot(ubiq, aes(GenomeSize, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  labs(x = "Estimated genome size (bp)",
       y = "Temperature range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ubiq_two <- ubiq %>%
  filter(Ubiquity >= 2)
ggplot(ubiq_two, aes(GenomeSize, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Estimated genome size (bp)",
       y = "Temperature range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggplot(ubiq_two, aes(Ubiquity, range)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Ubiquity",
       y = "Temperature range") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
summary(lm(range ~ Ubiquity, data = ubiq_two)) # R2 = 0.30, p < 0.001
summary(lm(range ~ GenomeSize, data = ubiq_two)) # R2 = 0.04, p = 0.08

# Do a third time but for pH? (range 3.1 to 8.6, 3.2 to 8.6 for Sylph)



# Big tree (subset the previous 993 tree down to the 915 in this rerun)
brady915_2 <- brady915 %>%
  mutate(GenomeID = gsub(".gz", "", Header))
toRemove <- as.data.frame(brady993_tree$tip.label) %>%
  set_names("GenomeID") %>%
  filter(GenomeID %notin% brady915_2$GenomeID)
brady993_tree <- read.tree("data/brady993_fasttree_nogap.tree")
length(brady993_tree$tip.label) # 993
brady915_tree <- ape::drop.tip(brady993_tree, toRemove$GenomeID)
length(brady915_tree$tip.label) # 915, good
brady915_tree$tip.label
brady915_tree$tip.label <- gsub(".fna", "", brady915_tree$tip.label)
brady915_tree$tip.label <- gsub("Brady_", "", brady915_tree$tip.label)
sum(sylph_detected$GenomeID %in% brady915_tree$tip.label) # 181
tipcols <- ifelse(brady915_tree$tip.label %in% sylph_detected$GenomeID,
                  "Detected",
                  "Not detected")
brady915_tree$tip.cols <- tipcols
tc <- data.frame(GenomeID = brady915_tree$tip.label,
                 Sylph = brady915_tree$tip.cols) %>%
  left_join(., tax, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
  mutate(Species = coalesce(Species, GenomeID)) %>%
  mutate(Species = gsub("Reference", "Bradyrhizobium diazoefficiens_F (Reference)", Species))
tcd <- tc %>%
  filter(Sylph == "Detected")
pdf("FinalFigs/FigureS3.pdf", width = 8.5, height = 17)
ggtree(brady915_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 0.5, vjust = 0.5, aes(color = Sylph, label = Species)) +
  #geom_tippoint(aes(shape = Sylph, color = Sylph)) +
  scale_color_manual(values = c("red", "grey50")) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()
png("FinalFigs/FigureS3.png", width = 8.5, height = 17, units = "in", res = 300)
ggtree(brady915_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 0.5, vjust = 0.5, aes(color = Sylph, label = Species)) +
  #geom_tippoint(aes(shape = Sylph, color = Sylph)) +
  scale_color_manual(values = c("red", "grey50")) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()



#### _Comp ####
# Composition and drivers, Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
sp_comp <- sylph_strains_331 %>%
  replace(is.na(.), 0)
input <- list()
input$map_loaded <- d_sylph_ai_sort
rownames(input$map_loaded) <- d_sylph_ai_sort$sampleID
sum(d_sylph_ai_sort$sampleID != names(sp_comp))
input$data_loaded <- sp_comp
input$taxonomy_loaded <- sp_comp %>%
  mutate(taxonomy1 = "Bacteria",
         taxonomy2 = "Proteobacteria",
         taxonomy3 = "Alphaproteobacteria",
         taxonomy4 = "Hyphomicrobiales",
         taxonomy5 = "Nitrobacteraceae",
         taxonomy6 = "Bradyrhizobium",
         taxonomy7 = rownames(.)) %>%
  dplyr::select(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
                taxonomy5, taxonomy6, taxonomy7)
sum(rownames(input$data_loaded) != rownames(input$taxonomy_loaded))
saveRDS(input, "data/input_sylph.rds")
input <- readRDS("data/input_sylph.rds")

# Prevalence
prev <- rowSums(input$data_loaded > 0)
range(prev)
sort(prev)

# Alpha
input$map_loaded$rich <- specnumber(input$data_loaded, 
                                    MARGIN = 2)
input$map_loaded$shannon <- vegan::diversity(input$data_loaded, 
                                             index = "shannon", 
                                             MARGIN = 2)
range(input$map_loaded$rich)
range(input$map_loaded$shannon)
pdf("InitialFigs/BradyComp_RichAI.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(se = F) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

# Check nitrate vs rich
summary(lm(rich ~ nitrate_nitrogen, data = input$map_loaded)) # p = 0.27
ggplot(input$map_loaded, aes(nitrate_nitrogen, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Nitrate",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))

# Dissimilarity matrices
bc <- calc_dm(input$data_loaded)
ja <- calc_dm(input$data_loaded, "jaccard")

# Import phylogenetic tree. Use Phyloseq to calculate weighted UniFrac distance
otu <- phyloseq::otu_table(input$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input$taxonomy_loaded))
map <- phyloseq::sample_data(input$map_loaded)
tree <- read_tree("data/brady181_fasttree_nogap.tree")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- root(tree, outgroup = "Nitrobacter.fna", resolve.root = TRUE)
is.rooted(tree) # TRUE. Tree is now rooted at Nitrobacter.
tree <- ape::drop.tip(tree, "Nitrobacter.fna")
is.rooted(tree) # Tree is still rooted even after dropping that tip.
tree$tip.label
tree$tip.label <- gsub(".fna", "", tree$tip.label)
tree$tip.label <- gsub("Brady_", "", tree$tip.label)
input.phy <- phyloseq::phyloseq(otu, tax, map, tree)
Wun <- distance(input.phy, 
                method = "wunifrac", 
                type = "samples")
un <- distance(input.phy, 
               method = "unifrac", 
               type = "samples")
set.seed(100)
mantel(bc, Wun, permutations = 2000)

set.seed(100)
mantel(Wun, un, permutations = 2000) # r = 0.6, p = 0.0005

# Check
hist(bc)
hist(ja)
hist(Wun)
hist(un)
range(Wun)

# Check environmental variables
d_env <- input$map_loaded %>%
  dplyr::select(sampleID, all_of(env_vars))
n_na <- c()
for (i in 1:ncol(d_env)) {
  n_na[i] <- sum(is.na(d_env[,i]))
}
n_na # Most NA is 67, which means 201 have, which is pretty good.
# envfit can handle NA anyway
# remove boron_hot_cacl2, exc_sodium, sulphur, bio12 because too correlated
# Also sand and iron
d_env <- input$map_loaded %>%
  dplyr::select(clay, conductivity, dtpa_copper, dtpa_manganese, 
                dtpa_zinc, exc_aluminium, exc_calcium, exc_magnesium, exc_potassium,
                latitude, longitude, nitrate_nitrogen, organic_carbon, ph, 
                phosphorus_colwell, silt, water_content, bio1, AI)
d_env2 <- d_env %>%
  dplyr::select(-latitude, -longitude)
names(d_env2)
names(d_env2) <- c("Clay", "Conductivity", "Cu", "Mn", "Zn", "Al", "Ca", "Mg", "K",
                   "NO3", "C", "pH", "P", "Silt", "H2O", "Temp.", "Aridity")
#plot(d_env$conductivity, d_env$boron_hot_cacl2)
#plot(d_env$conductivity, d_env$exc_sodium)
#plot(d_env$conductivity, d_env$sulphur)
#plot(d_env$bio12, d_env$AI)
#plot(d_env$bio1, d_env$AI)
m <- cor(d_env2, use = "pairwise.complete.obs")
pdf("FinalFigs/FigureS1.pdf", width = 8, height = 6)
corrplot(m, 
         method = "number",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5,
         number.cex = 0.5)
dev.off()
png("FinalFigs/FigureS1.png", width = 8, height = 6, units = "in", res = 300)
corrplot(m, 
         method = "number",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5,
         number.cex = 0.5)
dev.off()



# Ubiquity, genome size etc.
ubiq <- data.frame(GenomeID = rownames(input$data_loaded),
                   Ubiquity = rowSums(input$data_loaded > 0)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(Perc = Ubiquity/331*100) %>%
  mutate(Nfix = ifelse(Function == "N fix. Free" | Function == "N fix. Sym.",
                       "Yes", "No"))
table(ubiq$Function)
table(ubiq$Nfix) # 145 to 36
t.test(GenomeSize ~ Nfix, data = ubiq) # p = 0.52
ggplot(ubiq, aes(Nfix, GenomeSize)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()

t.test(Ubiquity ~ Nfix, data = ubiq) # p = 0.01
ggplot(ubiq, aes(Nfix, Ubiquity)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()
# Careful - this could be driven by the Strains! Rerun below

summary(lm(GenomeSize ~ Ubiquity, data = ubiq))
cor.test(ubiq$GenomeSize, ubiq$Ubiquity, method = "pearson")
figS4 <- ggplot(ubiq, aes(Ubiquity, GenomeSize)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  labs(x = "Prevalence (n samples detected)",
       y = "Estimated genome size (bp)") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
figS4
pdf("FinalFigs/FigureS4.pdf", width = 8, height = 6)
figS4
dev.off()
png("FinalFigs/FigureS4.png", width = 8, height = 6, units = "in", res = 300)
figS4
dev.off()

ubiq_gtdb <- ubiq %>%
  filter(Type %in% c("GTDB Isolate", "GTDB MAG", "StrainFinder Ref"))

t.test(GenomeSize ~ Nfix, data = ubiq_gtdb) # p = 0.07
ggplot(ubiq_gtdb, aes(Nfix, GenomeSize)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()

t.test(Ubiquity ~ Nfix, data = ubiq_gtdb) # p = 0.08
ggplot(ubiq_gtdb, aes(Nfix, Ubiquity)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()

summary(lm(genome_size ~ Ubiquity, data = ubiq_gtdb))
cor.test(ubiq_gtdb$genome_size, ubiq_gtdb$Ubiquity, method = "pearson")
ggplot(ubiq_gtdb, aes(Ubiquity, genome_size)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Prevalence (n samples detected)",
       y = "Genome size (bp)") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

ubiq_sort <- ubiq %>%
  arrange(desc(Ubiquity), desc(Type)) %>%
  mutate(rank = row_number())
pdf("InitialFigs/UbiquityRankAbund.pdf", width = 8, height = 6)
ggplot(ubiq_sort, aes(rank, Ubiquity)) +
  geom_point(aes(fill = Type), pch = 21, size = 3) +
  labs(x = "Rank",
       y = "Prevalence (n samples detected)") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()



#### __PCoA ####
# PCoAs (for each distance matrix)
pcoa <- cmdscale(bc, k = nrow(input$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)

pcoa <- cmdscale(ja, k = nrow(input$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)

pcoa <- cmdscale(Wun, k = nrow(input$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.001) %>%
  mutate(shortnames = c("Mn", "Zn", "Al", "Ca", "K", "Long.", "C", "P", "Silt",
                        "H2O", "Temp."))
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
input$map_loaded$ClimateClass <- input$map_loaded$`Climate Class`
micro.hulls <- ddply(input$map_loaded, "ClimateClass", find_hull)
fig4 <- ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  # geom_polygon(data = micro.hulls, aes(colour = ClimateClass, fill = ClimateClass),
  #              alpha = 0.1, show.legend = F, linewidth = NA) + # All overlapping
  geom_point(size = 3, pch = 16, alpha = 0.4) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  #ggtitle("Weighted UniFrac") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
fig4

pcoa <- cmdscale(un, k = nrow(input$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)

# NMDS
set.seed(100)
nmds <- metaMDS(Wun, trymax = 200)
stressplot(nmds)
nmds$stress
set.seed(100)
ef <- envfit(nmds, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(nmds)
plot(ef, cex = 0.5, p.max = 0.001)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(NMDS1 = NMDS1 * multiplier,
         NMDS2 = NMDS2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.001) %>%
  filter(variables != "dtpa_iron") %>%
  filter(variables != "clay") %>%
  mutate(shortnames = c("Mn", "Al", "Ca", "Mg", "K", "Long.", "C", "P",
                        "H2O", "Temp."))
input$map_loaded$Axis01 <- vegan::scores(nmds)[,1]
input$map_loaded$Axis02 <- vegan::scores(nmds)[,2]
fig4 <- ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.4, aes(colour = `Climate Class`)) +
  # geom_segment(data = vec.df,
  #              aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #              arrow = arrow(length = unit(0.35, "cm")),
  #              colour = "gray", alpha = 0.6,
  #              inherit.aes = FALSE) + 
  # geom_text_repel(data = vec.df,
  #                 aes(x = NMDS1, y = NMDS2, label = shortnames),
  #                 size = 3, color = "red") +
  labs(x = "NMDS1", 
       y = "NMDS2") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
fig4
pdf("FinalFigs/Figure4.pdf", width = 7, height = 5)
fig4
dev.off()
png("FinalFigs/Figure4.png", width = 7, height = 5, units = "in", res = 300)
fig4
dev.off()

m <- adonis2(Wun ~ input$map_loaded$`Climate Class`)
m # Sig but low R2 (0.03)
m <- betadisper(Wun, input$map_loaded$`Climate Class`)
anova(m) # Sig
m <- adonis2(Wun ~ input$map_loaded$vegetation_type)
m # Sig R2 = 0.10
m <- betadisper(Wun, input$map_loaded$vegetation_type)
anova(m) # Sig

ord <- calc_ordination(Wun, 'NMDS')
mctoolsr::plot_ordination(input, ord, color_cat = "vegetation_type", hulls = T)
mctoolsr::plot_ordination(input, ord, color_cat = "Climate Class", hulls = T)



#### __dbRDA ####
# Okay, now you can't have any NA
# Have to subset samples because all soil variables have NAs
# 216 samples retained
names(d_env2)
rownames(d_env2) <- rownames(d_env)
n_na <- c()
for (i in 1:ncol(d_env2)) {
  n_na[i] <- sum(is.na(d_env2[,i]))
}
n_na
d_env3 <- d_env2 %>%
  filter(is.na(pH) == FALSE) %>%
  dplyr::select(-Clay, -Silt, -H2O)
Wun_mat <- as.matrix(Wun)
Wun_sub_mat <- Wun_mat[rownames(d_env3), rownames(d_env3)]
Wun_sub <- as.dist(Wun_sub_mat)

input_sub <- filter_data(input,
                         filter_cat = "sampleID",
                         keep_vals = rownames(d_env3)) # 216 remaining

mod0 <- dbrda(Wun_sub ~ 1, d_env3)  # Model with intercept only
mod1 <- dbrda(Wun_sub ~ ., d_env3)  # Model with all explanatory variables
set.seed(100)
dbmod <- ordistep(mod0, scope = formula(mod1))
dbmod$anova # Temp, pH, Mn, Al, NO3, Zn



#### __Varpart ####
# varpart can take 2, 3, or 4 explanatory matrices
# Partition variation into geography vs. environment
sum(rownames(d_env3) != input_sub$map_loaded$sampleID)
var_env <- d_env3 %>%
  mutate(latitude = input_sub$map_loaded$latitude,
         longitude = input_sub$map_loaded$longitude)
mod <- varpart(Wun_sub, 
               ~ latitude + longitude, 
               ~ `Temp.` + pH + Mn + Al + NO3 + Zn,
               data = var_env)
mod
summary(mod)
par(mar = c(1, 1, 1, 1))
par(mfrow = c(1,1))
pdf("FinalFigs/FigureS4.pdf", width = 7, height = 5)
plot(mod, bg = c("#F8766D", "#619CFF"), Xnames = c('Geog.', 'Env.')) # 30%
dev.off()
png("FinalFigs/FigureS4.png", width = 7, height = 5, units = "in", res = 300)
plot(mod, bg = c("#F8766D", "#619CFF"), Xnames = c('Geog.', 'Env.')) # 30%
dev.off()

fig5a <- function() {
  par(
    mar = c(1, 2.8, 1, 0.45),
    mfrow = c(1,1)
  )
  plot(mod, bg = c("#F8766D", "#619CFF"), Xnames = c('Geog.', 'Env.'))
}
ggdraw(fig5a)
pdf("FinalFigs/Figure5.pdf", width = 7, height = 5)
plot_grid(fig5a, fig5b, ncol = 1, labels = "auto", vjust = 1, hjust = -1)
dev.off()
png("FinalFigs/Figure5.png", width = 7, height = 5, units = "in", res = 300)
plot_grid(fig5a, fig5b, ncol = 1, labels = "auto", vjust = 1, hjust = -1)
dev.off()



#### __Mantel ####
# Need to make geog.dist and env.dist
# For env.dist, only use the dbRDA-selected variables
# Geography run Mantel
# Environment run Partial Mantel controlling for geography

# Geographic distance
dist.geog <- geosphere::distm(cbind(input_sub$map_loaded$longitude, 
                                    input_sub$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input_sub$map_loaded$sampleID
colnames(dist.geog) <- input_sub$map_loaded$sampleID

# Environmental distance
d_env_sig <- d_env3 %>%
  dplyr::select(`Temp.`, pH, Mn, Al, NO3, Zn)
dist.env <- as.matrix(dist(d_env_sig, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID

set.seed(100)
mantel(Wun_sub, dist.geog, permutations = 2000)
set.seed(100)
mantel.partial(Wun_sub, dist.geog, dist.env, permutations = 2000)
set.seed(100)
mantel.partial(Wun_sub, dist.env, dist.geog, permutations = 2000) # r = 0.19, p = 0.0005

qplot(as.dist(dist.geog), Wun_sub, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Geographic Distance (m)",
       y = "Weighted UniFrac Distance") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

qplot(as.dist(dist.env), Wun_sub, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Environmental Distance",
       y = "Weighted UniFrac Distance") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

# Test some other variables
# pH
dist.env <- as.matrix(dist(input_sub$map_loaded$ph, 
                           method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID
set.seed(100)
mantel.partial(Wun_sub, dist.env, dist.geog, permutations = 2000) # r = 0.09, p = 0.003

# Temperature
# Geographic distance for full dataset
dist.geog <- geosphere::distm(cbind(input$map_loaded$longitude, 
                                    input$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input$map_loaded$sampleID
colnames(dist.geog) <- input$map_loaded$sampleID

dist.env <- as.matrix(dist(input$map_loaded$bio1, 
                           method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input$map_loaded$sampleID
colnames(dist.env) <- input$map_loaded$sampleID

set.seed(100)
mantel.partial(Wun, dist.env, dist.geog, permutations = 2000) # r = 0.31, p = 0.0005
set.seed(100)
mantel.partial(bc, dist.env, dist.geog, permutations = 2000) # r = 0.05, p = 0.0005

# Aridity
dist.env <- as.matrix(dist(input$map_loaded$AI, 
                           method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input$map_loaded$sampleID
colnames(dist.env) <- input$map_loaded$sampleID

set.seed(100)
mantel.partial(Wun, dist.env, dist.geog, permutations = 2000) # r = 0.08, p = 0.009



#### __GDM ####
# Again, test geog.dist and env.dist

# Remake as full, not just lower triangle
dist.env <- as.matrix(dist(d_env_sig, method = "euclidean", diag = T, upper = T))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID

# GDM - Weighted UniFrac
sum(rownames(as.matrix(Wun_sub)) != input_sub$map_loaded$sampleID)
gdm.sampID <- as.numeric(rownames(as.matrix(Wun_sub)))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(Wun_sub))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, dist.geog)
env.distance.v.mat.gdm <- cbind(gdm.sampID, dist.env)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(env.distance.v.mat.gdm)
input_sub$map_loaded$gdm.sampID <- input_sub$map_loaded$sampleID
gdm.data <- dplyr::select(input_sub$map_loaded, longitude,latitude)
sp_coords <- SpatialPointsDataFrame(coords = gdm.data, 
                                    data = gdm.data, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
utm_zone <- CRS("+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
utm_coords <- spTransform(sp_coords, utm_zone)
utm_x <- utm_coords@coords[, 1]
utm_y <- utm_coords@coords[, 2] 
df <- data.frame(gdm.sampID = input_sub$map_loaded$gdm.sampID,
                 X = utm_x,
                 Y = utm_y)
gdm.Wun <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "X",
                          YColumn = "Y",
                          predData = df,
                          distPreds = list(env.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.Wun, geo = TRUE) # Worked!
summary(gdm.1)
plot(gdm.1)
# Variable Importance
gdm.varImp(gdm.Wun, geo = TRUE)
#Geographic         5.321 (p = 0)
#matrix_1           97.708 (p = 0)
gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
par(mar = c(4,4,4,4))
par(mfrow = c(1,2))
plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance (m)", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_1"], gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Environmental distance", ylab="Partial ecological distance")
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
gdm.1.pred <- predict(gdm.1, gdm.Wun)
head(gdm.1.pred)
par(mfrow = c(1,1))
plot(gdm.Wun$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
# Make nice ggplot
gdm_df1 <- data.frame(xdist = gdm.1.splineDat$x[,"Geographic"],
                      ydist = gdm.1.splineDat$y[,"Geographic"],
                      variable = "Geographic (m)")
gdm_df2 <- data.frame(xdist = gdm.1.splineDat$x[,"matrix_1"],
                      ydist = gdm.1.splineDat$y[,"matrix_1"],
                      variable = "Environmental")
gdm_df <- rbind(gdm_df1, gdm_df2) %>%
  mutate(variable = factor(variable, levels = c("Geographic (m)", "Environmental")))
fig5b <- ggplot(gdm_df, aes(xdist, ydist, colour = variable)) +
  geom_line(linewidth = 2, show.legend = F) +
  labs(x = "Distance",
       y = "Partial ecological distance") +
  scale_colour_manual(values = c("#F8766D", "#619CFF")) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))
fig5b



#### 12. Brady 182 ####
# Pangenomics of the 182 detected genomes
#### _KO ####
ko <- read.delim("data/KO-PRESENCE-ABSENCE.txt")
brady_mat <- ko %>%
  dplyr::select(-key) %>%
  column_to_rownames(var = "KOfam") %>%
  set_names(gsub("Brady_", "", names(.))) %>%
  set_names(ifelse(grepl("GCA", names(.)) == TRUE,
                         gsub("_1", "\\.1", names(.)),
                         names(.))) %>%
  dplyr::select(all_of(get_taxa_name(p)))
names(brady_mat)
sum(names(brady_mat) %in% brady182_tree$tip.label)

# Get KOs in all 182, in 181/182, in 1
k182 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 182) %>%
  select(-sum) # 839
k181 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 181) %>%
  select(-sum) # 302
k1 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 1) %>%
  select(-sum) # 78

brady_mat_func <- brady_mat %>%
  filter(grepl(ignore.case = F, "Nif|photosynthetic", rownames(.))) %>%
  t() %>%
  as.data.frame()
sum(brady_mat_func$`nitrogenase iron protein NifH`)
pheatmap(brady_mat_func,
         color = c("grey80", "grey30"),
         legend = T,
         legend_breaks = c(0, 1),
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 315,
         fontsize_row = 3,
         fontsize_col = 6,
         border_color = "white",
         labels_row = rev(tree_data$Species),
         filename = "InitialFigs/Brady182_func.png",
         width = 4,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

#### _CAZyme ####
# AA 12 auxiliary activities
# CB 8 carbohydrate-binding
# CE 11 carbohydrate esterase
# GH 87 glycoside hydrolase
# GT 40 glycosyltransferase
# PL 15 polysaccharide lyase

cazy <- read.delim("data/CAZy-FREQUENCY.txt")

brady_mat <- cazy %>%
  mutate(CAZyme = gsub(".hmm", "", CAZyme)) %>%
  mutate(Class = substr(CAZyme, start = 1, stop = 2)) %>%
  mutate(Class = gsub("AA", "Auxiliary activities", Class)) %>%
  mutate(Class = gsub("CB", "Carbohydrate-binding", Class)) %>%
  mutate(Class = gsub("CE", "Carbohydrate esterase", Class)) %>%
  mutate(Class = gsub("GH", "Glycoside hydrolase", Class)) %>%
  mutate(Class = gsub("GT", "Glycosyltransferase", Class)) %>%
  mutate(Class = gsub("PL", "Polysaccharide lyase", Class)) %>%
  dplyr::select(-key, -CAZyme) %>%
  group_by(Class) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "Class") %>%
  set_names(gsub("Brady_", "", names(.))) %>%
  set_names(ifelse(grepl("GCA", names(.)) == TRUE,
                   gsub("_1", "\\.1", names(.)),
                   names(.))) %>%
  dplyr::select(all_of(get_taxa_name(p)))
names(brady_mat)
sum(names(brady_mat) %in% brady182_tree$tip.label)
brady_mat_cazy <- brady_mat %>%
  t() %>%
  as.data.frame()
pheatmap(brady_mat_cazy,
         legend = T,
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 315,
         fontsize_row = 3,
         fontsize_col = 6,
         border_color = "white",
         labels_row = rev(tree_data$Species),
         filename = "InitialFigs/Brady182_cazy.png",
         width = 4,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _Gene Clusters ####
# From anvio pangenomics analysis
# Pangenome summary file (anvi-summarize)
# Note: 50945 gene clusters with 1,387,358 genes
gc_all <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary.txt")
names(gc_all)
length(unique(gc_all$gene_cluster_id))

# Check if all genes in the gc have the same function
gc_check_fun <- gc_all %>%
  group_by(gene_cluster_id, COG20_FUNCTION_ACC) %>%
  summarise(n = n())
nrow(gc_check_fun) # 55179 is more than 50945, so no.



#### __Patterns ####
# Collapse to 1 row per gene cluster, with 1 COG category (most frequent) each
# Wrangle COG Categories. Note 2392 unknown function
gc <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary.txt") %>%
  arrange(gene_cluster_id, COG20_CATEGORY) %>%
  group_by(gene_cluster_id) %>%
  add_count(COG20_CATEGORY, sort = TRUE) %>% # Count occurrences of 'COG20_CATEGORY' within each group
  slice_max(n, n = 1, with_ties = FALSE) %>% # Select the most frequent
  #slice_head(n = 1) %>% # Old version, just getting first COG
  mutate(COG20_CATEGORY = ifelse(COG20_CATEGORY == "", "Unknown", COG20_CATEGORY)) %>%
  separate(COG20_CATEGORY, remove = F, sep = "\\|",
           into = c("COG1", "COG2", "COG3", "COG4", "COG5", "COG6", "COG7", "COG8"))
dup_fun <- gc %>%
  rowwise() %>%
  transmute(all_equal = n_distinct(c_across(c("COG1", "COG2", "COG3", "COG4", 
                                              "COG5", "COG6", "COG7", "COG8")), 
                                   na.rm = TRUE) == 1) %>%
  ungroup()
gc$all_equal <- dup_fun$all_equal
gc <- gc %>%
  mutate(COG20_Uniq = ifelse(all_equal == TRUE, COG1, "Multiple"))

# Prep data and plot % of gene clusters by COG
cog_df_all <- as.data.frame(table(gc$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate(All = Freq/nrow(gc)*100) %>%
  dplyr::select(COG, All)
gc_182 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 182)
cog_df_182 <- as.data.frame(table(gc_182$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_182" = Freq/nrow(gc_182)*100) %>%
  dplyr::select(COG, Genomes_182)
gc_2181 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 181 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_2181 <- as.data.frame(table(gc_2181$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_2_to_181" = Freq/nrow(gc_2181)*100) %>%
  dplyr::select(COG, Genomes_2_to_181)
gc_1 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 1)
cog_df_1 <- as.data.frame(table(gc_1$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_1" = Freq/nrow(gc_1)*100) %>%
  dplyr::select(COG, Genomes_1)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_182, by = c("COG")) %>%
  left_join(., cog_df_2181, by = c("COG")) %>%
  left_join(., cog_df_1, by = c("COG")) %>%
  #replace(is.na(.), 0) %>%
  column_to_rownames(var = "COG")
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
pheatmap(cog_comb,
         legend = T,
         legend_breaks = c(0.0255037, 20, 40, 60, 80, 88.99676),
         legend_labels = c("0", "20", "40", "60", "80", ""),
         main = "            % Gene Clusters by COG",
         cluster_rows = F,
         cluster_cols = F,
         labels_col = c("All (n = 50945)", 
                        "In 182 genomes (n = 1505)",
                        "In 2 to 181 genomes (n = 36654)",
                        "In 1 genome (n = 12786)"),
         angle_col = 315,
         display_numbers = T,
         number_color = "black",
         fontsize_number = 10,
         border_color = "white",
         filename = "InitialFigs/Brady182_GeneClusters_COGs.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Prep data and plot homogeneity by COG
cog_df_all <- gc %>%
  group_by(COG20_Uniq) %>%
  summarise(All_FH = mean(functional_homogeneity_index),
            All_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
cog_df_182 <- gc_182 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes182_FH = mean(functional_homogeneity_index),
            Genomes182_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_91181 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 181 & num_genomes_gene_cluster_has_hits >= 91)
cog_df_91181 <- gc_91181 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes91181_FH = mean(functional_homogeneity_index),
            Genomes91181_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_290 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 90 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_290 <- gc_290 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes290_FH = mean(functional_homogeneity_index),
            Genomes290_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_182, by = c("COG")) %>%
  left_join(., cog_df_91181, by = c("COG")) %>%
  left_join(., cog_df_290, by = c("COG")) %>%
  column_to_rownames(var = "COG") %>%
  dplyr::select(All_FH, Genomes182_FH, Genomes91181_FH, Genomes290_FH,
                All_GH, Genomes182_GH, Genomes91181_GH, Genomes290_GH)
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
ann_cols <- data.frame(row.names = colnames(cog_comb),
                       "Index" = c(rep("Functional", 4),
                                   rep("Geometric", 4)))
ann_colors <- list(Index = c(Functional = "#F8766D",
                             Geometric = "#619CFF"))
pheatmap(cog_comb,
         legend = T,
         main = "               Gene Clusters Homogeneity by COG",
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = 4,
         labels_col = c("All (n = 50945)",
                        "In 182 genomes (n = 1505)",
                        "In 91 to 181 genomes (n = 4943)",
                        "In 2 to 90 genomes (n = 31711)",
                        "All (n = 50945)",
                        "In 182 genomes (n = 1505)",
                        "In 91 to 181 genomes (n = 4943)",
                        "In 2 to 90 genomes (n = 31711)"),
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         display_numbers = T,
         number_format = "%.3f",
         number_color = "black",
         fontsize_number = 6,
         border_color = "white",
         filename = "InitialFigs/Brady182_GeneClusters_Homogeneity.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 13. Brady 181 ####
# Pangenomics of the 181 detected genomes
#### _ANI ####
# Want to quickly report the range of ANI among these 181
ani_181 <- read.delim("data/ANI_ani_181.txt") %>%
  column_to_rownames(var = "layers")
range(ani_181[ani_181 < 1]) # 0.795487 to 0.999997



#### _KO ####
ko <- read.delim("data/KO-PRESENCE-ABSENCE_181.txt", quote = "") # 3567 KOs
brady_mat <- ko %>%
  dplyr::select(-key) %>%
  column_to_rownames(var = "KOfam") %>%
  set_names(gsub("Brady_", "", names(.))) %>%
  set_names(ifelse(grepl("GCA", names(.)) == TRUE,
                   gsub("_1", "\\.1", names(.)),
                   names(.))) %>%
  dplyr::select(all_of(get_taxa_name(p)))
names(brady_mat)
sum(names(brady_mat) %in% brady181_tree$tip.label)

# Get KOs in all 181, in 180/181, in 1
k181 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 181) %>%
  select(-sum) # 1060
k180 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 180) %>%
  select(-sum) # 217
k1 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 1) %>%
  select(-sum) # 269

# Searching "nodulation" in KEGG yields:
# K12546 nodO; putative nodulation protein
# K14658 nodA; nodulation protein A [EC:2.3.1.-]
# K14660 nodE; nodulation protein E [EC:2.3.1.-]
# K14661 nodF; nodulation protein F [EC:2.3.1.-]
# K27744 NSP1; protein nodulation signaling pathway 1
# K27745 NSP2; protein nodulation signaling pathway 2
# But only A and F were found in the Brady (F only once though)

# Plus these don't have nod or nodulation in the name
# K14659 nodB; chitooligosaccharide deacetylase [EC:3.5.1.-]
# K14666 nodC; N-acetylglucosaminyltransferase [EC:2.4.1.-]
# K14657 nodD; LysR family transcriptional regulator, nod-box dependent transcriptional activator
# All 3 of these present!

# Plus Tao et al. 2021 mentions nodI and nodJ. They report nodABCIJ and nifBDEHKN.
# K09695 nodI; lipooligosaccharide transport system ATP-binding protein
# K09694 nodJ; lipooligosaccharide transport system permease protein

# Add those, and this needs to be Figure S7
brady_mat_func <- brady_mat %>%
  filter(grepl(ignore.case = F, "chitooligosaccharide deacetylase|N-acetylglucosaminyltransferase|LysR family transcriptional regulator, nod-box dependent transcriptional activator|nodulation|Nif|photosynthetic", rownames(.))) %>%
  t() %>%
  as.data.frame()
brady_func <- brady_mat_func %>%
  mutate(Nod = rowSums(brady_mat_func[, c(5, 3, 6, 2, 18)])) %>%
  mutate(Nfix = rowSums(brady_mat_func[, c(1, 4, 7:14)])) %>%
  mutate(Photo = rowSums(brady_mat_func[, 15:17])) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Nod, Nfix, Photo)
sum(brady_func$Nod > 0) # 175
sum(brady_func$Nfix > 0) # 146
sum(brady_func$Photo > 0) # 5

brady_mat_func <- brady_mat_func %>%
  dplyr::select(5, 3, 6, 2, 18,
                1, 4, 7:14,
                15:17)
names(brady_mat_func)[1:5] <- c("nodA", "nodB", "nodC", "nodD", "nodF")
names(brady_mat_func)[9:10] <- c("nitrogenase Mo-Fe protein NifN",
                                 "nitrogenase Mo-cofactor synthesis protein NifE")

sum(brady_mat_func$`nitrogen fixation protein NifZ`) # 145
sum(brady_mat_func$nodA) # 143
sum(brady_mat_func$nodF) # 1

ann_cols <- data.frame(row.names = colnames(brady_mat_func),
                       Function = c(rep("Nodulation", 5),
                                    rep("N fixation", 10),
                                    rep("Photosynthesis", 3)))
ann_rows <- data.frame(row.names = rownames(sylph_strains_331)) %>%
  mutate(`Genome Type` = row_dat$Type,
         `Genome Source` = row_dat$Source)
ann_colors <- list(`Genome Type` = c("Commercial" = "#EE3377",
                                     "GTDB Isolate" = "#66CCEE",
                                     "GTDB MAG" = "#332288",
                                     "StrainFinder" = "#EE7733",
                                     "StrainFinder Ref" = "yellow"),
                   `Genome Source` = c("Plant" = "#44AA99",
                                       "Soil" = "#DDCC77",
                                       "Other/NA" = "#DDDDDD"),
                   `Function` = c("Nodulation" = "#BBCCEE",
                                  "N fixation" = "#FFCCCC",
                                  "Photosynthesis" = "#CCDDAA"))

pheatmap(brady_mat_func,
         color = c("grey80", "grey30"),
         legend = T,
         legend_breaks = c(0, 1),
         legend_labels = c("0  ", "1  "),
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = ann_rows,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         fontsize_row = 3,
         fontsize_col = 5,
         fontsize = 6,
         border_color = "white",
         #gaps_col = c(1, 9),
         gaps_col = c(5, 15),
         labels_row = rev(tree_data_maxmin$SpeciesShort),
         filename = "FinalFigs/FigureS7.png",
         #filename = "InitialFigs/Brady181_func_wNod.png",
         width = 4,
         height = 9)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Remake using just the Tao 2021 genes plus photosynthesis. nifBDEHKN
brady_mat_func2 <- brady_mat %>%
  filter(grepl(ignore.case = F, "nodulation protein A|chitooligosaccharide deacetylase|N-acetylglucosaminyltransferase|lipooligosaccharide transport system ATP-binding protein|lipooligosaccharide transport system permease protein|nitrogen fixation protein NifB|nitrogenase molybdenum-iron protein alpha chain|nitrogenase molybdenum-cofactor synthesis protein NifE|nitrogenase iron protein NifH|nitrogenase molybdenum-iron protein beta chain|nitrogenase molybdenum-iron protein NifN|photosynthe", rownames(.))) %>%
  t() %>%
  as.data.frame()
names(brady_mat_func2)
brady_func <- brady_mat_func2 %>%
  mutate(Nod = rowSums(brady_mat_func2[, c(3, 2, 4, 1, 5)])) %>%
  mutate(Nfix = rowSums(brady_mat_func2[, c(9, 7, 8, 11, 10, 6)])) %>%
  mutate(Photo = rowSums(brady_mat_func2[, 12:14])) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Nod, Nfix, Photo)
sum(brady_func$Nod >= 4) # 143
sum(brady_func$Nfix >= 5) # 145
sum(brady_func$Photo >= 2) # 5

brady_mat_func2 <- brady_mat_func2 %>%
  dplyr::select(3, 2, 4, 1, 5,
                9, 7, 8, 11, 10, 6,
                12:14)
names(brady_mat_func2)[1:5] <- c("nodA", "nodB", "nodC", "nodI", "nodJ")
names(brady_mat_func2)[6:11] <- c("nifB", "nifD", "nifE", "nifH", "nifK", "nifN")

sum(brady_mat_func2$nifH) # 143
sum(brady_mat_func2$nodA) # 143

ann_cols <- data.frame(row.names = colnames(brady_mat_func2),
                       Function = c(rep("Nodulation", 5),
                                    rep("N fixation", 6),
                                    rep("Photosynthesis", 3)))
ann_rows <- data.frame(row.names = rownames(sylph_strains_331)) %>%
  mutate(`Genome Type` = row_dat$Type,
         `Genome Source` = row_dat$Source)
ann_colors <- list(`Genome Type` = c("Commercial" = "#EE3377",
                                     "GTDB Isolate" = "#66CCEE",
                                     "GTDB MAG" = "#332288",
                                     "StrainFinder" = "#EE7733",
                                     "StrainFinder Ref" = "yellow"),
                   `Genome Source` = c("Plant" = "#44AA99",
                                       "Soil" = "#DDCC77",
                                       "Other/NA" = "#DDDDDD"),
                   `Function` = c("Nodulation" = "#BBCCEE",
                                  "N fixation" = "#FFCCCC",
                                  "Photosynthesis" = "#CCDDAA"))
pheatmap(brady_mat_func2,
         color = c("grey80", "grey30"),
         legend = T,
         legend_breaks = c(0, 1),
         legend_labels = c("0  ", "1  "),
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = ann_rows,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         fontsize_row = 3,
         fontsize_col = 5,
         fontsize = 6,
         border_color = "white",
         gaps_col = c(5, 11),
         labels_row = rev(tree_data_maxmin$SpeciesShort),
         filename = "FinalFigs/FigureS7.png",
         width = 4,
         height = 9)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



# Plot number of KOs shared by each number of genomes
int_info_all <- data.frame(Genomes = seq(1:181)) %>%
  mutate(Genomes = factor(Genomes,
                          levels = seq(1:181))) %>%
  mutate(KOs = "NA")
for (i in 1:181) {
  k <- brady_mat %>%
    mutate(sum = rowSums(.)) %>%
    filter(sum == i) %>%
    select(-sum)
  int_info_all$KOs[i] <- nrow(k)
}
int_info_all$KOs <- as.integer(int_info_all$KOs)
pdf("InitialFigs/Genes_SharedKOs_181.pdf", width = 7, height = 4)
ggplot(int_info_all, aes(Genomes, KOs)) +
  geom_bar(stat = "identity") +
  #geom_text(aes(Genomes, KOs + 20, label = KOs), size = 2) +
  labs(x = "Number of genomes (1 to 181)",
       y = "Number of shared KOs") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank())
dev.off()



#### _CAZyme ####
# AA 12 auxiliary activities
# CB 8 carbohydrate-binding
# CE 11 carbohydrate esterase
# GH 87 glycoside hydrolase
# GT 40 glycosyltransferase
# PL 15 polysaccharide lyase

cazy <- read.delim("data/CAZy-FREQUENCY_181.txt")

brady_mat <- cazy %>%
  mutate(CAZyme = gsub(".hmm", "", CAZyme)) %>%
  mutate(Class = substr(CAZyme, start = 1, stop = 2)) %>%
  mutate(Class = gsub("AA", "Auxiliary activities", Class)) %>%
  mutate(Class = gsub("CB", "Carbohydrate-binding", Class)) %>%
  mutate(Class = gsub("CE", "Carbohydrate esterase", Class)) %>%
  mutate(Class = gsub("GH", "Glycoside hydrolase", Class)) %>%
  mutate(Class = gsub("GT", "Glycosyltransferase", Class)) %>%
  mutate(Class = gsub("PL", "Polysaccharide lyase", Class)) %>%
  dplyr::select(-key, -CAZyme) %>%
  group_by(Class) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "Class") %>%
  set_names(gsub("Brady_", "", names(.))) %>%
  set_names(ifelse(grepl("GCA", names(.)) == TRUE,
                   gsub("_1", "\\.1", names(.)),
                   names(.))) %>%
  dplyr::select(all_of(get_taxa_name(p)))
names(brady_mat)
sum(names(brady_mat) %in% brady181_tree$tip.label)
brady_mat_cazy <- brady_mat %>%
  t() %>%
  as.data.frame()
pheatmap(brady_mat_cazy,
         legend = T,
         scale = "column",
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 315,
         fontsize_row = 3,
         fontsize_col = 6,
         border_color = "white",
         labels_row = rev(tree_data$Species),
         filename = "InitialFigs/Brady181_cazy.png",
         width = 4,
         height = 10)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _Gene Clusters ####
# From anvio pangenomics analysis
# Pangenome summary file (anvi-summarize)
# Note: 79261 gene clusters with 1,416,019 genes
gc_all <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary_181.txt", quote = "")
names(gc_all)
length(unique(gc_all$gene_cluster_id)) # 79261

# Check if all genes in the gc have the same function
gc_check_fun <- gc_all %>%
  group_by(gene_cluster_id, COG20_FUNCTION_ACC) %>%
  summarise(n = n())
nrow(gc_check_fun) # 83313 is more than 79261, so no.



#### __Patterns ####
# Collapse to 1 row per gene cluster, sort by GH
gc <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary_181.txt", quote = "") %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  arrange(geometric_homogeneity_index)
gc_core <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 181)


# Collapse to 1 row per gene cluster, with 1 COG category (most frequent) each
# Wrangle COG Categories. Note 2392 unknown function
gc <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary_181.txt", quote = "") %>%
  arrange(gene_cluster_id, COG20_CATEGORY) %>%
  group_by(gene_cluster_id) %>%
  add_count(COG20_CATEGORY, sort = TRUE) %>% # Count # 'COG20_CATEGORY' in each group
  slice_max(n, n = 1, with_ties = FALSE) %>% # Select the most frequent
  mutate(COG20_CATEGORY = ifelse(COG20_CATEGORY == "", "Unknown", COG20_CATEGORY)) %>%
  separate(COG20_CATEGORY, remove = F, sep = "\\|",
           into = c("COG1", "COG2", "COG3", "COG4", "COG5", "COG6", "COG7", "COG8"))
dup_fun <- gc %>%
  rowwise() %>%
  transmute(all_equal = n_distinct(c_across(c("COG1", "COG2", "COG3", "COG4", 
                                              "COG5", "COG6", "COG7", "COG8")), 
                                   na.rm = TRUE) == 1) %>%
  ungroup()
gc$all_equal <- dup_fun$all_equal
gc <- gc %>%
  mutate(COG20_Uniq = ifelse(all_equal == TRUE, COG1, "Multiple"))

# Prep data and plot % of gene clusters by COG
cog_df_all <- as.data.frame(table(gc$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate(All = Freq/nrow(gc)*100) %>%
  dplyr::select(COG, All)
gc_181 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 181)
cog_df_181 <- as.data.frame(table(gc_181$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_181" = Freq/nrow(gc_181)*100) %>%
  dplyr::select(COG, Genomes_181)
gc_2180 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 180 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_2180 <- as.data.frame(table(gc_2180$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_2_to_180" = Freq/nrow(gc_2180)*100) %>%
  dplyr::select(COG, Genomes_2_to_180)
gc_1 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 1)
cog_df_1 <- as.data.frame(table(gc_1$COG20_Uniq)) %>%
  set_names(c("COG", "Freq")) %>%
  mutate("Genomes_1" = Freq/nrow(gc_1)*100) %>%
  dplyr::select(COG, Genomes_1)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_181, by = c("COG")) %>%
  left_join(., cog_df_2180, by = c("COG")) %>%
  left_join(., cog_df_1, by = c("COG")) %>%
  #replace(is.na(.), 0) %>%
  column_to_rownames(var = "COG")
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
pheatmap(cog_comb,
         legend = T,
         legend_breaks = c(0.001261655, 15, 30, 45, 60, 72.41355),
         legend_labels = c("0", "15", "30", "45", "60", ""),
         main = "            % Gene Clusters by COG",
         cluster_rows = F,
         cluster_cols = F,
         labels_col = c("All (n = 79261)", 
                        "In 181 genomes (n = 1935)",
                        "In 2 to 180 genomes (n = 35541)",
                        "In 1 genome (n = 41785)"),
         angle_col = 315,
         display_numbers = T,
         number_color = "black",
         fontsize_number = 10,
         border_color = "white",
         filename = "InitialFigs/Brady181_GeneClusters_COGs.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Prep data and plot homogeneity by COG
cog_df_all <- gc %>%
  group_by(COG20_Uniq) %>%
  summarise(All_FH = mean(functional_homogeneity_index),
            All_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
# 1?
cog_df_1 <- gc_1 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes1_FH = mean(functional_homogeneity_index),
            Genomes1_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
range(cog_df_1$Genomes1_FH)
range(cog_df_1$Genomes1_GH)
cog_df_181 <- gc_181 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes181_FH = mean(functional_homogeneity_index),
            Genomes181_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_91180 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 180 & num_genomes_gene_cluster_has_hits >= 91)
cog_df_91180 <- gc_91180 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes91180_FH = mean(functional_homogeneity_index),
            Genomes91180_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
gc_290 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 90 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_290 <- gc_290 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes290_FH = mean(functional_homogeneity_index),
            Genomes290_GH = mean(geometric_homogeneity_index)) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq)
cog_comb <- cog_df_all %>%
  left_join(., cog_df_181, by = c("COG")) %>%
  left_join(., cog_df_91180, by = c("COG")) %>%
  left_join(., cog_df_290, by = c("COG")) %>%
  column_to_rownames(var = "COG") %>%
  dplyr::select(Genomes181_FH, Genomes91180_FH, Genomes290_FH,
                Genomes181_GH, Genomes91180_GH, Genomes290_GH) %>%
  filter(rownames(.) != "Chromatin structure and dynamics")
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
ann_cols <- data.frame(row.names = colnames(cog_comb),
                       "Index" = c(rep("Functional", 3),
                                   rep("Geometric", 3)))
ann_colors <- list(Index = c(Functional = "#F8766D",
                             Geometric = "#619CFF"))
pheatmap(cog_comb,
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = 3,
         labels_col = c("In 181 genomes (n = 1935)",
                        "In 91 to 180 genomes (n = 4836)",
                        "In 2 to 90 genomes (n = 30705)",
                        "In 181 genomes (n = 1935)",
                        "In 91 to 180 genomes (n = 4836)",
                        "In 2 to 90 genomes (n = 30705)"),
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         display_numbers = T,
         number_format = "%.3f",
         number_color = "black",
         fontsize_number = 6,
         border_color = "white",
         filename = "InitialFigs/Brady181_GeneClusters_Homogeneity.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### __Heatmap ####
# Which are most variable? P/A or geometric?
# P/A index - how many genomes?
# How many total GC per COG class?
# Make Figure 7

# First combine the Unknown and Function unknown categories
gc <- gc %>%
  mutate(COG20_Uniq = gsub("Function unknown", "Unknown", COG20_Uniq))

gc_181 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits == 181)
cog_df_181 <- gc_181 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes181_GH = mean(geometric_homogeneity_index),
            Genomes181_nGenome = mean(num_genomes_gene_cluster_has_hits),
            Genomes181_nGC = n()) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq) %>%
  mutate(Genomes181_PA = Genomes181_nGenome/181,
         Genomes181_perGC = Genomes181_nGC/nrow(gc_181))
ggplot(cog_df_181, aes(Genomes181_nGC, Genomes181_GH)) +
  geom_point(size = 3, pch = 16, alpha = 0.8) +
  labs(x = "Number of gene clusters in COG",
       y = "Mean geometric homogeneity in COG") +
  ggtitle("181 genomes (n = 1935)") +
  theme_bw()

gc_91180 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 180 & num_genomes_gene_cluster_has_hits >= 91)
cog_df_91180 <- gc_91180 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes91180_GH = mean(geometric_homogeneity_index),
            Genomes91180_nGenome = mean(num_genomes_gene_cluster_has_hits),
            Genomes91180_nGC = n()) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq) %>%
  mutate(Genomes91180_PA = Genomes91180_nGenome/181,
         Genomes91180_perGC = Genomes91180_nGC/nrow(gc_91180))
ggplot(cog_df_91180, aes(Genomes91180_nGC, Genomes91180_GH)) +
  geom_point(size = 3, pch = 16, alpha = 0.8) +
  labs(x = "Number of gene clusters in COG",
       y = "Mean geometric homogeneity in COG") +
  ggtitle("91 to 180 genomes (n = 4836)") +
  theme_bw()

gc_290 <- gc %>%
  filter(num_genomes_gene_cluster_has_hits <= 90 & num_genomes_gene_cluster_has_hits >= 2)
cog_df_290 <- gc_290 %>%
  group_by(COG20_Uniq) %>%
  summarise(Genomes290_GH = mean(geometric_homogeneity_index),
            Genomes290_nGenome = mean(num_genomes_gene_cluster_has_hits),
            Genomes290_nGC = n()) %>%
  ungroup() %>%
  rename(COG = COG20_Uniq) %>%
  mutate(Genomes290_PA = Genomes290_nGenome/181,
         Genomes290_perGC = Genomes290_nGC/nrow(gc_290))
ggplot(cog_df_290, aes(Genomes290_nGC, Genomes290_GH)) +
  geom_point(size = 3, pch = 16, alpha = 0.8) +
  labs(x = "Number of gene clusters in COG",
       y = "Mean geometric homogeneity in COG") +
  ggtitle("2 to 90 genomes (n = 30705)") +
  theme_bw()

cog_comb <- cog_df_181 %>%
  left_join(., cog_df_91180, by = c("COG")) %>%
  left_join(., cog_df_290, by = c("COG")) %>%
  column_to_rownames(var = "COG") %>%
  dplyr::select(Genomes181_GH, Genomes181_PA, Genomes181_perGC,
                Genomes91180_GH, Genomes91180_PA, Genomes91180_perGC,
                Genomes290_GH, Genomes290_PA, Genomes290_perGC)
cog_comb_num <- cog_df_181 %>%
  left_join(., cog_df_91180, by = c("COG")) %>%
  left_join(., cog_df_290, by = c("COG")) %>%
  column_to_rownames(var = "COG") %>%
  dplyr::select(Genomes181_GH, Genomes181_PA, Genomes181_nGC,
                Genomes91180_GH, Genomes91180_PA, Genomes91180_nGC,
                Genomes290_GH, Genomes290_PA, Genomes290_nGC) %>%
  rename(Genomes181_perGC = Genomes181_nGC,
         Genomes91180_perGC = Genomes91180_nGC,
         Genomes290_perGC = Genomes290_nGC) %>%
  as.matrix()
colSums(cog_comb, na.rm = T)
min(cog_comb, na.rm = T)
max(cog_comb, na.rm = T)
ann_cols <- data.frame(row.names = colnames(cog_comb),
                       "Metric" = c("Geometric", "PA", "Num",
                                    "Geometric", "PA", "Num",
                                    "Geometric", "PA", "Num"),
                       "Genomes" = c(rep("181 genomes (n = 1935)", 3),
                                     rep("91 to 180 genomes (n = 4836)", 3),
                                     rep("2 to 90 genomes (n = 30705)", 3)))
ann_colors <- list(Metric = c(Geometric = "#DDAA33",
                              PA = "#BB5566",
                              Num = "#004488"),
                   Genomes = c("181 genomes (n = 1935)" = "#440154FF",
                               "91 to 180 genomes (n = 4836)" = "#21908CFF",
                               "2 to 90 genomes (n = 30705)" = "#FDE725FF"))
pheatmap(cog_comb,
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(3, 6),
         show_colnames = F,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         display_numbers = cog_comb_num,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6,
         border_color = "white",
         filename = "FinalFigs/Figure7.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Plot number of GC shared by each number of genomes
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
gc_genomeNum <- gc %>%
  dplyr::select(gene_cluster_id, num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  group_by(num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  summarise(nGC = n())
length(unique(gc_genomeNum$num_genomes_gene_cluster_has_hits))
pdf("InitialFigs/GeneClusters_NumGenomes.pdf", width = 8, height = 6)
ggplot(gc_genomeNum, aes(num_genomes_gene_cluster_has_hits, nGC, fill = COG20_Uniq)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of genomes (1 to 181)",
       y = "Number of gene clusters",
       fill = "COG Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        #axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank())
dev.off()

# Now plot proportions
gcSum_genomeNum <- gc_genomeNum %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(sum = sum(nGC)) %>%
  ungroup()
gcProp_genomeNum <- gc_genomeNum %>%
  left_join(., gcSum_genomeNum, by = "num_genomes_gene_cluster_has_hits") %>%
  mutate(Prop = nGC/sum)
pdf("InitialFigs/GeneClustersProp_NumGenomes.pdf", width = 8, height = 6)
ggplot(gcProp_genomeNum, aes(num_genomes_gene_cluster_has_hits, Prop, fill = COG20_Uniq)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of genomes (1 to 181)",
       y = "Proportion of gene clusters",
       fill = "COG Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,-5, "pt"),
        #axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank())
dev.off()

# Or plot as panels to see the distribution of each one
pdf("InitialFigs/GeneClustersProp_NumGenomes_Facet.pdf", width = 9, height = 6)
ggplot(gcProp_genomeNum, aes(num_genomes_gene_cluster_has_hits, Prop, fill = COG20_Uniq)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of genomes (1 to 181)",
       y = "Proportion of gene clusters",
       fill = "COG Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_wrap(~ COG20_Uniq, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 4))
dev.off()

# Or plot as heatmap
gcProp_genomeNum_wide <- gcProp_genomeNum %>%
  dplyr::select(-nGC, -sum) %>%
  pivot_wider(names_from = COG20_Uniq, values_from = Prop) %>%
  column_to_rownames(var = "num_genomes_gene_cluster_has_hits") %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %notin% c("Unknown", "Multiple")) %>%
  replace(is.na(.), 0)
pheatmap(gcProp_genomeNum_wide,
         color = c("white", colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100)),
         breaks = c(0, 
                    seq(2.393203e-05, max(gcProp_genomeNum_wide), 
                        length.out = 100)),
         cluster_rows = T,
         cluster_cols = F,
         #scale = "column",
         fontsize = 8,
         show_colnames = F)
gcProp_genomeNum_wide <- gcProp_genomeNum %>%
  dplyr::select(-nGC, -sum) %>%
  pivot_wider(names_from = COG20_Uniq, values_from = Prop) %>%
  column_to_rownames(var = "num_genomes_gene_cluster_has_hits") %>%
  t() %>%
  as.data.frame()
pheatmap(gcProp_genomeNum_wide,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         fontsize = 8,
         show_colnames = F)

#### __Per COG ####
# Need to normalize by the total number of GC in each COG
# Also need to double count those that are multiple
# Need to do COG by COG
# Number of GC in each COG in each NumGenomes
gc <- read.delim("~/Desktop/Fierer/Strains/Australia/Bradyrhizobium_Pan_gene_clusters_summary_181.txt", quote = "") %>%
  dplyr::select(gene_cluster_id, num_genomes_gene_cluster_has_hits,
                num_genes_in_gene_cluster, functional_homogeneity_index,
                geometric_homogeneity_index, COG20_CATEGORY_ACC, COG20_CATEGORY,
                COG20_FUNCTION) %>%
  arrange(gene_cluster_id, COG20_CATEGORY)
unique(gc$COG20_CATEGORY_ACC)
cogB <- gc %>%
  filter(grepl("B", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Chromatin structure and dynamics")
cogC <- gc %>%
  filter(grepl("C", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Energy Production and Conversion")
cogD <- gc %>%
  filter(grepl("D", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Cell cycle control, cell division, chromosome partitioning")
cogE <- gc %>%
  filter(grepl("E", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Amino acid transport and metabolism")
cogF <- gc %>%
  filter(grepl("F", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Nucleotide transport and metabolism")
cogG <- gc %>%
  filter(grepl("G", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Carbohydrate transport and metabolism")
cogH <- gc %>%
  filter(grepl("H", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Coenzyme transport and metabolism")
cogI <- gc %>%
  filter(grepl("I", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Lipid transport and metabolism")
cogJ <- gc %>%
  filter(grepl("J", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Translation, ribosomal structure and biogenesis")
cogK <- gc %>%
  filter(grepl("K", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Transcription")
cogL <- gc %>%
  filter(grepl("L", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Replication, recombination and repair")
cogM <- gc %>%
  filter(grepl("M", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Cell wall/membrane/envelope biogenesis")
cogN <- gc %>%
  filter(grepl("N", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Cell motility")
cogO <- gc %>%
  filter(grepl("O", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Posttranslational modification, protein turnover, chaperones")
cogP <- gc %>%
  filter(grepl("P", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Inorganic ion transport and metabolism")
cogQ <- gc %>%
  filter(grepl("Q", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Secondary metabolites biosynthesis, transport and catabolism")
cogR <- gc %>%
  filter(grepl("R", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "General function prediction only")
cogS <- gc %>%
  filter(grepl("S", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Function unknown")
cogT <- gc %>%
  filter(grepl("T", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Signal transduction mechanisms")
cogU <- gc %>%
  filter(grepl("U", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Intracellular trafficking, secretion, and vesicular transport")
cogV <- gc %>%
  filter(grepl("V", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Defense mechanisms")
cogW <- gc %>%
  filter(grepl("W", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Extracellular structures")
cogX <- gc %>%
  filter(grepl("X", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Mobilome: prophages, transposons")
cogZ <- gc %>%
  filter(grepl("Z", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "Cytoskeleton")
cogNA <- gc %>%
  filter(COG20_CATEGORY_ACC == "") %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.)) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(Count = n(),
            TotalCOG = mean(TotalCOG)) %>%
  ungroup() %>%
  mutate(Prop = Count/TotalCOG,
         COG = "NA")
cogProp <- rbind(cogB, cogC, cogD, cogE, cogF, cogG, cogH, cogI, cogJ, cogK,
                 cogL, cogM, cogN, cogO, cogP, cogQ, cogR, cogS, cogT, cogU,
                 cogV, cogW, cogX, cogZ, cogNA)
length(unique(cogProp$COG))
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
pdf("InitialFigs/GeneClustersPropCOG.pdf", width = 8, height = 6)
ggplot(cogProp, aes(num_genomes_gene_cluster_has_hits, Prop)) +
  geom_bar(stat = "identity", colour = "black", fill = "black") +
  labs(x = "Number of genomes (1 to 181)",
       y = "Proportion of gene clusters in COG",
       fill = "COG Category") +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~ COG, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 4))
dev.off()

cogProp_sing_core <- cogProp %>%
  filter(num_genomes_gene_cluster_has_hits == 1 | num_genomes_gene_cluster_has_hits == 181) %>%
  mutate(Type = ifelse(num_genomes_gene_cluster_has_hits == 1,
                       "Singleton", "Core")) %>%
  dplyr::select(COG, Prop, Type) %>%
  add_row(COG = "Chromatin structure and dynamics", Prop = 0, Type = "Core") %>%
  add_row(COG = "Cytoskeleton", Prop = 0, Type = "Core") %>%
  mutate(Type = factor(Type, levels = c("Singleton", "Core")))
pdf("InitialFigs/GeneClustersPropCOG_SingCore.pdf", width = 8, height = 6)
ggplot(cogProp_sing_core, aes(COG, Prop, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "COG Category",
       y = "Proportion of gene clusters in COG",
       fill = "Type") +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 5, 5, 60, "pt"))
dev.off()

cogProp_Type <- cogProp %>%
  mutate(Type = ifelse(num_genomes_gene_cluster_has_hits <= 9,
                       "Low",
                       ifelse(num_genomes_gene_cluster_has_hits > 9 & num_genomes_gene_cluster_has_hits < 172, "Medium", "High"))) %>%
  dplyr::select(COG, Prop, Type, TotalCOG) %>%
  add_row(COG = "Chromatin structure and dynamics", Prop = 0, Type = "Medium", TotalCOG = 1) %>%
  add_row(COG = "Chromatin structure and dynamics", Prop = 0, Type = "High", TotalCOG = 1) %>%
  add_row(COG = "Cytoskeleton", Prop = 0, Type = "Medium", TotalCOG = 4) %>%
  add_row(COG = "Cytoskeleton", Prop = 0, Type = "High", TotalCOG = 4) %>%
  mutate(Type = factor(Type, levels = c("Low", "Medium", "High"))) %>%
  group_by(COG, Type) %>%
  summarise(Prop = sum(Prop),
            nGC = mean(TotalCOG)) %>%
  ungroup()
num_lab <- cogProp_Type %>%
  group_by(COG) %>%
  summarise(y = max(Prop),
            nGC = mean(nGC))
ggplot(cogProp_Type, aes(COG, Prop, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = num_lab, aes(COG, y + 0.05, label = nGC),
            inherit.aes = F, size = 2) +
  labs(x = "COG Category",
       y = "Proportion of gene clusters in COG",
       fill = "Type") +
  scale_fill_manual(values = c("#DDAA33", "#BB5566", "#004488"),
                    labels = c("Low (1-9)", "Medium (10-171)", "High (172-181)")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 5, 5, 60, "pt"))

cogTot <- cogProp %>%
  group_by(COG) %>%
  summarise(nGC = mean(TotalCOG),
            y = max(Prop)) %>%
  ungroup()
ggplot(cogTot, aes(COG, nGC)) +
  geom_bar(stat = "identity") +
  labs(x = "COG Category",
       y = "# gene clusters in COG") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 60, "pt"))



#### ___Figure 7 ####
# For each COG, get mean Geometric Homogeneity, mean prop. genomes present, mean prop of GC
cogB <- gc %>%
  filter(grepl("B", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Chromatin structure and dynamics")
cogC <- gc %>%
  filter(grepl("C", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Energy Production and Conversion")
cogD <- gc %>%
  filter(grepl("D", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Cell cycle control, cell division, chromosome partitioning")
cogE <- gc %>%
  filter(grepl("E", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Amino acid transport and metabolism")
cogF <- gc %>%
  filter(grepl("F", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Nucleotide transport and metabolism")
cogG <- gc %>%
  filter(grepl("G", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Carbohydrate transport and metabolism")
cogH <- gc %>%
  filter(grepl("H", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Coenzyme transport and metabolism")
cogI <- gc %>%
  filter(grepl("I", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Lipid transport and metabolism")
cogJ <- gc %>%
  filter(grepl("J", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Translation, ribosomal structure and biogenesis")
cogK <- gc %>%
  filter(grepl("K", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Transcription")
cogL <- gc %>%
  filter(grepl("L", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Replication, recombination and repair")
cogM <- gc %>%
  filter(grepl("M", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Cell wall/membrane/envelope biogenesis")
cogN <- gc %>%
  filter(grepl("N", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Cell motility")
cogO <- gc %>%
  filter(grepl("O", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Posttranslational modification, protein turnover, chaperones")
cogP <- gc %>%
  filter(grepl("P", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Inorganic ion transport and metabolism")
cogQ <- gc %>%
  filter(grepl("Q", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Secondary metabolites biosynthesis, transport and catabolism")
cogR <- gc %>%
  filter(grepl("R", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "General function prediction only")
cogS <- gc %>%
  filter(grepl("S", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Function unknown")
cogT <- gc %>%
  filter(grepl("T", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Signal transduction mechanisms")
cogU <- gc %>%
  filter(grepl("U", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Intracellular trafficking, secretion, and vesicular transport")
cogV <- gc %>%
  filter(grepl("V", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Defense mechanisms")
cogW <- gc %>%
  filter(grepl("W", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Extracellular structures")
cogX <- gc %>%
  filter(grepl("X", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Mobilome: prophages, transposons")
cogZ <- gc %>%
  filter(grepl("Z", COG20_CATEGORY_ACC)) %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "Cytoskeleton")
cogNA <- gc %>%
  filter(COG20_CATEGORY_ACC == "") %>%
  group_by(gene_cluster_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(TotalCOG = nrow(.),
         PA = num_genomes_gene_cluster_has_hits/181,
         COG = "NA")

cogProp <- rbind(cogB, cogC, cogD, cogE, cogF, cogG, cogH, cogI, cogJ, cogK,
                 cogL, cogM, cogN, cogO, cogP, cogQ, cogR, cogS, cogT, cogU,
                 cogV, cogW, cogX, cogZ, cogNA)

cogProp_Type <- cogProp %>%
  mutate(Type = ifelse(num_genomes_gene_cluster_has_hits <= 9,
                       "Low",
                       ifelse(num_genomes_gene_cluster_has_hits > 9 & num_genomes_gene_cluster_has_hits < 172, "Medium", "High"))) %>%
  dplyr::select(COG, Type, geometric_homogeneity_index, PA, TotalCOG) %>%
  mutate(Type = factor(Type, levels = c("Low", "Medium", "High"))) %>%
  group_by(COG, Type) %>%
  summarise(Count = n(),
            nGC = mean(TotalCOG),
            GH = mean(geometric_homogeneity_index),
            PA = mean(PA)) %>%
  ungroup() %>%
  mutate(Prop = Count/nGC) %>%
  add_row(COG = "Chromatin structure and dynamics", Type = "Medium", Count = NA,
          nGC = 1, GH = NA, PA = NA, Prop = NA) %>%
  add_row(COG = "Chromatin structure and dynamics", Type = "High", Count = NA,
          nGC = 1, GH = NA, PA = NA, Prop = NA) %>%
  add_row(COG = "Cytoskeleton", Type = "Medium", Count = NA,
          nGC = 4, GH = NA, PA = NA, Prop = NA) %>%
  add_row(COG = "Cytoskeleton", Type = "High", Count = NA,
          nGC = 4, GH = NA, PA = NA, Prop = NA) %>%
  dplyr::select(-Count, -nGC) %>%
  mutate(COG = as.factor(COG),
         Type = as.factor(Type))
cogProp_Type_wide <- pivot_wider(cogProp_Type, 
                                 names_from = Type, 
                                 values_from = c(GH, PA, Prop), 
                                 names_sep = "_") %>%
  column_to_rownames(var = "COG")
colSums(cogProp_Type_wide, na.rm = T)
min(cogProp_Type_wide, na.rm = T)
max(cogProp_Type_wide, na.rm = T)
ann_cols <- data.frame(row.names = colnames(cogProp_Type_wide),
                       "Metric" = c("Geometric", "Geometric", "Geometric",
                                    "PA", "PA", "PA",
                                    "Prop", "Prop", "Prop"),
                       "Genomes" = c("Low (1-9)", "Medium (10-171)", "High (172-181)",
                                     "Low (1-9)", "Medium (10-171)", "High (172-181)",
                                     "Low (1-9)", "Medium (10-171)", "High (172-181)"))
ann_colors <- list(Metric = c(Geometric = "#440154FF",
                              PA = "#21908CFF",
                              Prop = "#FDE725FF"),
                   Genomes = c("Low (1-9)" = "#DDAA33",
                               "Medium (10-171)" = "#BB5566",
                               "High (172-181)" = "#004488"))
pheatmap(cogProp_Type_wide,
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(3, 6),
         show_colnames = F,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         angle_col = 315,
         display_numbers = T,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6,
         border_color = "white",
         filename = "FinalFigs/Figure7.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Or facet wrapped barplot instead of heatmap
cogProp_Type_long <- cogProp_Type %>%
  pivot_longer(cols = c(GH, PA, Prop)) %>%
  mutate(Type = factor(Type, levels = c("Low", "Medium", "High")))
# Remove presence absence panel, NA, cytoskeleton, and chromatin
cogProp_Type_long <- cogProp_Type_long %>%
  filter(COG %notin% c("Chromatin structure and dynamics",
                       "Cytoskeleton")) %>%
  filter(name != "PA") %>%
  droplevels()
facet_names <- c("GH" = "Mean Geometric Homogeneity Index",
                 #"PA" = "Mean Genome Proportion Presence/Absence",
                 "Prop" = "Proportion of Gene Clusters in COG")
num_lab2 <- num_lab %>%
  mutate(name = "Prop") %>%
  filter(COG %notin% c("Chromatin structure and dynamics",
                       "Cytoskeleton")) %>%
  droplevels()
fig7 <- ggplot() +
  geom_bar(data = subset(cogProp_Type_long, name %in% c("PA", "Prop")),
           aes(COG, value, fill = Type),
           stat = "identity", position = "dodge") +
  geom_point(data = subset(cogProp_Type_long, name %in% c("GH")),
             aes(COG, value, colour = Type),
             position = position_dodge(width = 0.7), size = 2, show.legend = F) +
  geom_text(data = num_lab2, aes(COG, y + 0.05, label = nGC),
            inherit.aes = F, size = 2) +
  labs(x = "COG Category",
       y = "Proportion",
       fill = "Genomes") +
  scale_colour_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  scale_fill_manual(values = c("#DDAA33", "#BB5566", "#004488"),
                    labels = c("Low (1-9)", "Medium (10-171)", "High (172-181)")) +
  facet_wrap(~ name, ncol = 1, scales = "free_y", labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme_bw() +
  theme(legend.position = "left",
        legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.margin = margin(0, -5, 0, 0, "pt"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title = element_blank(),
        plot.margin = margin(5, 5, 5, 10, "pt"))
fig7
pdf("FinalFigs/Figure7.pdf", width = 8, height = 6)
fig7
dev.off()
png("FinalFigs/Figure7.png", width = 8, height = 6, units = "in", res = 300)
fig7
dev.off()



#### 14. Plants ####
d_331_veg <- d_331 %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(SampleID, vegetation_dom_grasses, vegetation_dom_shrubs, 
                vegetation_dom_trees, vegetation_total_cover, vegetation_type) %>%
  filter(vegetation_dom_grasses != "")
d_species_uniq <- d_331_veg %>%
  separate_rows(vegetation_dom_grasses, sep = "; ") %>%
  separate_rows(vegetation_dom_shrubs, sep = "; ") %>%
  separate_rows(vegetation_dom_trees, sep = "; ") %>%
  dplyr::select(SampleID, vegetation_dom_grasses, vegetation_dom_shrubs, 
                vegetation_dom_trees) %>%
  filter(vegetation_dom_grasses != "") %>%
  pivot_longer(c(vegetation_dom_grasses, vegetation_dom_shrubs, 
                 vegetation_dom_trees), names_to = "Group", values_to = "Species") %>%
  separate(Species, sep = " ", into = c("Genus", "Species", "Abund")) %>%
  mutate(Species = tolower(as.character(Species))) %>%
  mutate(Spp = paste(Genus, Species, sep = " ")) %>%
  mutate(Genus = gsub("Acasia", "Acacia", Genus)) %>%
  mutate(Genus = gsub("Cloris", "Chloris", Genus)) %>%
  mutate(Genus = gsub("Cymbogon", "Cymbopogon", Genus)) %>%
  mutate(Genus = gsub("lomandra", "Lomandra", Genus)) %>%
  mutate(Genus = gsub("Noteleaea", "Notelaea", Genus)) %>%
  mutate(Spp = gsub("Cloris ventricosa", "Chloris ventricosa", Spp)) %>%
  mutate(Spp = gsub("Cymbogon refractus", "Cymbopogon refractus", Spp)) %>%
  mutate(Spp = gsub("Lepidosperma lateral", "Lepidosperma laterale", Spp)) %>%
  mutate(Spp = gsub("Lepidosperma lateralee", "Lepidosperma laterale", Spp)) %>%
  mutate(Spp = gsub("Noteleaea longofolia", "Notelaea longifolia", Spp)) %>%
  mutate(Spp = gsub("lomandra filiformis", "Lomandra filiformis", Spp)) %>%
  mutate(Spp = gsub("Eucalyptus cebra", "Eucalyptus crebra", Spp)) %>%
  mutate(Spp = gsub("Acasia parramattensis", "Acacia parramattensis", Spp)) %>%
  mutate(Spp = gsub("Corymbia gummifeia", "Corymbia gummifera", Spp)) %>%
  mutate(Spp = gsub("Elaeodendron austale", "Elaeodendron australe", Spp)) %>%
  mutate(Spp = gsub("Leptospermum trinenium", "Leptospermum trinervium", Spp)) %>%
  mutate(Spp = gsub("Lomandra logifolia", "Lomandra longifolia", Spp)) %>%
  mutate(Spp = gsub("Oplismenus imbecillus", "Oplismenus imbecillis", Spp)) %>%
  mutate(Spp = gsub("Syzygium austale", "Syzygium australe", Spp)) %>%
  mutate(Spp = gsub("Pittosporum revolutum\\(3%\\)", "Pittosporum revolutum", Spp)) %>%
  mutate(Spp = gsub("Rytidosperma \\(1%\\)", "Rytidosperma sp.", Spp)) %>%
  mutate(Spp = gsub("Rytidosperma \\(5%\\)", "Rytidosperma sp.", Spp)) %>%
  group_by(Spp) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  filter(Spp != "0 NA") %>%
  arrange(Spp)
# Need to get family level because first pass showed species and genus not working well!
d_species_uniq$Family <- c(rep("Fabaceae", 6), 
                           rep("Casuarinaceae", 2),
                           "Rhamnaceae",
                           rep("Myrtaceae", 4),
                           rep("Poaceae", 3),
                           "Rubiaceae",
                           rep("Poaceae", 2),
                           rep("Proteaceae", 3),
                           "Poaceae",
                           rep("Malvaceae", 2),
                           "Phyllanthaceae",
                           "Pittosporaceae",
                           "Myrtaceae",
                           rep("Cunoneaceae", 2),
                           rep("Poaceae", 2),
                           rep("Lamiaceae", 2),
                           "Myrtaceae",
                           "Cyperaceae",
                           "Poaceae",
                           "Poaceae",
                           "Cyperaceae",
                           rep("Fabaceae", 2),
                           "Celastraceae",
                           rep("Asphodelaceae", 3),
                           "Poaceae",
                           "Fabaceae",
                           "Sapindaceae",
                           "Elaeocarpaceae",
                           "Celastraceae",
                           rep("Poaceae", 2),
                           "Scrophulariaceae",
                           rep("Myrtaceae", 13),
                           "Santalaceae",
                           "Moraceae",
                           "Phyllanthaceae",
                           "Proteaceae",
                           rep("Proteaceae", 2),
                           "Dilleniaceae",
                           "Euphorbiaceae",
                           "Pittosporaceae",
                           "Poaceae",
                           "Fabaceae",
                           "Proteaceae",
                           "Juncaceae",
                           rep("Myrtaceae", 2),
                           "Poaceae",
                           "Cyperaceae",
                           rep("Myrtaceae", 2),
                           rep("Asparagaceae", 6),
                           "Celastraceae",
                           rep("Myrtaceae", 4),
                           "Meliaceae",
                           "Poaceae",
                           "Scrophulariaceae",
                           rep("Oleaceae", 2),
                           rep("Poaceae", 3),
                           "Asteraceae",
                           rep("Poaceae", 2),
                           "Poaceae",
                           "Proteaceae",
                           "Phyllanthaceae",
                           rep("Pittosporaceae", 2),
                           "Araliaceae",
                           "Cyperaceae",
                           "Fabaceae",
                           "Poaceae",
                           "Poaceae",
                           "Myrtaceae",
                           "Myrtaceae",
                           rep("Poaceae", 2),
                           "Cannabaceae")
unique(d_species_uniq$Spp) # 122 species across the 48/331 samples w veg data
unique(d_species_uniq$Genus) # 71 genera
unique(d_species_uniq$Family) # 29 families

d_species <- d_331_veg %>%
  separate_rows(vegetation_dom_grasses, sep = "; ") %>%
  separate_rows(vegetation_dom_shrubs, sep = "; ") %>%
  separate_rows(vegetation_dom_trees, sep = "; ") %>%
  dplyr::select(SampleID, vegetation_dom_grasses, vegetation_dom_shrubs, 
                vegetation_dom_trees) %>%
  filter(vegetation_dom_grasses != "") %>%
  pivot_longer(c(vegetation_dom_grasses, vegetation_dom_shrubs, 
                 vegetation_dom_trees), names_to = "Group", values_to = "Species") %>%
  separate(Species, sep = " ", into = c("Genus", "Species", "Abund")) %>%
  mutate(Species = tolower(as.character(Species))) %>%
  mutate(Spp = paste(Genus, Species, sep = " ")) %>%
  mutate(Genus = gsub("Acasia", "Acacia", Genus)) %>%
  mutate(Genus = gsub("Cloris", "Chloris", Genus)) %>%
  mutate(Genus = gsub("Cymbogon", "Cymbopogon", Genus)) %>%
  mutate(Genus = gsub("lomandra", "Lomandra", Genus)) %>%
  mutate(Genus = gsub("Noteleaea", "Notelaea", Genus)) %>%
  mutate(Spp = gsub("Cloris ventricosa", "Chloris ventricosa", Spp)) %>%
  mutate(Spp = gsub("Cymbogon refractus", "Cymbopogon refractus", Spp)) %>%
  mutate(Spp = gsub("Lepidosperma lateral", "Lepidosperma laterale", Spp)) %>%
  mutate(Spp = gsub("Lepidosperma lateralee", "Lepidosperma laterale", Spp)) %>%
  mutate(Spp = gsub("Noteleaea longofolia", "Notelaea longifolia", Spp)) %>%
  mutate(Spp = gsub("lomandra filiformis", "Lomandra filiformis", Spp)) %>%
  mutate(Spp = gsub("Eucalyptus cebra", "Eucalyptus crebra", Spp)) %>%
  mutate(Spp = gsub("Acasia parramattensis", "Acacia parramattensis", Spp)) %>%
  mutate(Spp = gsub("Corymbia gummifeia", "Corymbia gummifera", Spp)) %>%
  mutate(Spp = gsub("Elaeodendron austale", "Elaeodendron australe", Spp)) %>%
  mutate(Spp = gsub("Leptospermum trinenium", "Leptospermum trinervium", Spp)) %>%
  mutate(Spp = gsub("Lomandra logifolia", "Lomandra longifolia", Spp)) %>%
  mutate(Spp = gsub("Oplismenus imbecillus", "Oplismenus imbecillis", Spp)) %>%
  mutate(Spp = gsub("Syzygium austale", "Syzygium australe", Spp)) %>%
  mutate(Spp = gsub("Pittosporum revolutum\\(3%\\)", "Pittosporum revolutum", Spp)) %>%
  mutate(Spp = gsub("Rytidosperma \\(1%\\)", "Rytidosperma sp.", Spp)) %>%
  mutate(Spp = gsub("Rytidosperma \\(5%\\)", "Rytidosperma sp.", Spp)) %>%
  filter(Spp != "0 NA") %>%
  group_by(SampleID, Spp) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Spp)



#### _Chloroplasts ####
# Try Sylph on chloroplast database

# Note that BASE has some plant data for rows 113-159 (n = 46), can use to validate
# Download the metadata https://ngdc.cncb.ac.cn/cgir/genome (Hua et al. 2022)
# plant_meta <- read.csv("data/Genome_Infor_Result.csv")
# Nope, can only do 20 at a time. Use NCBI
plant_meta <- read.delim("data/ncbi_chloroplast_metadata.tsv") %>% # 14707
  separate(Organism.name, remove = F, into = c("Genus", "Species"), sep = " ")
length(unique(plant_meta$Organism.name)) # 14687
length(unique(plant_meta$Organism.taxid)) # 14687
table(plant_meta$Organelle.type) # 12989 chloro, 1718 plastid
table(plant_meta$Topology) # 14666 circular, 40 linear, 1 unknown
range(plant_meta$Length) # 11348 to 1352306
hist(plant_meta$Length)

# Need to see if the species in the veg data have chloroplast genomes
sum(d_species_uniq$Spp %in% plant_meta$Organism.name) # 13
sum(unique(d_species_uniq$Genus) %in% unique(plant_meta$Genus)) # 39
d_species_db <- d_species_uniq %>%
  filter(Spp %in% plant_meta$Organism.name)
d_species_uniq$match <- sapply(d_species_uniq$Spp, function(x) any(grepl(x, plant_meta$Organism.name)))
sum(d_species_uniq$match == TRUE) # Only 13 species were in the chloroplast database



# 99%
chloro_sylph_raw <- read.delim("data/sylph_profile_chloro14706_99.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file))
chloro_sylph <- read.delim("data/sylph_profile_chloro14706_99.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  left_join(., d_331_veg, by = "SampleID")
length(unique(chloro_sylph$SampleID)) # 18/331 samples
length(unique(chloro_sylph$GenomeID)) # 18/14606 genomes

# 98%
chloro_sylph_raw <- read.delim("data/sylph_profile_chloro14706_98.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file))
chloro_sylph <- read.delim("data/sylph_profile_chloro14706_98.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  left_join(., d_331_veg, by = "SampleID")
length(unique(chloro_sylph$SampleID)) # 42/331 samples
length(unique(chloro_sylph$GenomeID)) # 41/14606 genomes

# 97%
chloro_sylph_raw <- read.delim("data/sylph_profile_chloro14706_97.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file))
chloro_sylph <- read.delim("data/sylph_profile_chloro14706_97.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  left_join(., d_331_veg, by = "SampleID")
length(unique(chloro_sylph$SampleID)) # 55/331 samples
length(unique(chloro_sylph$GenomeID)) # 50/14606 genomes

# 96%
chloro_sylph_raw <- read.delim("data/sylph_profile_chloro14706_96.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file))
chloro_sylph <- read.delim("data/sylph_profile_chloro14706_96.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  left_join(., d_331_veg, by = "SampleID")
length(unique(chloro_sylph$SampleID)) # 62/331 samples
length(unique(chloro_sylph$GenomeID)) # 56/14606 genomes

# 95%
chloro_sylph_raw <- read.delim("data/sylph_profile_chloro14706_95.tsv") %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file))
chloro_sylph <- read.delim("data/sylph_profile_chloro14706_95.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/Australia_copy/QC_reads/", 
                         "", Sample_file)) %>%
  separate(SampleID, into = c("SampleID", "Junk1", "Junk2"), sep = "_") %>%
  dplyr::select(-Junk1, -Junk2) %>%
  mutate(GenomeID = gsub("/scratch/alpine/clbd1748/Australia_copy/chloroplasts/",
                         "", Genome_file)) %>%
  mutate(GenomeID = gsub(".fna",
                         "", GenomeID)) %>%
  left_join(., d_331_veg, by = "SampleID") %>%
  left_join(., plant_meta, by = c("GenomeID" = "RefSeq.accession"))
length(unique(chloro_sylph$SampleID)) # 67/331 samples
length(unique(chloro_sylph$GenomeID)) # 67/14606 genomes
length(unique(chloro_sylph$Species)) # 64 species
length(unique(chloro_sylph$Genus)) # 60 genera
# Note there was one sample with veg data
chloro_sylph_wVeg <- chloro_sylph %>%
  filter(SampleID %in% d_331_veg$SampleID)
  
# That had Corymbia gummifera (8%) which is in the db, and Sylph didn't detect it.
rich <- chloro_sylph %>%
  group_by(SampleID) %>%
  count()
table(rich$n)

# How many of these 67 were in the BASE veg data?
sum(unique(chloro_sylph$Organism.name) %in% d_species$Spp) # 1
sum(unique(chloro_sylph$Genus) %in% unique(d_species_uniq$Genus)) # 5
chloro_sylph_match <- chloro_sylph %>%
  filter(chloro_sylph$Organism.name %in% d_species$Spp)



#### _trnL ####
# Lucky Devin Leopold at Jonah Ventures already had made a trnL database!
trnL <- microseq::readFasta("~/Desktop/Fierer/Strains/trnL_taxDB_20250217.fasta")
trnL2 <- trnL %>%
  separate(Header, into = c("IDENTIFIER", "Domain", "Phylum", "Class", "Order",
                            "Family", "Genus", "Species"), sep = ";")
length(unique(trnL2$Species)) # 65208
# So there are some duplicates. Check
trnL_dup <- trnL2 %>%
  group_by(Species) %>%
  summarise(count = n())
# Good to know, but I think duplicates are okay
# Check the length distribution. Then format for phyloFlash
# Range is 100 to 289, don't filter
trnL_spp <- trnL2 %>%
  mutate(INTEGER1 = 1,
         INTEGER2 = nchar(Sequence)) %>%
  mutate(TAXONOMY = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";")) %>%
  mutate(Header = paste(IDENTIFIER, ".", INTEGER1, ".", INTEGER2, " ", TAXONOMY, 
                        sep = "")) %>%
  dplyr::select(Header, Sequence)
range(trnL_spp$INTEGER2) # 100 to 289
hist(trnL_spp$INTEGER2) # Most are around 150 bp
#microseq::writeFasta(trnL_spp, "~/Desktop/Fierer/Strains/SILVA_CustomDB_trnL.fasta")

# Need to see if the species and genera in the veg data have trnL seqs
sum(d_species_uniq$Spp %in% trnL2$Species) # 47/122 (this is 34 more than chloroplast)
sum(unique(d_species_uniq$Genus) %in% unique(trnL2$Genus)) # 65/71
sum(unique(d_species_uniq$Family) %in% unique(trnL2$Family)) # 28/29
d_species_db <- d_species_uniq %>%
  filter(Spp %in% trnL2$Species)
d_species_uniq$match <- sapply(d_species_uniq$Spp, function(x) any(grepl(x, trnL2$Species)))
sum(d_species_uniq$match == TRUE)
# 48 (i.e, using grep found 1 more than %in%)


# Analyze phyloFlash results
setwd("~/Desktop/Fierer/Strains/Australia/pf_trnL/")
# Some output files had zero bytes so removed
files <- list.files()
files
# 304/331 samples have phyloFlash data! 27 had no plants detected!

merged_data <- files %>%
  lapply(function(file) {
    sample_id <- gsub("\\.phyloFlash.NTUfull_abundance\\.csv", "", basename(file))
    data <- read.csv(file, header = FALSE) %>%
      set_names(c("NTU", sample_id)) %>%
      return(data)
  }) %>%
  reduce(function(x, y) merge(x, y, by = "NTU", all = TRUE))
merged_data[is.na(merged_data)] <- 0
sample_id <- data.frame(fullname = list.files()) %>%
  separate(fullname, into = c("sampleID", "Junk1", "Junk2", sep = "_"))
names(merged_data) <- c("NTU", sample_id$sampleID)
setwd("~/Documents/GitHub/AussieStrains/")

# Now see if those 47 species, 65 genera and 28 families with trnL seqs were detected
# In the whole dataset:
veg_pf <- merged_data %>%
  mutate(Tot = rowSums(.[2:ncol(.)])) %>%
  arrange(desc(Tot)) %>%
  separate(NTU, into = c("Domain", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"), sep = ";", remove = F) %>%
  mutate(Spp = paste(Genus, Species, sep = " "))
sum(d_species_uniq$Spp %in% unique(veg_pf$Species)) # 6 of 47 possible identified
sum(unique(d_species_uniq$Genus) %in% unique(veg_pf$Genus)) # 34  of 65 possible identified
sum(unique(d_species_uniq$Family) %in% unique(veg_pf$Family)) # 27 of 28 possible identified

# In the 36/48 samples with veg and pf data
Fam <- d_species_uniq %>%
  dplyr::select(Spp, Family)
d_species_wPF <- d_species %>%
  filter(SampleID %in% names(merged_data)) %>%
  left_join(., Fam, by = "Spp")
d_species_wPF_uniq <- d_species_uniq %>%
  filter(SampleID %in% names(merged_data))
length(unique(d_species$SampleID)) # 48
length(unique(d_species_wPF$SampleID)) # 36 overlap
length(unique(d_species_wPF_uniq$SampleID)) # Only 26, don't use
length(unique(d_species_wPF$Spp)) # 110 of the 122 in the overlapping dataset
length(unique(d_species_wPF$Genus)) # 63 of the 71 in the overlapping dataset
length(unique(d_species_wPF$Family)) # 24 of the 29 families in the overlapping dataset
sum(unique(d_species_wPF$Spp) %in% unique(trnL2$Species)) # 41
sum(unique(d_species_wPF$Genus) %in% unique(trnL2$Genus)) # 57
sum(unique(d_species_wPF$Family) %in% unique(trnL2$Family)) # 23
veg_pf <- merged_data %>%
  dplyr::select(NTU, all_of(unique(d_species_wPF$SampleID))) %>%
  mutate(Tot = rowSums(.[2:ncol(.)])) %>%
  arrange(desc(Tot)) %>%
  separate(NTU, into = c("Domain", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"), sep = ";", remove = F) %>%
  mutate(Spp = paste(Genus, Species, sep = " "))
sum(unique(d_species_wPF$Spp) %in% unique(veg_pf$Species)) # 6 of 41 possible identified
sum(unique(d_species_wPF$Genus) %in% unique(veg_pf$Genus)) # 28  of 57 possible identified
sum(unique(d_species_wPF$Family) %in% unique(veg_pf$Family)) # 22 of 23 possible identified

# Andrew sent over more plant data
ala <- read.csv("data/records_plant_climate_ALA.csv")
# This doesn't have any species though, just more broad categories and info

# Check num seqs per sample
seqs <- data.frame(SampleID = names(merged_data[2:ncol(merged_data)]),
                   Seqs = colSums(merged_data[2:ncol(merged_data)]))
hist(seqs$Seqs)
range(seqs$Seqs)

# Convert the merged data to mctoolrs format
# Change parentheses to NA
length(unique(merged_data$NTU))
seqtab_wTax <- merged_data %>%
  mutate(RowID = row_number()) %>%
  mutate(OTU = paste("OTU", RowID, sep = "_")) %>%
  separate(NTU, into = c("Domain", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"), sep = ";") %>%
  mutate(across(everything(), ~ ifelse(grepl("\\(|\\)", .), "NA", .))) %>%
  mutate(taxonomy = paste(Domain, Phylum, Class, Order, Family, Genus, Species,
                          sep = ";")) %>%
  dplyr::select(-RowID) %>%
  dplyr::select(OTU, everything(), taxonomy)
names(seqtab_wTax)
#out_fp <- "data/seqtab_wTax_mctoolsr_trnL.txt"
#write("#Exported for mctoolsr", out_fp)
#suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Import
tax_table_fp <- "data/seqtab_wTax_mctoolsr_trnL.txt"
map_fp <- "~/Documents/GitHub/AussieStrains/data/metadata_331.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 304 samples loaded
input$map_loaded$sampleID <- rownames(input$map_loaded)
input <- filter_data(input, 
                     filter_cat = "sampleID", 
                     keep_vals = d_sylph$sampleID) # 244 of 268 have plant data
# Check singletons and doubletons
singdoub <- data.frame("count" = rowSums(input$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 3665. Don't remove, since will aggregate
input_filt <- filter_taxa_from_input(input,
                                     taxa_IDs_to_remove = singdoub$ASV) # 3665 removed

# Rarefaction
rarecurve(t(input$data_loaded), step = 100, label = F)
rarecurve(t(input_filt$data_loaded), step = 100, label = F)

# That is very unrealistic. Need to go ahead and aggregate to species, genus, family etc.
tax_sum_phyla <- summarize_taxonomy(input = input, 
                                    level = 2, 
                                    report_higher_tax = F) # 1, Streptophyta
tax_sum_class <- summarize_taxonomy(input = input, 
                                    level = 3, 
                                    report_higher_tax = F) # 20
tax_sum_order <- summarize_taxonomy(input = input, 
                                    level = 4, 
                                    report_higher_tax = F) # 120
tax_sum_family <- summarize_taxonomy(input = input, 
                                     level = 5, 
                                     report_higher_tax = F) # 477
tax_sum_family2 <- summarize_taxonomy(input = input,
                                      relative = F,
                                      level = 5, 
                                      report_higher_tax = F) # 477
tax_sum_genus <- summarize_taxonomy(input = input, 
                                     level = 6, 
                                     report_higher_tax = F) # 2926
tax_sum_species <- summarize_taxonomy(input = input, 
                                      level = 7, 
                                      report_higher_tax = F) # 6421
tax_sum_species2 <- summarize_taxonomy(input = input, 
                                       relative = F,
                                       level = 7, 
                                       report_higher_tax = F) # 6421
pdf("InitialFigs/trnL_rarefaction.pdf", width = 7, height = 5)
rarecurve(t(tax_sum_species2), step = 100, label = F)
dev.off()

sort(colSums(input$data_loaded))
pdf("InitialFigs/trnL_seqs.pdf", width = 7, height = 5)
plot(sort(colSums(input$data_loaded)), xlab = "Sample rank", ylab = "# seqs")
dev.off()
mean(colSums(input$data_loaded)) # 555.6434
se(colSums(input$data_loaded)) # 49.8105

sort(colSums(input_filt$data_loaded))
plot(sort(colSums(input_filt$data_loaded)))
mean(colSums(input_filt$data_loaded)) # 536.0287
se(colSums(input_filt$data_loaded)) # 48.09166

# Add richness and Shannon
input_filt$map_loaded$rich <- specnumber(input_filt$data_loaded, MARGIN = 2)
input_filt$map_loaded$shannon <- diversity(input_filt$data_loaded, 
                                                index = "shannon", MARGIN = 2)

# Species and family (input)
input$map_loaded$rich <- specnumber(tax_sum_species, MARGIN = 2)
range(input$map_loaded$rich)
input$map_loaded$richFam <- specnumber(tax_sum_family, MARGIN = 2)
range(input$map_loaded$richFam) # 1 to 207

# Plot seqs and rich against env
input$map_loaded$trnL_seqs <- colSums(input$data_loaded)
ggplot(input$map_loaded, aes(vegetation_type, trnL_seqs)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "trnL sequences") +
  theme_bw()
pdf("InitialFigs/trnL_speciesRichness.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(vegetation_type, rich)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "Plant species richness") +
  theme_bw()
dev.off()

input_filt$map_loaded$trnL_seqs <- colSums(input_filt$data_loaded)
ggplot(input_filt$map_loaded, aes(vegetation_type, trnL_seqs)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "trnL sequences") +
  theme_bw()
ggplot(input_filt$map_loaded, aes(vegetation_type, rich)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "trnL richness") +
  theme_bw()

ggplot(input_filt$map_loaded, aes(AI, trnL_seqs)) +
  geom_point() +
  labs(x = "Aridity index", y = "trnL sequences") +
  theme_bw()
ggplot(input_filt$map_loaded, aes(AI, rich)) +
  geom_point() +
  labs(x = "Aridity index", y = "trnL richness") +
  theme_bw()

ggplot(input_filt$map_loaded, aes(AI, rich)) +
  geom_point() +
  labs(x = "Aridity index", y = "trnL richness") +
  theme_bw()



#### __Composition ####
bc <- calc_dm(tax_sum_species)
nmds <- calc_ordination(bc, "nmds", input$map_loaded)
pdf("InitialFigs/trnL_nmds.pdf", width = 7, height = 5)
mctoolsr::plot_ordination(input, nmds, color_cat = "vegetation_type", hulls = T)
dev.off()

# Rarefy and try again
set.seed(530)
input_rare <- single_rarefy(input, 1000) # n = 55
sort(colSums(input_rare$data_loaded))
tax_sum_speciesR <- summarize_taxonomy(input = input_rare, 
                                       level = 7, 
                                       report_higher_tax = F)
bc <- calc_dm(tax_sum_speciesR)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
mctoolsr::plot_ordination(input_rare, nmds, color_cat = "vegetation_type", hulls = T)

# But we already knew species could be bad
# Need to do at family level, with and without rarefaction
bc <- calc_dm(tax_sum_family)
nmds <- calc_ordination(bc, "nmds", input$map_loaded)
pdf("InitialFigs/trnL_nmds_Fam.pdf", width = 7, height = 5)
mctoolsr::plot_ordination(input, nmds, color_cat = "vegetation_type", hulls = T)
dev.off()
set.seed(530)
input_rare <- single_rarefy(input, 1000) # n = 55
sort(colSums(input_rare$data_loaded))
tax_sum_familyR <- summarize_taxonomy(input = input_rare, 
                                       level = 5, 
                                       report_higher_tax = F)
bc <- calc_dm(tax_sum_familyR)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
mctoolsr::plot_ordination(input_rare, nmds, color_cat = "vegetation_type", hulls = T)

set.seed(530)
input_rare <- single_rarefy(input, 1500) # n = 27
sort(colSums(input_rare$data_loaded))
tax_sum_familyR <- summarize_taxonomy(input = input_rare, 
                                      level = 5, 
                                      report_higher_tax = F)
bc <- calc_dm(tax_sum_familyR)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
mctoolsr::plot_ordination(input_rare, nmds, color_cat = "vegetation_type", hulls = T)

set.seed(530)
input_rare <- single_rarefy(input, 2000) # n = 16
sort(colSums(input_rare$data_loaded))
tax_sum_familyR <- summarize_taxonomy(input = input_rare, 
                                      level = 5, 
                                      report_higher_tax = F)
bc <- calc_dm(tax_sum_familyR)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
mctoolsr::plot_ordination(input_rare, nmds, color_cat = "vegetation_type", hulls = T)

set.seed(530)
input_rare <- single_rarefy(input, 500) # n = 97
sort(colSums(input_rare$data_loaded))
tax_sum_familyR <- summarize_taxonomy(input = input_rare, 
                                      level = 5, 
                                      report_higher_tax = F)
bc <- calc_dm(tax_sum_familyR)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
mctoolsr::plot_ordination(input_rare, nmds, color_cat = "vegetation_type", hulls = T)

cliffplot_taxa_bars(input, 5, "vegetation_type")
cliffplot_taxa_bars(input, 4, "vegetation_type")
cliffplot_taxa_bars(input, 3, "vegetation_type") # Mostly Magnoliopsida



#### __Filt ####
# Filter taxa and samples do at family level
# Family needs to have at least 15 reads
# Sample needs to have at least 100 reads
input_filt <- input
tax_sum_family2 <- summarize_taxonomy(input = input_filt,
                                      relative = F,
                                      level = 5, 
                                      report_higher_tax = F)
input_filt$data_loaded <- tax_sum_family2
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 15) %>%
  mutate(ASV = rownames(.)) # 266
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 266 removed
min100 <- data.frame("count" = colSums(input_filt$data_loaded)) %>%
  filter(count < 100) %>%
  mutate(SampleID = rownames(.)) # 114
input_filt <- filter_data(input_filt,
                          "sampleID",
                          filter_vals = min100$SampleID) # 130 samples remaining

fam_rel <- as.data.frame(decostand(t(input_filt$data_loaded), method = "total", MARGIN = 1)) %>%
  t() %>%
  as.data.frame()
fam_prev <- data.frame("Family" = rownames(fam_rel),
                       "Present" = rowSums(fam_rel > 0))
input_filt$map_loaded$trnL_seqs <- colSums(input_filt$data_loaded)
range(input_filt$map_loaded$trnL_seqs) # 104 to 3676
input_filt$map_loaded$rich <- specnumber(input_filt$data_loaded, MARGIN = 2)
range(input_filt$map_loaded$rich) # 32 to 176
table(input_filt$map_loaded$vegetation_type)

pdf("InitialFigs/trnL_filt_seqs.pdf", width = 7, height = 5)
ggplot(input_filt$map_loaded, aes(vegetation_type, trnL_seqs, colour = vegetation_type)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 2, width = 0.2) +
  labs(x = "Vegetation", y = "# seqs", colour = NULL) +
  scale_colour_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                                 "darkgoldenrod", "chartreuse3", "brown")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("InitialFigs/trnL_filt_FamRichness.pdf", width = 7, height = 5)
ggplot(input_filt$map_loaded, aes(vegetation_type, rich, colour = vegetation_type)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 2, width = 0.2) +
  labs(x = "Vegetation", y = "Plant family richness", colour = NULL) +
  scale_colour_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                                "darkgoldenrod", "chartreuse3", "brown")) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

bc <- calc_dm(fam_rel)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
pdf("InitialFigs/trnL_filt_FamNMDS.pdf", width = 7, height = 5)
mctoolsr::plot_ordination(input_filt, nmds, color_cat = "vegetation_type", hulls = T) +
  scale_color_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                               "darkgoldenrod", "chartreuse3", "brown")) +
  scale_fill_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                               "darkgoldenrod", "chartreuse3", "brown"))
dev.off()

# What if high NA is driving this? Remove and rerun.
fam_rel_nona <- fam_rel[!(rownames(fam_rel) == "NA"), ]
bc <- calc_dm(fam_rel_nona)
nmds <- calc_ordination(bc, "nmds", input_rare$map_loaded)
pdf("InitialFigs/trnL_filt_FamNMDS_nona.pdf", width = 7, height = 5)
mctoolsr::plot_ordination(input_filt, nmds, color_cat = "vegetation_type", hulls = T) +
  scale_color_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                                "darkgoldenrod", "chartreuse3", "brown")) +
  scale_fill_manual(values = c("antiquewhite", "darkgreen", "gold", "purple",
                               "darkgoldenrod", "chartreuse3", "brown"))
dev.off()

pdf("InitialFigs/trnL_filt_FamTop11.pdf", width = 7, height = 5)
plot_taxa_bars(fam_rel, input_filt$map_loaded, "vegetation_type", 11, data_only = FALSE) +
  scale_fill_brewer(palette = "Paired") +
  labs(fill = "Family") +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal()
dev.off()



#### _BOLD ####
# If you go to ID a plant, it says the public plant database, from combined matK, rbcL, ITS is 453,647 non-redundant.
# https://bench.boldsystems.org/index.php/Public_SearchTerms
# If I search "Tracheophyta" in the data portal (vascular plants), this yields:
# Found 384,277 published records,
# with 384,277 records with sequences,
# forming 0 BINs (clusters),
# with specimens from 213 countries,
# deposited in 386 institutions.
# Of these records, 374,967 have species names, and represent 122,083 species.
# Try downloading all of those!

# Import the download
# Split by gene, then reformat for PhyloFlash
# See section 4.3 https://hrgv.github.io/phyloFlash/install.html
# Fasta headers should have the format {IDENTIFIER}.{INTEGER}.{INTEGER} {TAXONOMY-STRING} where:
#   IDENTIFIER is a unique sequence identifier which does not have spaces or periods
# The difference between the two INTEGERs should be the length of the sequence, e.g. 1.1700 for a 1700 bp sequence
# TAXONOMY-STRING is in SILVA or NCBI format, delimited by semicolons with no spaces (but spaces in taxon names allowed)
# There is a single space before the TAXONOMY-STRING
# Make full taxonomy string: for example the MFD is like this:
# FLASV2.1.1310 Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Xanthobacteraceae;Bradyrhizobium;MFD_s_2

bold <- microseq::readFasta("~/Desktop/Fierer/Strains/BOLD_plants.fasta")
bold2 <- bold %>%
  separate(Header, into = c("IDENTIFIER", "Species", "Gene"), sep = "\\|")
length(unique(bold2$Species)) # 123789
length(unique(bold2$Gene)) # 38??
table(bold2$Gene)

ITS2 <- bold2 %>%
  filter(Gene == "ITS2")
length(unique(ITS2$Species)) # 46061, 38810 after filtering
ITS2_spp <- ITS2 %>%
  group_by(Species) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  separate(Species, remove = F, into = c("Genus", "Spp"), sep = " ") %>%
  filter(is.na(Spp) == FALSE) %>%
  mutate(INTEGER1 = 1,
         INTEGER2 = nchar(Sequence)) %>%
  filter(INTEGER2 >= 150 & INTEGER2 <= 400) %>%
  mutate(TAXONOMY = gsub(" var. ", "_", Species)) %>%
  mutate(TAXONOMY = gsub(" ", ";", TAXONOMY)) %>%
  mutate(Domain = "Tracheophyta",
         Phylum = "Tracheophyta",
         Class = "Tracheophyta",
         Order = "Tracheophyta",
         Family = "Tracheophyta") %>%
  mutate(TAXONOMY = paste(Domain, Phylum, Class, Order, Family, TAXONOMY, sep = ";")) %>%
  mutate(Header = paste(IDENTIFIER, ".", INTEGER1, ".", INTEGER2, " ", TAXONOMY, 
                        sep = "")) %>%
  dplyr::select(Header, Sequence)
range(ITS2_spp$INTEGER2) # 108 to 967
hist(ITS2_spp$INTEGER2)

ITS <- bold2 %>%
  filter(Gene == "ITS") %>%
  separate(Species, remove = F, into = c("Genus", "S"), sep = " ")
length(unique(ITS$Species)) # 52434, 50703 after filtering
ITS_spp <- ITS %>%
  group_by(Species) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  separate(Species, remove = F, into = c("Genus", "Spp"), sep = " ") %>%
  filter(is.na(Spp) == FALSE) %>%
  mutate(INTEGER1 = 1,
         INTEGER2 = nchar(Sequence)) %>%
  filter(INTEGER2 >= 500 & INTEGER2 <= 900) %>%
  mutate(TAXONOMY = gsub(" var. ", "_", Species)) %>%
  mutate(TAXONOMY = gsub(" ", ";", TAXONOMY)) %>%
  mutate(Domain = "Tracheophyta",
         Phylum = "Tracheophyta",
         Class = "Tracheophyta",
         Order = "Tracheophyta",
         Family = "Tracheophyta") %>%
  mutate(TAXONOMY = paste(Domain, Phylum, Class, Order, Family, TAXONOMY, sep = ";")) %>%
  mutate(Header = paste(IDENTIFIER, ".", INTEGER1, ".", INTEGER2, " ", TAXONOMY, 
                        sep = "")) %>%
  dplyr::select(Header, Sequence)
range(ITS_spp$INTEGER2) # 134 to 3344
hist(ITS_spp$INTEGER2)
#microseq::writeFasta(ITS_spp, "~/Downloads/SILVA_CustomDB_ITS.fasta")

matK <- bold2 %>%
  filter(Gene == "matK")
length(unique(matK$Species)) # 54690
# 53628 removing those without species classification
# 53021 removing those >= 2000 bp
matK_spp <- matK %>%
  group_by(Species) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  separate(Species, remove = F, into = c("Genus", "Spp"), sep = " ") %>%
  filter(is.na(Spp) == FALSE) %>%
  mutate(INTEGER1 = 1,
         INTEGER2 = nchar(Sequence)) %>%
  filter(INTEGER2 < 2000) %>%
  mutate(TAXONOMY = gsub(" var. ", "_", Species)) %>%
  mutate(TAXONOMY = gsub(" ", ";", TAXONOMY)) %>%
  mutate(Domain = "Tracheophyta",
         Phylum = "Tracheophyta",
         Class = "Tracheophyta",
         Order = "Tracheophyta",
         Family = "Tracheophyta") %>%
  mutate(TAXONOMY = paste(Domain, Phylum, Class, Order, Family, TAXONOMY, sep = ";")) %>%
  mutate(Header = paste(IDENTIFIER, ".", INTEGER1, ".", INTEGER2, " ", TAXONOMY, 
                        sep = "")) %>%
  dplyr::select(Header, Sequence)
range(matK_spp$INTEGER2) # 131 to 5602
hist(matK_spp$INTEGER2)
#microseq::writeFasta(matK_spp, "~/Downloads/SILVA_CustomDB_matK.fasta")

rbcL <- bold2 %>%
  filter(Gene == "rbcL")
length(unique(rbcL$Species)) # 38132
# 53628 removing those without species classification
# 53021 removing those >= 2000 bp
rbcL_spp <- rbcL %>%
  group_by(Species) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  separate(Species, remove = F, into = c("Genus", "Spp"), sep = " ") %>%
  filter(is.na(Spp) == FALSE) %>%
  mutate(INTEGER1 = 1,
         INTEGER2 = nchar(Sequence)) %>%
  filter(INTEGER2 < 2000) %>%
  mutate(TAXONOMY = gsub(" var. ", "_", Species)) %>%
  mutate(TAXONOMY = gsub(" ", ";", TAXONOMY)) %>%
  mutate(Domain = "Tracheophyta",
         Phylum = "Tracheophyta",
         Class = "Tracheophyta",
         Order = "Tracheophyta",
         Family = "Tracheophyta") %>%
  mutate(TAXONOMY = paste(Domain, Phylum, Class, Order, Family, TAXONOMY, sep = ";")) %>%
  mutate(Header = paste(IDENTIFIER, ".", INTEGER1, ".", INTEGER2, " ", TAXONOMY, 
                        sep = "")) %>%
  dplyr::select(Header, Sequence)
range(rbcL_spp$INTEGER2) # 131 to 5602
hist(rbcL_spp$INTEGER2)
#microseq::writeFasta(rbcL_spp, "~/Downloads/SILVA_CustomDB_rbcL.fasta")



#### __ITS ####
# Need to see if the species in the veg data have ITS seqs
sum(d_species_uniq$Spp %in% ITS$Species) # 48 (this is 35 more than chloroplast)
sum(unique(d_species_uniq$Genus) %in% unique(ITS$Genus)) # 67 genera
d_species_db <- d_species_uniq %>%
  filter(Spp %in% ITS$Species)
d_species_uniq$match <- sapply(d_species_uniq$Spp, function(x) any(grepl(x, ITS$Species)))
# 52 (using grep found 4 more)


# Analyze phyloFlash results
setwd("~/Desktop/Fierer/Strains/Australia/pf_ITS/")
files <- list.files()
# There are only 156 files. The other 175 had errors.
# Error was: Can't use string as a HASH ref while "strict refs" in use
# For now, proceed, see how these data look
merged_data <- files %>%
  lapply(function(file) {
    sample_id <- gsub("\\.phyloFlash.NTUfull_abundance\\.csv", "", basename(file))
    data <- read.csv(file, header = FALSE) %>%
      set_names(c("NTU", sample_id)) %>%
      return(data)
  }) %>%
  reduce(function(x, y) merge(x, y, by = "NTU", all = TRUE))
merged_data[is.na(merged_data)] <- 0
sample_id <- data.frame(fullname = list.files()) %>%
  separate(fullname, into = c("sampleID", "Junk1", "Junk2", sep = "_"))
names(merged_data) <- c("NTU", sample_id$sampleID)
setwd("~/Documents/GitHub/AussieStrains/")

# Now see if those 48 species with ITS seqs were detected
d_species_wPF <- d_species %>%
  filter(SampleID %in% names(merged_data))
d_species_wPF_uniq <- d_species_uniq %>%
  filter(SampleID %in% names(merged_data))
length(unique(d_species$SampleID)) # 48
length(unique(d_species_wPF$SampleID)) # 39 overlap
length(unique(d_species_wPF$Spp)) # 111 of the 122 in the overlapping dataset
veg_pf <- merged_data %>%
  dplyr::select(NTU, all_of(unique(d_species_wPF$SampleID))) %>%
  mutate(Tot = rowSums(.[2:ncol(.)])) %>%
  filter(Tot > 0) %>%
  arrange(desc(Tot)) %>%
  separate(NTU, into = c("Domain", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"), sep = ";", remove = F) %>%
  mutate(Spp = paste(Genus, Species, sep = " "))
sum(d_species_wPF_uniq$Spp %in% ITS$Species) # 38
sum(unique(d_species_wPF_uniq$Genus) %in% unique(ITS$Genus)) # 54
sum(d_species_uniq$Spp %in% veg_pf$Spp) # 3 of 38 possible identified
sum(unique(d_species_uniq$Genus) %in% unique(veg_pf$Genus)) # 18 of 54 possible

# Change parentheses to NA
length(unique(merged_data$NTU)) # 23391
seqtab_wTax <- merged_data %>%
  mutate(RowID = row_number()) %>%
  mutate(OTU = paste("OTU", RowID, sep = "_")) %>%
  separate(NTU, into = c("Domain", "Phylum", "Class", "Order", "Family", 
                         "Genus", "Species"), sep = ";") %>%
  mutate(across(everything(), ~ ifelse(grepl("\\(|\\)", .), "NA", .))) %>%
  mutate(taxonomy = paste(Domain, Phylum, Class, Order, Family, Genus, Species,
                          sep = ";")) %>%
  dplyr::select(-RowID) %>%
  dplyr::select(OTU, everything(), taxonomy)
names(seqtab_wTax)
#out_fp <- "data/seqtab_wTax_mctoolsr_ITS.txt"
#write("#Exported for mctoolsr", out_fp)
#suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Import
tax_table_fp <- "data/seqtab_wTax_mctoolsr_ITS.txt"
map_fp <- "~/Documents/GitHub/AussieStrains/data/metadata_331.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 156 samples loaded
input$map_loaded$sampleID <- rownames(input$map_loaded)
input <- filter_data(input, 
                     filter_cat = "sampleID", 
                     keep_vals = d_sylph$sampleID) # 128 of 268 have plant data
# Check singletons and doubletons
singdoub <- data.frame("count" = rowSums(input$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 15154 Don't remove, since will aggregate
input_filt <- filter_taxa_from_input(input,
                                     taxa_IDs_to_remove = singdoub$ASV) # 15154 removed

# Rarefaction
rarecurve(t(input$data_loaded), step = 100, label = F)
rarecurve(t(input_filt$data_loaded), step = 100, label = F)

# That is very unrealistic. Need to go ahead and aggregate to species, genus, family etc.
tax_sum_genus <- summarize_taxonomy(input = input, 
                                    level = 6, 
                                    report_higher_tax = F) # 5090 genera
tax_sum_species <- summarize_taxonomy(input = input, 
                                      level = 7, 
                                      report_higher_tax = F) # 6421 species
sort(colSums(input_filt$data_loaded))
plot(sort(colSums(input_filt$data_loaded)))
mean(colSums(input_filt$data_loaded)) # 440.0625
se(colSums(input_filt$data_loaded)) # 72.91068

# Add richness and Shannon
input_filt$map_loaded$rich <- specnumber(input_filt$data_loaded, MARGIN = 2)
input_filt$map_loaded$shannon <- diversity(input_filt$data_loaded, 
                                           index = "shannon", MARGIN = 2)

# Plot seqs and rich against env
ITS_seqs <- input_filt$map_loaded
input_filt$map_loaded$ITS_seqs <- colSums(input_filt$data_loaded)
ggplot(input_filt$map_loaded, aes(vegetation_type, ITS_seqs)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "ITS sequences") +
  theme_bw()
ggplot(input_filt$map_loaded, aes(vegetation_type, rich)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.2) +
  labs(x = "Vegetation", y = "ITS richness") +
  theme_bw()

ggplot(input_filt$map_loaded, aes(AI, ITS_seqs)) +
  geom_point() +
  labs(x = "Aridity index", y = "ITS sequences") +
  theme_bw()
ggplot(input_filt$map_loaded, aes(AI, rich)) +
  geom_point() +
  labs(x = "Aridity index", y = "ITS richness") +
  theme_bw()

ggplot(input_filt$map_loaded, aes(ITS_seqs, rich)) +
  geom_point() +
  labs(x = "ITS sequences", y = "ITS richness") +
  theme_bw()
