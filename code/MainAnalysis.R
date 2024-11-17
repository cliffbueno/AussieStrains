# Analysis for soil bacterial strain manuscript
# By Cliff Bueno de Mesquita, Fierer Lab, Fall 2024

#### 1. Setup ####
# Libraries
library(tidyverse)
library(ggside)
library(ozmaps)
library(sf)
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

# Functions
`%notin%` <- Negate(`%in%`)
source("~/Documents/GitHub/SunflowerGxE/code/cliffplot_taxa_bars.R")

# Working directory
setwd("~/Documents/GitHub/AussieStrains/")

d <- read.csv("data/metadata104.csv") %>%
  dplyr::select(-X, -Prefix)

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
              #"water_content", # 51 NA
              "bio1", "bio12", "AI")
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
              "Sand", "Silt", "S", "Temp", "Precip", "Aridity")) %>%
  dplyr::select(Clay, Silt, Sand, `Org. C`, NO3, `Colwell P`, pH, Conductivity,
                B, Cu, Fe, Mn, Zn, Al, Ca, Mg, K, Na, S,
                Lat, Long, Temp, Precip, Aridity)

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



#### 4. Ref Genomes ####
# Subset the GTDB full database to the genera of interest
bacGT <- read.delim("~/Desktop/Strains/bac120_metadata_r220.tsv") # takes a while to load
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
# Don't forget dos2unix!!
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
meta_key <- meta %>%
  dplyr::select(sampleID, vegetation_type, AI) %>%
  mutate(SampleID = as.character(sampleID)) %>%
  dplyr::select(-sampleID)
meta_AI <- meta_key %>%
  dplyr::select(-vegetation_type)



#### __95% ANI ####
sylph_profile <- read.delim("data/sylph_profile_ani95.tsv")
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
genus_info <- read_xlsx("data/genus_info_75point1.xlsx") %>%
  left_join(., sylph_profile_genera, by = "Genus") %>%
  left_join(., sylph_profile_genera_samples, by = "Genus") %>%
  replace_na(list(n_Genomes_Sylph95 = 0,
                  n_Samples_Sylph95 = 0))
ggplot(genus_info, aes(n_Genomes_downloaded, n_Genomes_Sylph95)) +
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
  filter(Genus == "Bradyrhizobium")
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



# Relative abundance
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



#### 5. Ref Genomics ####
# Analyze the different genomes from the Sylph output
# Use coverM to calculate coverage and relative abundance of the genomes with Sylph hits


#### 6. Strain Genomics ####
# Run StrainFinder on KBase. Use the coverage information to inform the parameters
# Needed if different reference genomes not detected or if too few samples

# Assess gene content

# Import DRAM annotation for the strain genomes
# Compute shared and unique KOs

#### 7. Pop. Genomics ####
# Assess strain population genomics from .vcf files
