# Analysis for Bradyrhizobium strain-level diversity and ecology manuscript
# For submission to ISME J
# By Cliff Bueno de Mesquita, Fierer Lab, Fall 2024 to Spring 2025

# Document outline (see document outline on right)
# 1. Setup
# 2. Environment
# 3. Dominant taxa (mTAGs)
# 4. Reference genomes
# 5. Strain Info (ANI, 16S)
# 6. Strain Distributions/Ecology
# 7. Strain Pangenomics



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
save_pheatmap_pdf <- function(x, filename, width = 6, height = 8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Working directory
setwd("~/Documents/GitHub/AussieStrains/")

# Metadata (first pass was 104 samples, later will use all 331 samples)
d <- read.csv("data/metadata104.csv") %>%
  dplyr::select(-X, -Prefix) %>%
  mutate(sampleID = as.character(sampleID))

# Output just sampleID and coordinates
c <- d %>%
  dplyr::select(sampleID, latitude, longitude)
#write.csv(c, "data/coords104.csv")

meta <- d %>%
  mutate(Sample = paste("X", sampleID, sep = "")) %>%
  dplyr::select(Sample, everything())

# Save metadata as text for mctoolsr import
length(unique(meta$Sample))
#write.table(meta, "data/metadata104.txt", row.names = F, sep = "\t")



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
#pdf("InitialFigs/Env_Corrplot.pdf", width = 7, height = 5)
corrplot(m, 
         method = "square",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5)
#dev.off()

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
#pdf("InitialFigs/Env_PCA.pdf", width = 7, height = 5)
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
#dev.off()



#### 3. Taxa ####
# Use mTAGs to select abundant and ubiquitous soil bacteria
# mTAGs identifies and classifies SSU sequences from the metaG
# Note: Originally considered other taxa but then decided to focus on Bradyrhizobium

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
map_fp <- "data/metadata104.txt"
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
#pdf("InitialFigs/Alpha_Rich_AI.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, alpha = 0.75, pch = 16, aes(colour = vegetation_type)) +
  geom_smooth(method = "lm") +
  labs(x = "Aridity index",
       y = "Richness (# mTAG OTUs)",
       colour = "Vegetation") +
  scale_colour_viridis_d() +
  theme_bw()
#dev.off()
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
              "Sand", "Silt", "S", "H2O", "Temp", "Precip", "Aridity")) %>%
  dplyr::select(Clay, Silt, Sand, `Org. C`, NO3, `Colwell P`, pH, Conductivity,
                B, Cu, Fe, Mn, Zn, Al, Ca, Mg, K, Na, S,
                Lat, Long, H2O, Temp, Precip, Aridity)
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
#pdf("InitialFigs/Beta_Bray_PCoA.pdf", width = 7, height = 5)
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
#dev.off()

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
#View(input_filt_rare$taxonomy_loaded)

# Get focal taxa and samples they are most abundant in to try StrainFinder
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
length(unique(genus_otus$Genus))
# 2781 OTUs across those 40 genera! Most Bryobacter. Bradyrhizobium with 57.

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
#pdf("InitialFigs/TopGenera_AI_75point1_used33.pdf", width = 12, height = 8)
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
#dev.off()

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
#pdf("InitialFigs/otus_brady.pdf", width = 8, height = 8)
g
#dev.off()

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
#pdf("InitialFigs/otus_strepto.pdf", width = 8, height = 8)
g
#dev.off()

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
#pdf("InitialFigs/otus_udaeo.pdf", width = 8, height = 8)
g
#dev.off()

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
#pdf("InitialFigs/otus_myco.pdf", width = 8, height = 8)
g
#dev.off()



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
# 35/40. Missing Acidibacter, Tundrisphaera, Burkholderia-Caballeronia-Paraburkholderia, Actinomadura, P3OB-42. Don't want to focus on those anyways since they are messy.

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
# This got a little messy since there were download issues so had to do multiple batches

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

# Unzipped, then extracted fastas, got up to 21467 of 21471 .fna files
# Need to get the last 4! Do manually.
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
         height = 8
         #,
         #filename = "InitialFigs/Sylph_Brady.png"
         )
#dev.off()
#dev.set(dev.next())

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
         height = 8
         #,
         #filename = "InitialFigs/Sylph_Strepto.png"
         )
#dev.off()
#dev.set(dev.next())

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
         height = 8
         #,
         #filename = "InitialFigs/Sylph_Udaeo.png"
         )
#dev.off()
#dev.set(dev.next())

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
         height = 8
         #,
         #filename = "InitialFigs/Sylph_Myco.png"
         )
#dev.off()
#dev.set(dev.next())



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
         height = 8
         #,
         #filename = "InitialFigs/compare_sylph_coverm.png",
         )
#dev.off()
#dev.set(dev.next())

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
# Use this! (Dylan Chivian from KBase confirmed it's a good idea)
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
#write_xlsx(myco_info, "data/Mycobacteria_StrainFinder_Info.xlsx", format_headers = F)

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



#### 5. Strain Info ####
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
#pdf("InitialFigs/Strains_RelAbund.pdf", width = 7, height = 5)
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
#dev.off()
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
#pdf("InitialFigs/Strains_MissingBac120.pdf", width = 7, height = 5)
ggplot(bac120, aes(strainID, number_missing_genes)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, pch = 16, width = 0.25) +
  labs(x = "Genome",
       y = "Number of missing bac120 genes") +
  ggtitle("Bradyrhizobium genomes, n = 1081") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
#dev.off()
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
brady181_tree <- read.tree("data/brady181_fasttree_nogap.tree")
brady181_tree$tip.label
brady181_tree$tip.label <- gsub(".fna", "", brady181_tree$tip.label)
brady181_tree$tip.label <- gsub("Brady_", "", brady181_tree$tip.label)
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

# Now just EMP 515-806 part of the 16S
combined_16S <- microseq::readFasta("data/combined_16S.fasta")
EMP <- combined_16S %>%
  mutate(Sequence = substr(Sequence, start = 515, stop = 806))
#microseq::writeFasta(EMP, "data/EMP_16S.fasta")
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

# Figure S8
# Need to plot % identity on X and # comparisons on y
facet_names <- c("16S EMP" = "515-806 16S rRNA gene (n OTUs)",
                 "16S" = "Full length 16S rRNA gene (n OTUs)",
                 "ANI" = "ANI (n genomes)")
figS8 <- ggplot(otu_comb, aes(Cutoff, Count)) +
  geom_bar(stat = "identity") +
  labs(x = "Similarity cutoff (%)",
       y = NULL) +
  facet_wrap(~ Metric, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 8))
figS8
pdf("FinalFigs/FigureS8.pdf", width = 7, height = 5)
figS8
dev.off()
png("FinalFigs/FigureS8.png", width = 7, height = 5, units = "in", res = 300)
figS8
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



#### 6. Strain Distributions ####
# Reran Sylph but removed Strain 3 and 4 and GTDB with > 4 missing bac120
# So, now 915 input genomes
# Run on all 331 BASE metagenomes that met criteria!
# Make tree with those detected by Sylph
# Rerun ecological analyses with new input data, new UniFrac
# Import the needed data so you can restart the analysis just from here
d_331 <- read.csv("data/metadata_331.csv") %>%
  mutate("Climate Class" = ifelse(AI < 0.5, "Arid to semi-arid",
                                  ifelse(AI >= 0.5 & AI < 0.65, "Dry sub-humid",
                                         "Humid")))
range(d_331$bio1)
range(d_331$bio12)
bacGT <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv") # takes a while to load
input <- readRDS("data/input_sylph.rds")
row_dat <- readRDS("data/row_data.rds")
ubiq <- data.frame(GenomeID = rownames(input$data_loaded),
                   Ubiquity = rowSums(input$data_loaded > 0)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(Perc = Ubiquity/331*100) %>%
  mutate(Nfix = ifelse(Function == "N fix. Free" | Function == "N fix. Sym.",
                       "Yes", "No"))
tree_data_maxmin <- readRDS("data/tree_data_maxmin.rds")
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
d_brady <- d_331 %>%
  filter(sampleID %in% ref_ani$sampleID) %>%
  mutate(sampleID = as.character(sampleID)) %>%
  arrange(sampleID)

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
# com_samples <- sylph_strains_331 %>% # Need to make sylph_strains_331 first
#   filter(rownames(.) %in% com_ubiq$GenomeID) %>%
#   select_if(function(x){!all(is.na(x))}) %>%
#   rownames_to_column(var = "GenomeID")
# names(com_samples)
#write.csv(com_samples, "data/commercial_detection.csv", row.names = F)

# Check where commercial samples were found
# d_com <- d_331 %>%
#   filter(sampleID %in% names(com_samples))
# coords_trans <- st_as_sf(d_com, 
#                          coords = c('longitude', 'latitude'), 
#                          crs=4326)
# sf_oz <- ozmap("states")
# map <- ggplot(data = sf_oz) + 
#   geom_sf(fill = "grey90", color = "white") +
#   geom_sf(data = coords_trans,
#           aes(fill = AI, shape = `Climate Class`), 
#           size = 3, alpha = 1, color = "black", stroke = 0.3) +
#   scale_shape_manual(values = c(21, 24, 22)) +
#   annotation_scale() +
#   annotation_north_arrow(pad_x = unit(1, "cm"), pad_y = unit(1, "cm"),
#                          height = unit(1, "cm"), width = unit(1, "cm"),) +
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) +
#   guides(shape = guide_legend(order = 2),
#          fill = guide_colorbar(order = 1)) +
#   xlim(111, 155) +
#   ylim(43, 12) +
#   labs(fill = "Aridity\nindex") +
#   theme_minimal()
# pdf("InitialFigs/CommercialMap.pdf", width = 7, height = 5)
# map
# dev.off()



#### _Setup ####
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
# row_dat <- data.frame(GenomeID = rownames(sylph_strains_331)) %>%
#   left_join(., gtdb89_use, by = c("GenomeID" = "ncbi_genbank_assembly_accession")) %>%
#   replace_na(list(genome_size = 8085095, Type = "StrainFinder", Source = "Soil")) %>%
#   mutate(Type = ifelse(GenomeID == "Reference", "StrainFinder Ref", Type)) %>%
#   mutate(Type = ifelse(GenomeID %in% com_ubiq$GenomeID, "Commercial", Type)) %>% # Add later
#   left_join(., checkm_181, by = "GenomeID") %>%
#   mutate(GenomeSize = 100 * genome_size / checkm_completeness) %>%
#   left_join(., brady_func, by = "GenomeID") %>%
#   mutate(Function = ifelse(Nfix >= 5 & Nod >= 4, "N fix. Sym.",
#                            ifelse(Nfix >= 5 & Nod < 4, "N fix. Free",
#                                   ifelse(Photo > 0, "Photosyn.", " "))))
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
hm <- pheatmap(sylph_strains_331,
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
               show_colnames = F)
save_pheatmap_pdf(hm, "FinalFigs/FigureS4forPPT.pdf")

# Remake with presence/absence and show black, will make easier to see
# Use this as the main figure and have rel abund one as supp
sylph_strains_331_pa <- sylph_strains_331
sylph_strains_331_pa[sylph_strains_331_pa > 0] <- 1
hm <- pheatmap(sylph_strains_331_pa,
               breaks = seq(0, 1, length.out = 3),
               color = c("black"),
               legend = F,
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
               show_colnames = F)
save_pheatmap_pdf(hm, "FinalFigs/Figure3forPPT.pdf")
  


### _Tree ####
#### __915 ####
# Show big tree as Figure S3
# Made this for 993 genomes before more stringent filtering to 915
# Can prune to 915
brady993_tree <- read.tree("data/brady993_fasttree_nogap.tree")
brady993_tree$tip.label
brady993_tree$tip.label <- gsub(".fna", "", brady993_tree$tip.label)
toPrune <- data.frame(GenomeID = brady993_tree$tip.label) %>%
  filter(GenomeID %notin% c(brady915_ncbi$Header, "Brady_Reference"))
brady915_tree <- drop.tip(brady993_tree,
                          tip = toPrune$GenomeID)
brady915_tree$tip.label <- gsub("Brady_", "", brady915_tree$tip.label)
ggtree(brady915_tree)
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
pdf("FinalFigs/FigureS3.pdf", width = 8.5, height = 17)
ggtree(brady915_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 0.5, vjust = 0.5, aes(color = Sylph, label = Species)) +
  scale_color_manual(values = c("red", "grey50")) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()
png("FinalFigs/FigureS3.png", width = 8.5, height = 17, units = "in", res = 300)
ggtree(brady915_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 0.5, vjust = 0.5, aes(color = Sylph, label = Species)) +
  scale_color_manual(values = c("red", "grey50")) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
dev.off()



#### __181 ####
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
#pdf("InitialFigs/brady181_Fasttree_Bac120_2026AA.pdf", width = 8.5, height = 11)
ggtree(brady181_tree, linewidth = 0.1)  %<+% tc +
  geom_tiplab(size = 1, vjust = 0.5, aes(label = Species)) +
  geom_treescale() +
  geom_nodelab(aes(label = label), size = 0.75) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = c(0.7, 0.3))
#dev.off()

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
  left_join(., row_dat, by = "GenomeID") %>%
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
  geom_treescale() +
  theme(plot.margin = margin(0,-20,0,-15))
p

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

pdf("FinalFigs/FigureS6.pdf", width = 8, height = 8)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()
png("FinalFigs/FigureS6.png", width = 8, height = 8, units = "in", res = 300)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()

# Try to make 3 groupings
# Genomes detected in at least 2 samples
# All samples only in a specific climate class
ubiq <- data.frame(GenomeID = rownames(input$data_loaded),
                   Ubiquity = rowSums(input$data_loaded > 0)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(Perc = Ubiquity/331*100) %>%
  mutate(Nfix = ifelse(Function == "N fix. Free" | Function == "N fix. Sym.",
                       "Yes", "No"))
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
        legend.position.inside = c(1, 0),
        legend.justification = c(1, 0),
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
  summarise(minT = min(bio1),
            maxT = max(bio1)) %>%
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
  geom_treescale() +
  theme(plot.margin = margin(0,-20,0,-15))
p

s <- ggplot(tree_data, aes(x = bio1, y = GenomeID, colour = Type)) +
  geom_segment(data = tree_data_maxmin, 
               aes(x = minT, xend = maxT, y = GenomeID, yend = GenomeID,
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

pdf("FinalFigs/FigureS7.pdf", width = 8, height = 8)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()
png("FinalFigs/FigureS7.png", width = 8, height = 8, units = "in", res = 300)
plot_grid(p, s, align = "h", rel_widths = c(0.5, 0.5))
dev.off()

# Try to make 2 groupings
# Genomes detected in at least 2 samples
# Narrow range vs. wide range
ubiq <- data.frame(GenomeID = rownames(input$data_loaded),
                   Ubiquity = rowSums(input$data_loaded > 0)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(Perc = Ubiquity/331*100) %>%
  mutate(Nfix = ifelse(Function == "N fix. Free" | Function == "N fix. Sym.",
                       "Yes", "No"))
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
ggplot(ubiq, aes(GenomeSize.x, range)) +
  geom_point(aes(fill = Type.x), pch = 21, size = 3) +
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
ggplot(ubiq_two, aes(GenomeSize.x, range)) +
  geom_point(aes(fill = Type.x), pch = 21, size = 3) +
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
  geom_point(aes(fill = Type.x), pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Ubiquity",
       y = "Temperature range",
       fill = "Type") +
  scale_fill_manual(values = c("#EE3377",
                               "#66CCEE",
                               "#332288",
                               "#EE7733",
                               "yellow")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
summary(lm(range ~ Ubiquity, data = ubiq_two)) # R2 = 0.30, p < 0.001
summary(lm(range ~ GenomeSize.x, data = ubiq_two)) # R2 = 0.04, p = 0.08

# Do a third time but for pH? (range 3.1 to 8.6, 3.2 to 8.6 for Sylph)



#### _Comp ####
# Composition and drivers, Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
# sp_comp <- sylph_strains_331 %>%
#   replace(is.na(.), 0)
# input <- list()
# input$map_loaded <- d_sylph_ai_sort
# rownames(input$map_loaded) <- d_sylph_ai_sort$sampleID
# sum(d_sylph_ai_sort$sampleID != names(sp_comp))
# input$data_loaded <- sp_comp
# input$taxonomy_loaded <- sp_comp %>%
#   mutate(taxonomy1 = "Bacteria",
#          taxonomy2 = "Proteobacteria",
#          taxonomy3 = "Alphaproteobacteria",
#          taxonomy4 = "Hyphomicrobiales",
#          taxonomy5 = "Nitrobacteraceae",
#          taxonomy6 = "Bradyrhizobium",
#          taxonomy7 = rownames(.)) %>%
#   dplyr::select(taxonomy1, taxonomy2, taxonomy3, taxonomy4,
#                 taxonomy5, taxonomy6, taxonomy7)
# sum(rownames(input$data_loaded) != rownames(input$taxonomy_loaded))
#saveRDS(input, "data/input_sylph.rds")
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
#pdf("InitialFigs/BradyComp_RichAI.pdf", width = 7, height = 5)
ggplot(input$map_loaded, aes(AI, rich)) +
  geom_point(size = 3, pch = 16, alpha = 0.5) +
  geom_smooth(se = F) +
  labs(x = "drier <= Aridity index => wetter",
       y = "Number of genomes detected") +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw() +  
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
#dev.off()

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
mantel(bc, Wun, permutations = 2000) # r = 0.23, p = 0.0005

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
# Reload those data here
row_dat <- readRDS("data/row_data.rds")
ubiq <- data.frame(GenomeID = rownames(input$data_loaded),
                   Ubiquity = rowSums(input$data_loaded > 0)) %>%
  left_join(., row_dat, by = "GenomeID") %>%
  mutate(Perc = Ubiquity/331*100) %>%
  mutate(Nfix = ifelse(Function == "N fix. Free" | Function == "N fix. Sym.",
                       "Yes", "No"))
table(ubiq$Function)
table(ubiq$Nfix) # 145 yes to 36 no
t.test(GenomeSize ~ Nfix, data = ubiq) # No difference, p = 0.52
ggplot(ubiq, aes(Nfix, GenomeSize)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()

t.test(Ubiquity ~ Nfix, data = ubiq) # p = 0.01
ggplot(ubiq, aes(Nfix, Ubiquity)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25) +
  theme_bw()
# Careful - this could be driven by the StrainFinder MAGs! Rerun below

summary(lm(GenomeSize ~ Ubiquity, data = ubiq))
cor.test(ubiq$GenomeSize, ubiq$Ubiquity, method = "pearson")
figS5 <- ggplot(ubiq, aes(Ubiquity, GenomeSize)) +
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
figS5
pdf("FinalFigs/FigureS5.pdf", width = 8, height = 6)
figS5
dev.off()
png("FinalFigs/FigureS5.png", width = 8, height = 6, units = "in", res = 300)
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
#pdf("InitialFigs/UbiquityRankAbund.pdf", width = 8, height = 6)
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
#dev.off()



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
ggplot(input$map_loaded, aes(Axis01, Axis02)) +
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
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))

pcoa <- cmdscale(un, k = nrow(input$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, d_env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5)

#### _NMDS ####
# Just do on Weighted UniFrac for the manuscript. Figure S9.
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
figS9 <- ggplot(input$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.4) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = NMDS1, y = NMDS2, label = shortnames),
                  size = 3, color = "red") +
  geom_text(data = NULL, aes(x = -0.5, y = -0.3, label = "stress = 0.13"), size = 3,
            check_overlap = T) +
  labs(x = "NMDS1", 
       y = "NMDS2") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(size = 12), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
figS9
pdf("FinalFigs/FigureS9.pdf", width = 7, height = 5)
figS9
dev.off()
png("FinalFigs/FigureS9.png", width = 7, height = 5, units = "in", res = 300)
fig4
dev.off()

m <- adonis2(Wun ~ input$map_loaded$`Climate Class`)
m # Sig but very low R2 (0.03), and not much clustering. Could be driven by dispersion.
m <- betadisper(Wun, input$map_loaded$`Climate Class`)
anova(m) # Sig
m <- adonis2(Wun ~ input$map_loaded$vegetation_type)
m # Sig R2 = 0.10
m <- betadisper(Wun, input$map_loaded$vegetation_type)
anova(m) # Sig



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
# Here let's partition variation into geography vs. environment
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

fig5a <- function() {
  par(
    mar = c(1, 2.8, 1, 0.45),
    mfrow = c(1,1)
  )
  plot(mod, bg = c("#F8766D", "#619CFF"), Xnames = c('Geog.', 'Env.'))
}
ggdraw(fig5a)

# Need to run the GDM section first to make fig5b, then can come back here and combine
pdf("FinalFigs/Figure4.pdf", width = 7, height = 5)
plot_grid(fig5a, fig5b, ncol = 1, labels = "auto", vjust = 1, hjust = -1)
dev.off()
png("FinalFigs/Figure4.png", width = 7, height = 5, units = "in", res = 300)
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
mantel.partial(Wun_sub, dist.geog, dist.env, permutations = 2000) # p = 0.75
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
# Note: Temp available for all samples, so can use geographic distance for full dataset
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
# R2 lower than temperature



#### __GDM ####
# Again, test geog.dist and env.dist
# This will complement the variation partitioning

# Remake as full, not just lower triangle
dist.env <- as.matrix(dist(d_env_sig, method = "euclidean", diag = T, upper = T))
rownames(dist.env) <- input_sub$map_loaded$sampleID
colnames(dist.env) <- input_sub$map_loaded$sampleID

# Remake with n = 216
dist.geog <- geosphere::distm(cbind(input_sub$map_loaded$longitude, 
                                    input_sub$map_loaded$latitude),
                              fun = distHaversine)
rownames(dist.geog) <- input_sub$map_loaded$sampleID
colnames(dist.geog) <- input_sub$map_loaded$sampleID


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



#### 7. Strain Pangenomics ####
# Pangenomics of the 181 detected strains

#### _ANI ####
# Want to quickly report the range of ANI among these 181
ani_181 <- read.delim("data/ANI_ani_181.txt") %>%
  column_to_rownames(var = "layers")
range(ani_181[ani_181 < 1]) # 0.795487 to 0.999997



#### _KO new ####
# In order to get KO IDs, had to do individually with anvi-export-functions
# Note: This yields more KOs than the bulk way!
setwd("data/KO_tables")
# Import them all as a list
ko_list <- list()
ko_files <- list.files()
strainIDs <- substr(ko_files, start = 1, stop = nchar(ko_files) - 4)
for (i in 1:181) {
  ko_list[[i]] <- read.delim(ko_files[i], quote = "") %>%
    group_by(accession) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(presence = 1) %>%
    dplyr::select(accession, presence) %>%
    set_names(c("KO", strainIDs[i]))
}
# Merge them all
kos_merged <- ko_list[[1]]
for (i in 2:181) {
  kos_merged <- kos_merged %>%
    full_join(., ko_list[[i]], by = "KO")
}
kos_merged[is.na(kos_merged)] <- 0
setwd("~/Documents/GitHub/AussieStrains/")

brady_mat <- kos_merged %>%
  column_to_rownames(var = "KO") %>%
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
  select(-sum) # 1168
k180 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 180) %>%
  select(-sum) # 240
k1 <- brady_mat %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 1) %>%
  select(-sum) # 311

## For heatmap, use the Tao et al. 2021 nod and nif genes. nodABCIJ and nifBDEHKN.
## And Chris Greening energy acquisition strategies (Bay et al. 2021 genes; see Greening section)
# K14658 nodA; nodulation protein A [EC:2.3.1.-]
# K14659 nodB; chitooligosaccharide deacetylase [EC:3.5.1.-]
# K14666 nodC; N-acetylglucosaminyltransferase [EC:2.4.1.-]
# K09695 nodI; lipooligosaccharide transport system ATP-binding protein
# K09694 nodJ; lipooligosaccharide transport system permease protein

# K02585 nifB; nitrogen fixation protein NifB
# K02586 nifD; nitrogenase molybdenum-iron protein alpha chain [EC:1.18.6.1]
# K02587 nifE; nitrogenase molybdenum-cofactor synthesis protein NifE
# K02588 nifH; nitrogenase iron protein NifH
# K02591 nifK; nitrogenase molybdenum-iron protein beta chain [EC:1.18.6.1]
# K02592 nifN; nitrogenase molybdenum-iron protein NifN

# Plus photosynthesis genes
# K13991 puhA; photosynthetic reaction center H subunit
# K08929 pufM; photosynthetic reaction center M subunit
# K08928 pufL; photosynthetic reaction center L subunit

# Plus trace gas metabolism genes (CO, CH4, H2 (nope))
# K03520 coxL, cutL; aerobic carbon-monoxide dehydrogenase large subunit [EC:1.2.5.3]
# K10944 pmoA-amoA; methane/ammonia monooxygenase subunit A [EC:1.14.18.3 1.14.99.39]
greening <- read_xlsx("data/GreeningGenes.xlsx") %>%
  filter(KO %in% kos_merged$KO)
missing <- read_xlsx("data/GreeningGenes.xlsx") %>%
  filter(KO %notin% kos_merged$KO)
# # K00192 cdhA; anaerobic carbon-monoxide dehydrogenase, CODH/ACS complex subunit alpha [EC:1.2.7.4] is missing
# All of the hydrogenases are missing
# Formate dehydrogenase missing
# K00093 mdh; methanol dehydrogenase [EC:1.1.1.244] is missing

# Methylotrophy
# K23995 xoxF; lanthanide-dependent methanol dehydrogenase [EC:1.1.2.10] (Bao et al. 2014)

# Denitrification suite (own pathway, but also trace gas)
# And remember, can be complete or partial (see Shan et al. 2021)
# K00370 narG, narZ, nxrA; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]
# K02567 napA; nitrate reductase (cytochrome) [EC:1.9.6.1]
# K00368 nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]
# K15864 nirS; nitrite reductase (NO-forming) / hydroxylamine reductase [EC:1.7.2.1 1.7.99.1]
# K03385 nrfA; nitrite reductase (cytochrome c-552) [EC:1.7.2.2] (missing)
# K04561 norB; nitric oxide reductase subunit B [EC:1.7.2.5]
# K00376 nosZ; nitrous-oxide reductase [EC:1.7.2.4]

brady_mat_func <- brady_mat %>%
  filter(row.names(.) %in% c("K02111", "K02274", "K00404", "K02297", "K00425", "K00335", "K00239",
                             "K14658", "K14659", "K14666", "K09695", "K09694",
                             "K02585", "K02586", "K02587", "K02588", "K02591", "K02592",
                             "K00370", "K02567", "K00368", "K15864", "K04561", "K00376",
                             "K03520", "K10944", "K23995", 
                             "K17218", "K17230", "K17224", 
                             "K01601",
                             "K13991", "K08929", "K08928")) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(c("K02111", "K02274", "K00404", "K02297", "K00425", "K00335", "K00239",
                  "K14658", "K14659", "K14666", "K09695", "K09694",
                  "K02585", "K02586", "K02587", "K02588", "K02591", "K02592",
                  "K00370", "K02567", "K00368", "K15864", "K04561", "K00376",
                  "K03520", "K10944", "K23995", 
                  "K17218", "K17230", "K17224", 
                  "K01601",
                  "K13991", "K08929", "K08928")) %>%
  set_names(c("atpA", "coxA", "ccoN", "cyoA", "cydA", "nuoF", "sdhA",
              "nodA", "nodB", "nodC", "nodI", "nodJ", 
              "nifB", "nifD", "nifE", "nifH", "nifK", "nifN",
              "narG", "napA", "nirK", "nirS", "norB", "nosZ",
              "coxL", "pmoA-amoA", "xoxF",
              "sqr", "fccA", "soxB", 
              "rbcL",
              "puhA", "pufM", "pufL"))

# Run on an old brady_mat_func
# brady_func <- brady_mat_func %>%
#   mutate(Nod = rowSums(brady_mat_func[, c(1:5)])) %>%
#   mutate(Nfix = rowSums(brady_mat_func[, c(6:11)])) %>%
#   mutate(Photo = rowSums(brady_mat_func[, 12:14])) %>%
#   mutate(CO = brady_mat_func[, 15]) %>%
#   mutate(CH4 = brady_mat_func[, 16]) %>%
#   mutate(Meth = brady_mat_func[, 17]) %>%
#   mutate(DNF = rowSums(brady_mat_func[, 18:23])) %>%
#   mutate(DNF1 = rowSums(brady_mat_func[, c(19,20,22,23)])) %>%
#   rownames_to_column(var = "GenomeID") %>%
#   dplyr::select(GenomeID, Nod, Nfix, Photo, CO, CH4, Meth, DNF, DNF1)
# 
# # Groups
# sum(brady_func$Nod >= 4) # 143
# sum(brady_func$Nfix >= 5) # 145
# sum(brady_func$Photo >= 2) # 5
# sum(brady_func$CO > 0) # 181
# sum(brady_func$CH4 > 0) # 4
# sum(brady_func$Meth > 0) # 177
# sum(brady_func$DNF >= 4) # 13
# sum(brady_func$DNF1 == 4) # 11 with napA + nirK + norB + nosZ complete DNF!

# Individual
sum(brady_mat_func$nodA) # 143
sum(brady_mat_func$nifH) # 143
sum(brady_mat_func$puhA) # 5
sum(brady_mat_func$coxL) # 181
sum(brady_mat_func$`pmoA-amoA`) # 4
sum(brady_mat_func$xoxF) # 177
sum(brady_mat_func$narG) # 3
sum(brady_mat_func$napA) # 49
sum(brady_mat_func$nirK) # 142
sum(brady_mat_func$nirS) # 2
sum(brady_mat_func$norB) # 38
sum(brady_mat_func$nosZ) # 15

ann_cols <- data.frame(row.names = colnames(brady_mat_func),
                       Function = c(rep("Aerobic respiration", 7),
                                    rep("Nodulation", 5),
                                    rep("N fixation", 6),
                                    rep("Denitrification", 6),
                                    rep("C1 metabolism", 3),
                                    rep("Sulfur metabolism", 3),
                                    "Calvin-Benson cycle",
                                    rep("Photosynthesis", 3)))
sum(rownames(brady_mat_func) != rownames(sylph_strains_331))
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
                   `Function` = c("Aerobic respiration" = "#A6CEE3",
                                  "Nodulation" = "#1F78B4",
                                  "N fixation" = "#FB9A99",
                                  "Denitrification" = "#E31A1C",
                                  "C1 metabolism" = "#CAB2D6",
                                  "Sulfur metabolism" = "#6A3D9A",
                                  "Calvin-Benson cycle" = "#B15928",
                                  "Photosynthesis" = "#B2DF8A"
                                  ))
tree_data_maxmin <- tree_data_maxmin %>%
  mutate(SpeciesShort = gsub(" \\(Reference\\)", "", SpeciesShort))
hm <- pheatmap(brady_mat_func,
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
         gaps_col = c(7, 12, 18, 24, 27, 30, 31),
         labels_row = rev(tree_data_maxmin$SpeciesShort))
save_pheatmap_pdf(hm, "FinalFigs/Figure5forPPT.pdf")



#### __Greening ####
# See functions in Figure 1 Bay et al. 2021 Nature Microbiology
# K02111 ATPF1A, atpA; F-type H+/Na+-transporting ATPase subunit alpha [EC:7.1.2.2 7.2.2.1]
# K02274 coxA, ctaD; cytochrome c oxidase subunit I [EC:7.1.1.9]
# K00404 ccoN; cytochrome c oxidase cbb3-type subunit I [EC:7.1.1.9]
# K02297 cyoA; cytochrome o ubiquinol oxidase subunit II [EC:7.1.1.3]
# K00425 cydA; cytochrome bd ubiquinol oxidase subunit I [EC:7.1.1.7]
# K00335 nuoF; NADH-quinone oxidoreductase subunit F [EC:7.1.1.2]
# K00351 nqrF; Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
# K00239 sdhA, frdA; succinate dehydrogenase flavoprotein subunit [EC:1.3.5.1]
# K18556 frdA; NADH-dependent fumarate reductase subunit A [EC:1.3.1.6]
# K05299 fdhA; formate dehydrogenase (NADP+) alpha subunit [EC:1.17.1.10]
# K10944 pmoA-amoA; methane/ammonia monooxygenase subunit A [EC:1.14.18.3 1.14.99.39]
# K16157 mmoX; methane monooxygenase component A alpha chain [EC:1.14.13.25]
# K03520 coxL, cutL; aerobic carbon-monoxide dehydrogenase large subunit [EC:1.2.5.3]
# K18008 hydA; [NiFe] hydrogenase small subunit [EC:1.12.2.1]
# K17218 sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]
# K17230 fccA; cytochrome subunit of sulfide dehydrogenase
# K17224 soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]
# K11180 dsrA; dissimilatory sulfite reductase alpha subunit [EC:1.8.1.22]
# K16950 asrA; anaerobic sulfite reductase subunit A
# Skip NxrA
# K20932 K20932; hydrazine synthase subunit [EC:1.7.2.7]
# K00370 narG, narZ, nxrA; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]
# K02567 napA; nitrate reductase (cytochrome) [EC:1.9.6.1]
# K00368 nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]
# K15864 nirS; nitrite reductase (NO-forming) / hydroxylamine reductase [EC:1.7.2.1 1.7.99.1]
# K03385 nrfA; nitrite reductase (cytochrome c-552) [EC:1.7.2.2]
# K00376 nosZ; nitrous-oxide reductase [EC:1.7.2.4]
# K02588 nifH; nitrogenase iron protein NifH
# K14138 acsB; acetyl-CoA synthase [EC:2.3.1.169]
# K01601 rbcL, cbbL; ribulose-bisphosphate carboxylase large chain [EC:4.1.1.39]
# K15231 aclB; ATP-citrate lyase beta-subunit [EC:2.3.3.8]
# K14468 mcr; malonyl-CoA reductase / 3-hydroxypropionate dehydrogenase (NADP+) [EC:1.2.1.75 1.1.1.298]
# K14534 abfD; 4-hydroxybutyryl-CoA dehydratase / vinylacetyl-CoA-Delta-isomerase [EC:4.2.1.120 5.3.3.3]
# K18593 E6.2.1.56; 4-hydroxybutyrate---CoA ligase (ADP-forming) [EC:6.2.1.56]
# K02689 psaA; photosystem I P700 chlorophyll a apoprotein A1 [EC:1.97.1.12]
# K02703 psbA; photosystem II P680 reaction center D1 protein [EC:1.10.3.9]

# Not mentioned but check
# K10535 hao; hydroxylamine dehydrogenase [EC:1.7.2.6] (not present. Not AOB)
# Instead of the Bay photosystems, use the others.
bay <- data.frame("KO" = c("K02111", "K02274", "K00404", "K02297", "K00425", "K00335", 
                           "K00351", "K00239", "K18556", "K05299", "K10944", "K16157", "K03520", 
                           "K18008", "K17218", "K17230", "K17224", "K11180", "K16950", 
                           "K20932", "K00370", "K02567", "K00368", "K15864", 
                           "K03385", "K00376", "K02588", "K14138", "K01601", "K15231", 
                           "K14468", "K14534", "K18593", "K02689", "K02703", "K13991"),
                  "Name" = c("atpA", "coxA", "ccoN", "cyoA", "cydA", "nuoF", "nqrF",
                             "sdhA", "frdA", "fdhA", "pmoA-amoA", "mmoX", "coxL", "hydA",
                             "sqr", "fccA", "soxB", "dsrA", "asrA", "hzsA", "narG", 
                             "napA", "nirK", "nirS", "nrfA", "nosZ", "nifH", "acsB",
                             "rbcL", "aclB", "mcr", "abfD", "hbsT", "psaA", "psbA", "puhA"),
                  "Category" = c("Oxidative phosphorylation", "Aerobic respiration",
                                 "Aerobic respiration", "Aerobic respiration", 
                                 "Aerobic respiration", "NADH oxidation", "NADH oxidation",
                                 "Succinate oxidation", "Fumarate reduction",
                                 "Formate oxidation", "CH4 oxidation", "CH4 oxidation",
                                 "CO oxidation", "H2 oxidation", "Sulfide oxidation",
                                 "Sulfide oxidation", "Thiosulfate oxidation", 
                                 "Sulfide ox./Sulfite red.", "Sulfide ox./Sulfite red.", 
                                 "Anammox", "Nitrate reduction",
                                 "Nitrate reduction",
                                 "Nitrite reduction to NO",
                                 "Nitrite reduction to NO", "Nitrite reduction to NH3",
                                 "N2O reduction", "N fixation", "Wood-Ljungdahl Pathway",
                                 "Calvin-Benson cycle", "Reductive TCA cycle",
                                 "3-hydroxypropionate cycle", "4-hydroxybutyrate cycles",
                                 "4-hydroxybutyrate cycles", "Photosystem I", 
                                 "Photosystem II", "Photosynthesis")) %>%
  mutate(Presence = ifelse(KO %in% kos_merged$KO, "Yes", "No")) %>%
  filter(Presence == "Yes")

brady_mat_func <- brady_mat %>%
  filter(row.names(.) %in% bay$KO) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(bay$KO)) %>%
  set_names(c(bay$Name))

brady_cat <- data.frame(
  "Category" = c("Oxidative phosphorylation - AtpA", 
                 "Aerobic respiration - CoXA, CcoN, CyoA, CydA",
                 "NADH oxidation - NuoF",
                 "Succinate oxidation - SdhA",
                 "Methane oxidation - PmoA",
                 "Carbon monoxide oxidation - coxL",
                 "Sulfide oxidation - Sqr, FccA",
                 "Thiosulfate oxidation - SoxB",
                 "Nitrate reduction - NarG, NapA",
                 "Nitrite reduction to nitric oxide - NirS, NirK",
                 "Nitrous oxide reduction - NosZ",
                 "Nitrogen fixation - NifBDEHKN",
                 "Calvin-Benson cycle - RbcL",
                 "Photosynthesis - PuhA, PufM, PufL"),
  "GenomeCount" = c(sum(brady_mat_func$atpA),
                    sum(brady_mat_func$coxA == 1 | brady_mat_func$ccoN == 1 |
                          brady_mat_func$cyoA == 1 | brady_mat_func$cydA == 1),
                    sum(brady_mat_func$nuoF),
                    sum(brady_mat_func$sdhA),
                    sum(brady_mat_func$`pmoA-amoA`),
                    sum(brady_mat_func$coxL),
                    sum(brady_mat_func$sqr == 1 | brady_mat_func$fccA == 1),
                    sum(brady_mat_func$soxB),
                    sum(brady_mat_func$narG == 1 | brady_mat_func$napA == 1),
                    sum(brady_mat_func$nirK == 1 | brady_mat_func$nirS == 1),
                    sum(brady_mat_func$nosZ),
                    145,
                    sum(brady_mat_func$rbcL),
                    sum(brady_mat_func$puhA))
  )
brady_cat$Category <- factor(brady_cat$Category, levels = rev(brady_cat$Category))
ggplot(brady_cat, aes(Category, GenomeCount)) +
  geom_bar(stat = "identity") +
  geom_text(aes(Category, GenomeCount + 7, label = GenomeCount), size = 3) +
  labs(y = "# strains present") +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.margin = margin(5,15,5,5))



#### _Gene Clusters ####
# From anvio pangenomics analysis
# Pangenome summary file (anvi-summarize)
# Note: 79261 gene clusters with 1,416,019 genes
# Huge file (1.09 Gb), can't put on GitHub. It's on FigShare. 10.6084/m9.figshare.29175329
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
# pheatmap(cog_comb,
#          legend = T,
#          legend_breaks = c(0.001261655, 15, 30, 45, 60, 72.41355),
#          legend_labels = c("0", "15", "30", "45", "60", ""),
#          main = "            % Gene Clusters by COG",
#          cluster_rows = F,
#          cluster_cols = F,
#          labels_col = c("All (n = 79261)", 
#                         "In 181 genomes (n = 1935)",
#                         "In 2 to 180 genomes (n = 35541)",
#                         "In 1 genome (n = 41785)"),
#          angle_col = 315,
#          display_numbers = T,
#          number_color = "black",
#          fontsize_number = 10,
#          border_color = "white",
#          filename = "InitialFigs/Brady181_GeneClusters_COGs.png",
#          width = 6,
#          height = 8)
# dev.off()
# dev.set(dev.next())
# dev.set(dev.next())

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
# pheatmap(cog_comb,
#          legend = T,
#          cluster_rows = F,
#          cluster_cols = F,
#          gaps_col = 3,
#          labels_col = c("In 181 genomes (n = 1935)",
#                         "In 91 to 180 genomes (n = 4836)",
#                         "In 2 to 90 genomes (n = 30705)",
#                         "In 181 genomes (n = 1935)",
#                         "In 91 to 180 genomes (n = 4836)",
#                         "In 2 to 90 genomes (n = 30705)"),
#          annotation_col = ann_cols,
#          annotation_colors = ann_colors,
#          angle_col = 315,
#          display_numbers = T,
#          number_format = "%.3f",
#          number_color = "black",
#          fontsize_number = 6,
#          border_color = "white",
#          filename = "InitialFigs/Brady181_GeneClusters_Homogeneity.png",
#          width = 8,
#          height = 6)
# dev.off()
# dev.set(dev.next())
# dev.set(dev.next())



#### __Heatmap ####
# Which are most variable? P/A or geometric?
# P/A index - how many genomes?
# How many total GC per COG class?
# Make Figure S11

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
# pheatmap(cog_comb,
#          legend = T,
#          cluster_rows = F,
#          cluster_cols = F,
#          gaps_col = c(3, 6),
#          show_colnames = F,
#          annotation_col = ann_cols,
#          annotation_colors = ann_colors,
#          angle_col = 315,
#          display_numbers = cog_comb_num,
#          number_format = "%.2f",
#          number_color = "black",
#          fontsize = 7,
#          fontsize_number = 6,
#          border_color = "white",
#          filename = "InitialFigs/oldFigure7.png",
#          width = 8,
#          height = 6)
# dev.off()
# dev.set(dev.next())
# dev.set(dev.next())

# Plot number of GC shared by each number of genomes
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
gc_genomeNum <- gc %>%
  dplyr::select(gene_cluster_id, num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  group_by(num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  summarise(nGC = n())
length(unique(gc_genomeNum$num_genomes_gene_cluster_has_hits))
#pdf("InitialFigs/GeneClusters_NumGenomes.pdf", width = 8, height = 6)
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
#dev.off()

# Don't split by COG. Log10 y.
gc_genomeNum <- gc %>%
  dplyr::select(gene_cluster_id, num_genomes_gene_cluster_has_hits) %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(nGC = n()) %>%
  mutate(Type = ifelse(num_genomes_gene_cluster_has_hits == 1, "Singleton",
                       ifelse(num_genomes_gene_cluster_has_hits == 181, "Core", "Accessory"))) %>%
  mutate(Type = factor(Type, levels = c("Singleton", "Accessory", "Core")))
s10 <- ggplot(gc_genomeNum, aes(num_genomes_gene_cluster_has_hits, log10(nGC), fill = Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of strains (1 to 181)",
       y = "log number of gene clusters",
       fill = NULL) +
  scale_fill_manual(values = c("red", "grey50", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.01, 0.5)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        #axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,-10))
pdf("FinalFigs/FigureS10.pdf", width = 8, height = 6)
s10
dev.off()
png("FinalFigs/FigureS10.png", width = 8, height = 6, units = "in", res = 300)
s10
dev.off()

# Now plot proportions
gc_genomeNum <- gc %>%
  dplyr::select(gene_cluster_id, num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  group_by(num_genomes_gene_cluster_has_hits, COG20_Uniq) %>%
  summarise(nGC = n())
gcSum_genomeNum <- gc_genomeNum %>%
  group_by(num_genomes_gene_cluster_has_hits) %>%
  summarise(sum = sum(nGC)) %>%
  ungroup()
gcProp_genomeNum <- gc_genomeNum %>%
  left_join(., gcSum_genomeNum, by = "num_genomes_gene_cluster_has_hits") %>%
  mutate(Prop = nGC/sum)
#pdf("InitialFigs/GeneClustersProp_NumGenomes.pdf", width = 8, height = 6)
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
#dev.off()

# Or plot as panels to see the distribution of each one
#pdf("InitialFigs/GeneClustersProp_NumGenomes_Facet.pdf", width = 9, height = 6)
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
#dev.off()

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
#pdf("InitialFigs/GeneClustersPropCOG.pdf", width = 8, height = 6)
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
#dev.off()

cogProp_sing_core <- cogProp %>%
  filter(num_genomes_gene_cluster_has_hits == 1 | num_genomes_gene_cluster_has_hits == 181) %>%
  mutate(Type = ifelse(num_genomes_gene_cluster_has_hits == 1,
                       "Singleton", "Core")) %>%
  dplyr::select(COG, Prop, Type) %>%
  add_row(COG = "Chromatin structure and dynamics", Prop = 0, Type = "Core") %>%
  add_row(COG = "Cytoskeleton", Prop = 0, Type = "Core") %>%
  mutate(Type = factor(Type, levels = c("Singleton", "Core")))
#pdf("InitialFigs/GeneClustersPropCOG_SingCore.pdf", width = 8, height = 6)
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
#dev.off()

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



#### ___Figure S11 ####
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
         border_color = "white"
         # ,
         # filename = "FinalFigs/oldFigure7.png",
         # width = 8,
         # height = 6
         )
# dev.off()
# dev.set(dev.next())
# dev.set(dev.next())

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
figS11 <- ggplot() +
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
figS11
pdf("FinalFigs/FigureS11.pdf", width = 8, height = 6)
figS11
dev.off()
png("FinalFigs/FigureS11.png", width = 8, height = 6, units = "in", res = 300)
figS11
dev.off()



#### End ####