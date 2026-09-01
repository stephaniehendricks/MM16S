#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al

#https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#https://rpubs.com/mrgambero/lesson_20_ancom

#clear working environment
rm(list=ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC", force = TRUE)

#download packages
library(phyloseq)
library(dplyr)
library(ANCOMBC)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(DT)
library(cowplot) 

#set working directory
setwd("~/Desktop/StephV3V4/gzfiles/silva138.1rescript")

# FORMATTING DATA ----
# 1. making data into phyloseq object ----
physeq <- qza_to_phyloseq(
  features = "rep_seq_feature_table2.qza",
  tree = "rooted-tree.qza",
  taxonomy = "classified_taxonomy_table.qza",
  metadata = "MM16S_metadata.tsv") #this is the QIIME2 metadata file
physeq #otu_table()   OTU Table:         [ 3007 taxa and 69 samples ]

# 2. clean data ----
#extract metadata from the phyloseq object
metadata <- sample_data(physeq)
#remove samples (outliers)
samples_to_remove <- c("BS14_6dpf_L", "TS14_2dpf_L")
metadata <- metadata[!rownames(metadata) %in% samples_to_remove, ]
physeq <- prune_samples(!sample_names(physeq) %in% samples_to_remove, physeq)
#change names from sterile and filtered to LMR and HMR
metadata$type <- factor(
  ifelse(grepl("sterile", metadata$type), "LMR", "HMR"),
  levels = c("LMR", "HMR"))
#change name from water to seawater
metadata$medium <- as.character(metadata$medium)
metadata$medium[metadata$medium == "water"] <- "seawater"
#add ˚C to temp and treat
metadata$tempC <- paste0(metadata$temperature, "\u00B0C")  # \u00B0 = degree symbol
metadata$treatC <- paste(metadata$tempC, metadata$type, sep = " ")
#create a new treat column in metadata by combining column temperature and type in metadata
metadata$treat <- paste(metadata$temperature, metadata$type, sep = "")
metadata$all <- paste(metadata$day, metadata$temperature, metadata$type, metadata$medium, sep = ".")
#update metadata in phyloseq object
sample_data(physeq) <- sample_data(metadata)
#ensure column was added to metadata
(sample_data(physeq))
nsamples(physeq) #67 

#remove chloroplast, mitochondria, and eukaryota
clean <- subset_taxa(physeq, 
                    !Order %in% "Chloroplast" &
                    !Family %in% "Mitochondria" &
                    !grepl("Eukaryota", Kingdom) &
                    !is.na(Phylum) & Phylum != "")
clean #otu_table()   OTU Table:         [ 2832 taxa and 67 samples ]
physeq <- clean

# 3. SEAWATER SAMPLES ONLY at genus level ----
options(contrasts = c("contr.treatment", "contr.poly"))

pseq_genus_seawater <- physeq %>%
  subset_samples(medium == "seawater") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

metadata <- data.frame(sample_data(pseq_genus_seawater))

metadata$medium <- factor(metadata$medium, levels = c("seawater", "larval"))
metadata$type <- factor(metadata$type, levels = c("LMR", "HMR"))
metadata$temperature <- factor(metadata$temperature, levels = c(14, 18))

sample_data(pseq_genus_seawater) <- sample_data(metadata)

contrasts(sample_data(pseq_genus_seawater)$type)
contrasts(sample_data(pseq_genus_seawater)$temperature)

#testing for differential abundance of bacterial genera across temp and type 
out_genus_seawater <- ancombc(
  data = pseq_genus_seawater,
  tax_level = "Genus",
  formula = "temperature + type", #main effects 
  p_adj_method = "BH", #adjusts p-values using the Benjamini-Hochberg FDR method
  prv_cut = 0.10, #removes taxa present in less than 10% of samples
  lib_cut = 1000, #removes samples with library size < 1000 reads
  group = NULL, #main comparison, which variable to perform pairwise comparisons on
  struc_zero = FALSE, #detects structural zeros (taxa that are truly absent in a group)
  neg_lb = TRUE, #ensures negative lower bounds for confidence intervals are allowed
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05, #significance threshold for declaring taxa differentially abundant
  global = FALSE, #calculates a global test statistic across all treatment levels, not just pairwise comparisons
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus_seawater)$temperature)
table(sample_data(pseq_genus_seawater)$type)

#extract results
res_genus_seawater <- out_genus_seawater$res
head(res_genus_seawater)

#log fold changes table
tab_lfc <- res_genus_seawater$lfc
tab_lfc %>% datatable(caption = "Log Fold Changes")
head(tab_lfc)

#standard errors
tab_se <- res_genus_seawater$se
tab_se %>% datatable(caption = "SEs") 
head(tab_se)

#test statistics
tab_w <- res_genus_seawater$W
tab_w %>% datatable(caption = "Test Statistics") 
head(tab_w)

#p-values
tab_p <- res_genus_seawater$p_val
tab_p %>% datatable(caption = "P-values")
head(tab_p)

#adjusted p-values
tab_q <- res_genus_seawater$q
tab_q %>% datatable(caption = "Adjusted p-values") 
head(tab_q)

#differentially abundant taxa
tab_diff <- res_genus_seawater$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa")
head(tab_diff)

# → visualization for type ----
res_genus_seawater <- out_genus_seawater$res

df_fig_type <- data.frame(
  taxon_id = res_genus_seawater$lfc$taxon,
  typeHMR = as.numeric(res_genus_seawater$lfc$typeHMR),
  typeHMRSE = as.numeric(res_genus_seawater$se$typeHMR),
  typeHMR_q = as.numeric(res_genus_seawater$q$typeHMR),
  typeHMR_diff = as.logical(res_genus_seawater$diff_abn$typeHMR))

#determine direction of the type effect
#positive LFC = higher in HMR
#negative LFC = higher in LMR
df_fig_type <- df_fig_type %>%
  mutate(type_direction = ifelse(
    typeHMR > 0, "HMR", "LMR"))
head(df_fig_type)

#significant taxa
df_fig_type %>%
  filter(typeHMR_diff == TRUE) %>%
  select(
    taxon_id,
    typeHMR,
    typeHMRSE,
    typeHMR_q,
    type_direction) %>%
  arrange(typeHMR_q)

#top 5 pos and neg LFC
#top 5 significant taxa enriched in HMR samples
top_pos <- df_fig_type %>%
  filter(typeHMR_diff == TRUE, typeHMR > 0) %>%
  slice_max(order_by = typeHMR, n = 5)

#top 5 significant taxa enriched in LMR samples
top_neg <- df_fig_type %>%
  filter(typeHMR_diff == TRUE, typeHMR < 0) %>%
  slice_min(order_by = typeHMR, n = 5)

#combine and order for plotting
top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(type_direction = factor(type_direction, levels = c("LMR", "HMR"))) %>%
  arrange(typeHMR) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#check which taxa are being plotted
top_pos_neg %>%
  select(
    taxon_id,
    typeHMR,
    typeHMRSE,
    typeHMR_q,
    type_direction)

res_genus_seawater$lfc %>% filter(taxon == "Lentibacter")
res_genus_seawater$q %>% filter(taxon == "Lentibacter")
res_genus_seawater$diff_abn %>% filter(taxon == "Lentibacter")
res_genus_seawater$se %>% filter(taxon == "Lentibacter")

res_genus_seawater$lfc %>% filter(taxon == "Ekhidna")
res_genus_seawater$q %>% filter(taxon == "Ekhidna")
res_genus_seawater$diff_abn %>% filter(taxon == "Ekhidna")
res_genus_seawater$se %>% filter(taxon == "Ekhidna")

res_genus_seawater$lfc %>% filter(taxon == "Psychrobium")
res_genus_seawater$q %>% filter(taxon == "Psychrobium")
res_genus_seawater$diff_abn %>% filter(taxon == "Psychrobium")
res_genus_seawater$se %>% filter(taxon == "Psychrobium")

res_genus_seawater$lfc %>% filter(taxon == "Bermanella")
res_genus_seawater$q %>% filter(taxon == "Bermanella")
res_genus_seawater$diff_abn %>% filter(taxon == "Bermanella")
res_genus_seawater$se %>% filter(taxon == "Bermanella")

#plot aesthetics
colorsW <- c("Lentibacter" = "#BEE3BA", #rhodo
             "uncultured_60" = "#CCCCCC", #other
             "Mycobacterium" = "#CCCCCC", #other
             "Peregrinibacteria" = "#CCCCCC", #other
             "Vicingus" = "#CCCCCC", #other
             "Pseudophaeobacter" = "#BEE3BA", #rhodo
             "Ekhidna" = "#CCCCCC", #other
             "Psychrobium" = "#FFB347", #shew
             "Aestuariibacter" = "#6495ED", #altero
             "Bermanella" = "#CCCCCC") #other

top_pos_neg$taxon_id <- as.character(top_pos_neg$taxon_id)
top_pos_neg$taxon_id[top_pos_neg$taxon_id == "Candidatus_Peregrinibacteria"] <- "Peregrinibacteria"

top_pos_neg <- top_pos_neg %>%
  arrange(typeHMR, taxon_id) %>%   
  mutate(taxon_id = factor(taxon_id, levels = unique(taxon_id)))

# → plot ----
quartz()
p_type_all_sw_plot <- ggplot(top_pos_neg, aes(x = taxon_id, y = typeHMR, fill = taxon_id)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = typeHMR - typeHMRSE, ymax = typeHMR + typeHMRSE),
    width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
    axis.title = element_text(size = 11, colour = "black"),
    legend.position = "none") +
  coord_flip()
p_type_all_sw_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_type_all_sw_plot.jpg", plot = p_type_all_sw_plot, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
#ggsave("p_type_all_sw_plot.svg", plot = p_type_all_sw_plot, scale = 1, width = 3, height = 2.5,
#       units = "in", dpi = 300)

# → visualization for temp ----
df_fig_temp <- data.frame(
  taxon_id = res_genus_seawater$lfc$taxon,
  temperature18 = as.numeric(res_genus_seawater$lfc$temperature18),
  temperature18SE = as.numeric(res_genus_seawater$se$temperature18),
  temperature18_q = as.numeric(res_genus_seawater$q$temperature18),
  temperature18_diff = as.logical(res_genus_seawater$diff_abn$temperature18))

#determine direction of the temp effect
#positive LFC = higher in 18
#negative LFC = higher in 14
df_fig_temp <- df_fig_temp %>%
  mutate(temp_direction = ifelse(
    temperature18 > 0, "18", "14"))
head(df_fig_temp)

#significant taxa
df_fig_temp %>%
  filter(temperature18_diff == TRUE) %>%
  select(
    taxon_id,
    temperature18,
    temperature18SE,
    temperature18_q,
    temp_direction) %>%
  arrange(temperature18_q)

#top 5 pos and neg LFC
#top 5 significant taxa enriched in 18 samples
top_pos <- df_fig_temp %>%
  filter(temperature18_diff == TRUE, temperature18 > 0) %>%
  slice_max(order_by = temperature18, n = 5)

#top 5 significant taxa enriched in 14 samples
top_neg <- df_fig_temp %>%
  filter(temperature18_diff == TRUE, temperature18 < 0) %>%
  slice_min(order_by = temperature18, n = 5)

#combine and order for plotting
top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(temp_direction = factor(temp_direction, levels = c("14", "18"))) %>%
  arrange(temperature18) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#check which taxa are being plotted
top_pos_neg %>%
  select(
    taxon_id,
    temperature18,
    temperature18SE,
    temperature18_q,
    temp_direction)

#stop plotting bc there is only 1 result 

res_genus_seawater$lfc %>% filter(taxon == "Nautella")
res_genus_seawater$q %>% filter(taxon == "Nautella")
res_genus_seawater$diff_abn %>% filter(taxon == "Nautella")
res_genus_seawater$se %>% filter(taxon == "Nautella")

# 4. LARVAL SAMPLES ONLY at genus level ----
options(contrasts = c("contr.treatment", "contr.poly"))

pseq_genus_larval <- physeq %>%
  subset_samples(medium == "larval") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

metadata <- data.frame(sample_data(pseq_genus_larval))

metadata$medium <- factor(metadata$medium, levels = c("seawater", "larval"))
metadata$type <- factor(metadata$type, levels = c("LMR", "HMR"))
metadata$temperature <- factor(metadata$temperature, levels = c(14, 18))

sample_data(pseq_genus_larval) <- sample_data(metadata)

contrasts(sample_data(pseq_genus_larval)$type)
contrasts(sample_data(pseq_genus_larval)$temperature)

#testing for differential abundance of bacterial genera across temp and type 
out_genus_larval <- ancombc(
  data = pseq_genus_larval,
  tax_level = "Genus",
  formula = "temperature + type", 
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = NULL, 
  struc_zero = FALSE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = FALSE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus_larval)$temperature)
table(sample_data(pseq_genus_larval)$type)

#extract results
res_genus_larval <- out_genus_larval$res
head(res_genus_larval)

#log fold changes table
tab_lfc <- res_genus_larval$lfc
tab_lfc %>% datatable(caption = "Log Fold Changes")
head(tab_lfc)

#standard errors
tab_se <- res_genus_larval$se
tab_se %>% datatable(caption = "SEs") 
head(tab_se)

#test statistics
tab_w <- res_genus_larval$W
tab_w %>% datatable(caption = "Test Statistics") 
head(tab_w)

#p-values
tab_p <- res_genus_larval$p_val
tab_p %>% datatable(caption = "P-values") 
head(tab_p)

#adjusted p-values
tab_q <- res_genus_larval$q
tab_q %>% datatable(caption = "Adjusted p-values") 
head(tab_q)

#differentially abundant taxa
tab_diff <- res_genus_larval$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa")
head(tab_diff)

# → visualization for type ----
df_fig_type <- data.frame(
  taxon_id = res_genus_larval$lfc$taxon,
  typeHMR = as.numeric(res_genus_larval$lfc$typeHMR),
  typeHMRSE = as.numeric(res_genus_larval$se$typeHMR),
  typeHMR_q = as.numeric(res_genus_larval$q$typeHMR),
  typeHMR_diff = as.logical(res_genus_larval$diff_abn$typeHMR))

#determine direction of the type effect
#positive LFC = higher in HMR
#negative LFC = higher in LMR
df_fig_type <- df_fig_type %>%
  mutate(type_direction = ifelse(
    typeHMR > 0, "HMR", "LMR"))
head(df_fig_type)

#significant taxa
df_fig_type %>%
  filter(typeHMR_diff == TRUE) %>%
  select(
    taxon_id,
    typeHMR,
    typeHMRSE,
    typeHMR_q,
    type_direction) %>%
  arrange(typeHMR_q)

#top 5 pos and neg LFC
#top 5 significant taxa enriched in HMR samples
top_pos <- df_fig_type %>%
  filter(typeHMR_diff == TRUE, typeHMR > 0) %>%
  slice_max(order_by = typeHMR, n = 5)

#top 5 significant taxa enriched in LMR samples
top_neg <- df_fig_type %>%
  filter(typeHMR_diff == TRUE, typeHMR < 0) %>%
  slice_min(order_by = typeHMR, n = 5)

#combine and order for plotting
top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(type_direction = factor(type_direction, levels = c("LMR", "HMR"))) %>%
  arrange(typeHMR) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#check which taxa are being plotted
top_pos_neg %>%
  select(
    taxon_id,
    typeHMR,
    typeHMRSE,
    typeHMR_q,
    type_direction)

res_genus_larval$lfc %>% filter(taxon == "Flavicella")
res_genus_larval$q %>% filter(taxon == "Flavicella")
res_genus_larval$diff_abn %>% filter(taxon == "Flavicella")
res_genus_larval$se %>% filter(taxon == "Flavicella")

res_genus_larval$lfc %>% filter(taxon == "Ekhidna")
res_genus_larval$q %>% filter(taxon == "Ekhidna")
res_genus_larval$diff_abn %>% filter(taxon == "Ekhidna")
res_genus_larval$se %>% filter(taxon == "Ekhidna")

res_genus_larval$lfc %>% filter(taxon == "Psychrobium")
res_genus_larval$q %>% filter(taxon == "Psychrobium")
res_genus_larval$diff_abn %>% filter(taxon == "Psychrobium")
res_genus_larval$se %>% filter(taxon == "Psychrobium")

res_genus_larval$lfc %>% filter(taxon == "Bermanella")
res_genus_larval$q %>% filter(taxon == "Bermanella")
res_genus_larval$diff_abn %>% filter(taxon == "Bermanella")
res_genus_larval$se %>% filter(taxon == "Bermanella")

#plot aesthetics
colorsW <- c("Flavicella" = "#FFD1DF", #flavo
             "Marinicella" = "#CCCCCC", #other
             "uncultured_15" = "#CCCCCC", #other
             "uncultured_11" = "#CCCCCC", #other
             "Marinagarivorans" = "#CCCCCC", #other
             "Pleionea" = "#CCCCCC", #other
             "Fulvivirga" = "#CCCCCC", #other
             "Psychrobium" = "#FFB347", #shew
             "Ekhidna" = "#CCCCCC", #other
             "Bermanella" = "#CCCCCC") #other

top_pos_neg$taxon_id <- as.character(top_pos_neg$taxon_id)

top_pos_neg <- top_pos_neg %>%
  arrange(typeHMR, taxon_id) %>%   
  mutate(taxon_id = factor(taxon_id, levels = unique(taxon_id)))

# → plot ----
quartz()
p_type_all_larval_plot <- ggplot(top_pos_neg, aes(x = taxon_id, y = typeHMR, fill = taxon_id)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = typeHMR - typeHMRSE, ymax = typeHMR + typeHMRSE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip()
p_type_all_larval_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_type_all_larval_plot.jpg", plot = p_type_all_larval_plot, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
#ggsave("p_type_all_larval_plot.svg", plot = p_type_all_larval_plot, scale = 1, width = 3, height = 2.5,
#       units = "in", dpi = 300)

# → visualization for temp ----
df_fig_temp <- data.frame(
  taxon_id = res_genus_larval$lfc$taxon,
  temperature18 = as.numeric(res_genus_larval$lfc$temperature18),
  temperature18SE = as.numeric(res_genus_larval$se$temperature18),
  temperature18_q = as.numeric(res_genus_larval$q$temperature18),
  temperature18_diff = as.logical(res_genus_larval$diff_abn$temperature18))

#determine direction of the temp effect
#positive LFC = higher in 18
#negative LFC = higher in 14
df_fig_temp <- df_fig_temp %>%
  mutate(temp_direction = ifelse(
    temperature18 > 0, "18", "14"))
head(df_fig_temp)

#significant taxa
df_fig_temp %>%
  filter(temperature18_diff == TRUE) %>%
  select(
    taxon_id,
    temperature18,
    temperature18SE,
    temperature18_q,
    temp_direction) %>%
  arrange(temperature18_q)

#top 5 pos and neg LFC
#top 5 significant taxa enriched in 18 samples
top_pos <- df_fig_temp %>%
  filter(temperature18_diff == TRUE, temperature18 > 0) %>%
  slice_max(order_by = temperature18, n = 5)

#top 5 significant taxa enriched in 14 samples
top_neg <- df_fig_temp %>%
  filter(temperature18_diff == TRUE, temperature18 < 0) %>%
  slice_min(order_by = temperature18, n = 5)

#combine and order for plotting
top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(temp_direction = factor(temp_direction, levels = c("14", "18"))) %>%
  arrange(temperature18) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#check which taxa are being plotted
top_pos_neg %>%
  select(
    taxon_id,
    temperature18,
    temperature18SE,
    temperature18_q,
    temp_direction)

res_genus_larval$lfc %>% filter(taxon == "Neisseria")
res_genus_larval$q %>% filter(taxon == "Neisseria")
res_genus_larval$diff_abn %>% filter(taxon == "Neisseria")
res_genus_larval$se %>% filter(taxon == "Neisseria")

res_genus_larval$lfc %>% filter(taxon == "Nautella")
res_genus_larval$q %>% filter(taxon == "Nautella")
res_genus_larval$diff_abn %>% filter(taxon == "Nautella")
res_genus_larval$se %>% filter(taxon == "Nautella")

#plot aesthetics
colorsW <- c("Prevotella_7" = "#CCCCCC", #other
             "Haemophilus" = "#CCCCCC", #other
             "Neisseria" = "#CCCCCC", #other     
             "Porphyromonas" = "#CCCCCC", #other
             "uncultured_16" = "#CCCCCC", #other
             "Marivita" = "#BEE3BA", #rhodo
             "uncultured_60" = "#CCCCCC", #other
             "Thalassospira" = "#CCCCCC", #other
             "Nautella" = "#BEE3BA") #rhodo

top_pos_neg$taxon_id <- as.character(top_pos_neg$taxon_id)

top_pos_neg <- top_pos_neg %>%
  arrange(temperature18, taxon_id) %>%   
  mutate(taxon_id = factor(taxon_id, levels = unique(taxon_id)))

# → plot ----
quartz()
p_temp_all_larval_plot <- ggplot(top_pos_neg, aes(x = taxon_id, y = temperature18, fill = taxon_id)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = temperature18 - temperature18SE, ymax = temperature18 + temperature18SE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip()
p_temp_all_larval_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_temp_all_larval_plot.jpg", plot = p_temp_all_larval_plot, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
#ggsave("p_temp_all_larval_plot.svg", plot = p_temp_all_larval_plot, scale = 1, width = 3, height = 2.5,
#       units = "in", dpi = 300)

# 5. SEAWATER VS LARVAL SAMPLES (medium) at family level ----
options(contrasts = c("contr.treatment", "contr.poly"))

pseq_family_medium_6dpf <- physeq %>%
  subset_taxa(!is.na(Family) & Family != "") %>%
  subset_samples(day == "6")

metadata <- data.frame(sample_data(pseq_family_medium_6dpf))

metadata$medium <- factor(metadata$medium, levels = c("seawater", "larval"))
metadata$type <- factor(metadata$type, levels = c("LMR", "HMR"))
metadata$temperature <- factor(metadata$temperature, levels = c(14, 18))

sample_data(pseq_family_medium_6dpf) <- sample_data(metadata)

contrasts(sample_data(pseq_family_medium_6dpf)$medium)
contrasts(sample_data(pseq_family_medium_6dpf)$type)
contrasts(sample_data(pseq_family_medium_6dpf)$temperature)

#testing for differences between larval and seawater while controlling for temperature and type
out_family_medium <- ancombc(
  data = pseq_family_medium_6dpf,
  tax_level = "Family",
  formula = "medium + temperature + type", #main effects         
  p_adj_method = "BH", #adjusts p-values using the Benjamini-Hochberg FDR method
  prv_cut = 0.10, #removes taxa present in less than 10% of samples
  lib_cut = 1000, #removes samples with library size < 1000 reads
  group = NULL, #main comparison, which variable to perform pairwise comparisons on, only need to state when doing global = TRUE
  struc_zero = FALSE, #detects structural zeros (taxa that are truly absent in a group)
  neg_lb = TRUE, #ensures negative lower bounds for confidence intervals are allowed
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05, #significance threshold for declaring taxa differentially abundant
  global = FALSE, #calculates a global test statistic across all treatment levels, not just pairwise comparisons
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_family_medium_6dpf)$medium)

#medium1       = larval vs seawater
#type1         = HMR vs LMR
#temperature1  = 18°C vs 14°C

#LFC directions
#medium1:
#+ → higher in larval
#− → higher in seawater

#type1:
#+ → higher in HMR
#− → higher in LMR

#temperature1:
#+ → higher at 18°C
#− → higher at 14°C

res_family_medium <- out_family_medium$res
head(res_family_medium)

#log fold changes table
tab_lfc <- res_family_medium$lfc
tab_lfc %>% datatable(caption = "Log Fold Changes")
head(tab_lfc)

#standard errors
tab_se <- res_family_medium$se
tab_se %>% datatable(caption = "SEs") 
head(tab_se)

#test statistics
tab_w <- res_family_medium$W
tab_w %>% datatable(caption = "Test Statistics") 
head(tab_w)

#p-values
tab_p <- res_family_medium$p_val
tab_p %>% datatable(caption = "P-values")
head(tab_p)

#adjusted p-values
tab_q <- res_family_medium$q
tab_q %>% datatable(caption = "Adjusted p-values") 
head(tab_q)

#differentially abundant taxa
tab_diff <- res_family_medium$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa")
head(tab_diff)

res_family_medium$lfc %>% filter(taxon == "Alteromonadaceae")
res_family_medium$q %>% filter(taxon == "Alteromonadaceae")
res_family_medium$diff_abn %>% filter(taxon == "Alteromonadaceae")
res_family_medium$se %>% filter(taxon == "Alteromonadaceae")

res_family_medium$lfc %>% filter(taxon == "Flavobacteriaceae")
res_family_medium$q %>% filter(taxon == "Flavobacteriaceae")
res_family_medium$diff_abn %>% filter(taxon == "Flavobacteriaceae")
res_family_medium$se %>% filter(taxon == "Flavobacteriaceae")

res_family_medium$lfc %>% filter(taxon == "Pseudoalteromonadaceae")
res_family_medium$q %>% filter(taxon == "Pseudoalteromonadaceae")
res_family_medium$diff_abn %>% filter(taxon == "Pseudoalteromonadaceae")
res_family_medium$se %>% filter(taxon == "Pseudoalteromonadaceae")

res_family_medium$lfc %>% filter(taxon == "Rhodobacteraceae")
res_family_medium$q %>% filter(taxon == "Rhodobacteraceae")
res_family_medium$diff_abn %>% filter(taxon == "Rhodobacteraceae")
res_family_medium$se %>% filter(taxon == "Rhodobacteraceae")

res_family_medium$lfc %>% filter(taxon == "Saccharospirillaceae")
res_family_medium$q %>% filter(taxon == "Saccharospirillaceae")
res_family_medium$diff_abn %>% filter(taxon == "Saccharospirillaceae")
res_family_medium$se %>% filter(taxon == "Saccharospirillaceae")

res_family_medium$lfc %>% filter(taxon == "Shewanellaceae")
res_family_medium$q %>% filter(taxon == "Shewanellaceae")
res_family_medium$diff_abn %>% filter(taxon == "Shewanellaceae")
res_family_medium$se %>% filter(taxon == "Shewanellaceae")

res_family_medium$lfc %>% filter(taxon == "Vibrionaceae")
res_family_medium$q %>% filter(taxon == "Vibrionaceae")
res_family_medium$diff_abn %>% filter(taxon == "Vibrionaceae")
res_family_medium$se %>% filter(taxon == "Vibrionaceae")

# 6. SEAWATER VS LARVAL SAMPLES (medium) at genus level ----
options(contrasts = c("contr.treatment", "contr.poly"))

pseq_genus_medium_6dpf <- physeq %>%
  subset_taxa(!is.na(Genus) & Genus != "") %>%
  subset_samples(day == "6")

metadata <- data.frame(sample_data(pseq_genus_medium_6dpf))

metadata$medium <- factor(metadata$medium, levels = c("seawater", "larval"))
metadata$type <- factor(metadata$type, levels = c("LMR", "HMR"))
metadata$temperature <- factor(metadata$temperature, levels = c(14, 18))

sample_data(pseq_genus_medium_6dpf) <- sample_data(metadata)

contrasts(sample_data(pseq_genus_medium_6dpf)$medium)
contrasts(sample_data(pseq_genus_medium_6dpf)$type)
contrasts(sample_data(pseq_genus_medium_6dpf)$temperature)

#testing for differences between larval and seawater while controlling for temperature and type
out_genus_medium <- ancombc(
  data = pseq_genus_medium_6dpf,
  tax_level = "Genus",
  formula = "medium + temperature + type",     
  p_adj_method = "BH", 
  prv_cut = 0.10, 
  lib_cut = 1000, 
  group = NULL, 
  struc_zero = FALSE, 
  neg_lb = TRUE, 
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05, 
  global = FALSE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus_medium_6dpf)$medium)

#medium1       = larval vs seawater
#type1         = HMR vs LMR
#temperature1  = 18°C vs 14°C

#LFC directions
#medium1:
#+ → higher in larval
#− → higher in seawater

#type1:
#+ → higher in HMR
#− → higher in LMR

#temperature1:
#+ → higher at 18°C
#− → higher at 14°C

res_genus_medium <- out_genus_medium$res
head(res_genus_medium)

#log fold changes table
tab_lfc <- res_genus_medium$lfc
tab_lfc %>% datatable(caption = "Log Fold Changes")
head(tab_lfc)

#standard errors
tab_se <- res_genus_medium$se
tab_se %>% datatable(caption = "SEs") 
head(tab_se)

#test statistics
tab_w <- res_genus_medium$W
tab_w %>% datatable(caption = "Test Statistics") 
head(tab_w)

#p-values
tab_p <- res_genus_medium$p_val
tab_p %>% datatable(caption = "P-values")
head(tab_p)

#adjusted p-values
tab_q <- res_genus_medium$q
tab_q %>% datatable(caption = "Adjusted p-values") 
head(tab_q)

#differentially abundant taxa
tab_diff <- res_genus_medium$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa")
head(tab_diff)

# → visualization for medium ----
df_fig_medium <- data.frame(
  taxon_id = res_genus_medium$lfc$taxon,
  mediumlarval = as.numeric(res_genus_medium$lfc$mediumlarval),
  mediumlarvalSE = as.numeric(res_genus_medium$se$mediumlarval),
  mediumlarval_q = as.numeric(res_genus_medium$q$mediumlarval),
  mediumlarval_diff = as.logical(res_genus_medium$diff_abn$mediumlarval))

#determine direction of medium
#positive LFC = higher in larval
#negative LFC = higher in seawater
df_fig_medium <- df_fig_medium %>%
  mutate(medium_direction = ifelse(
    mediumlarval > 0, "larval", "seawater"))
head(df_fig_medium)

#significant taxa
df_fig_medium %>%
  filter(mediumlarval_diff == TRUE) %>%
  select(
    taxon_id,
    mediumlarval,
    mediumlarvalSE,
    mediumlarval_q,
    medium_direction) %>%
  arrange(mediumlarval_q)

#top 5 pos and neg LFC
#top 5 significant taxa enriched in larval samples
top_pos <- df_fig_medium %>%
  filter(mediumlarval_diff == TRUE, mediumlarval > 0) %>%
  slice_max(order_by = mediumlarval, n = 5)

#top 5 significant taxa enriched in seawater samples
top_neg <- df_fig_medium %>%
  filter(mediumlarval_diff == TRUE, mediumlarval < 0) %>%
  slice_min(order_by = mediumlarval, n = 5)

#combine and order for plotting
top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(medium_direction = factor(medium_direction, levels = c("seawater", "larval"))) %>%
  arrange(mediumlarval) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#check which taxa are being plotted
top_pos_neg %>%
  select(
    taxon_id,
    mediumlarval,
    mediumlarvalSE,
    mediumlarval_q,
    medium_direction)

#plot aesthetics
colorsW <- c(
  "Aureispira" = "#CCCCCC", #other
  "Marinicella" = "#CCCCCC", #other
  "Reichenbachiella" = "#CCCCCC", #other
  "Lewinella" = "#CCCCCC", #other
  "uncultured_9" = "#CCCCCC", #other
  "Marivita" = "#BEE3BA", #rhodo
  "uncultured_60" = "#CCCCCC", #other
  "Lentibacter" = "#BEE3BA", #rhodo
  "Maricaulis" = "#CCCCCC", #other
  "Flavicella" = "#FFD1DF" #flavo
)

# → plot ----
quartz()
p_medium_plot <- ggplot(top_pos_neg, aes(x = taxon_id, y = mediumlarval, fill = taxon_id)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mediumlarval - mediumlarvalSE, ymax = mediumlarval + mediumlarvalSE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  scale_y_continuous(breaks = c(-3, 0, 3), limits = c(-5, 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip()
p_medium_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_medium_plot.jpg", plot = p_medium_plot, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
#ggsave("p_medium_plot.svg", plot = p_medium_plot, scale = 1, width = 4, height = 4.25,
#       units = "in", dpi = 300)

############## fin
