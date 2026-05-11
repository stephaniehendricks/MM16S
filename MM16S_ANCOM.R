#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al

#https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#https://rpubs.com/mrgambero/lesson_20_ancom

#clear working environment
rm(list=ls())

knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)

options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC", force = TRUE)

remove.packages("readxl")
install.packages("readxl", type = "source")
remove.packages(c("haven", "readxl"))
install.packages(c("haven", "readxl"), type = "source")

#download packages
library(phyloseq)
library(dplyr)
library(ANCOMBC)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(DT)
library(cowplot) 

#use MM16S_manuscript.R to load in ps.rarefied
# TYPE ----
# 1. seawater by type ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "seawater") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

metadata$type <- factor(metadata$type, levels = c("LMR", "HMR"))  
levels(metadata$type)
metadata$tempC <- factor(paste0(metadata$temperature, "\u00B0C")) 
metadata$medium <- factor(metadata$medium)
metadata$medium <- factor(metadata$medium, levels = c("seawater", "larval"))
metadata$treat <- paste(metadata$temperature, metadata$type, sep = "") 
metadata$treat <- factor(metadata$treat)
metadata$treat <- factor(metadata$treat, levels = c("14LMR", "14HMR", "18LMR", "18HMR"))
levels(metadata$treat)

sample_data(pseq_genus) <- sample_data(metadata)
metadata <- sample_data(pseq_genus)
unique(metadata$medium)

#testing for differential abundance of bacterial families across medium 
out <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "temperature + type", #main effects + interaction
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "type", #main comparison
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus)$temperature)
table(sample_data(pseq_genus)$type)

#extract results
res <- out$res
res_global <- out$res_global

#log fold changes table
tab_lfc <- res$lfc
colnames(tab_lfc)
#add taxon names 
col_name = c("Taxon", "Intercept", "type", "medium")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#standard errors
tab_se <- res$se
colnames(tab_se) <- col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#p-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#bias-corrected log abundances
samp_frac = out$samp_frac
#replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
#add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
#adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
#show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% datatable(caption = "Bias-corrected log observed abundances")

#visualization (bar and heatmap)
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

colnames(df_se)[4] <- "typeHMRSE"

df_fig_type <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `typeHMR`, `typeHMRSE`) %>%  #LFC + SE
  dplyr::filter(`typeHMR` != 0) %>%
  dplyr::arrange(desc(`typeHMR`)) %>%
  dplyr::mutate(direct = ifelse(`typeHMR` > 0, "Positive LFC", "Negative LFC"))

df_fig_type$taxon_id <- factor(df_fig_type$taxon_id, levels = df_fig_type$taxon_id)
df_fig_type$direct <- factor(df_fig_type$direct, levels = c("Positive LFC", "Negative LFC"))

p_type_sw <- ggplot(df_fig_type, aes(x = taxon_id, y = `typeHMR`, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `typeHMR` - `typeHMRSE` , ymax = `typeHMR` + `typeHMRSE`),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log Fold Change", title = "HMR vs LMR Microbiome") +
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_type_sw

#split top positive and negative LFC
top_pos <- df_fig_type %>% filter(typeHMR > 0) %>% slice_max(typeHMR, n = 5)
top_neg <- df_fig_type %>% filter(typeHMR < 0) %>% slice_min(typeHMR, n = 5)
top_pos_neg <- bind_rows(top_pos, top_neg)

#reorder for plotting 
top_pos_neg <- top_pos_neg %>%
  arrange(typeHMR) %>%   #sorts negative to positive
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(type_direction = ifelse(typeHMR > 0, "HMR", "LMR")) %>%
  arrange(typeHMR) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id)) %>%
  mutate(type_direction = factor(type_direction, levels = c("LMR", "HMR")))

#bars extending to the right (positive) - taxa enriched in HMR
#bars extending to the left (negative) - taxa enriched in LMR
levels(sample_data(pseq_genus)$type) #LMR is reference, HMR is comparison

colorsW <- c("Lentibacter" = "#BEE3BA", #rhodo
             "uncultured_60" = "#CCCCCC", #other
             "Mycobacterium" = "#CCCCCC", #other
             "Peregrinibacteria" = "#CCCCCC", #other
             "Aliivibrio" = "#CDC0F8", #vibrio
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

## → plot p_type_all_sw ----
p_type_all_sw <- ggplot(top_pos_neg, aes(x = taxon_id, y = typeHMR, fill = taxon_id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = typeHMR - typeHMRSE, ymax = typeHMR + typeHMRSE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip()
p_type_all_sw
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_type_all_sw_plot.jpg", plot = p_type_all_sw, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
ggsave("p_type_all_sw_plot.svg", plot = p_type_all_sw, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)

# 2. larval by type ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "larval") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

metadata <- sample_data(pseq_genus)
unique(metadata$medium)

#testing for differential abundance of bacterial families across temp and type 
out <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "temperature + type", #main effects + interaction
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "type", #main comparison
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus)$temperature)
table(sample_data(pseq_genus)$type)

#extract results
res <- out$res
res_global <- out$res_global

#log fold changes table
tab_lfc <- res$lfc
colnames(tab_lfc)

#add taxon names 
col_name = c("Taxon", "Intercept", "type", "medium")
colnames(tab_lfc) = col_name
tab_lfc %>% datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#standard errors
tab_se <- res$se
colnames(tab_se) <- col_name
tab_se %>% datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#p-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#bias-corrected log abundances
samp_frac = out$samp_frac
#replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
#add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
#adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
#show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% datatable(caption = "Bias-corrected log observed abundances")

#visualization
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

colnames(df_se)[4] <- "typeHMRSE"

df_fig_type <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `typeHMR`, `typeHMRSE`) %>%  #LFC + SE
  dplyr::filter(`typeHMR` != 0) %>%
  dplyr::arrange(desc(`typeHMR`)) %>%
  dplyr::mutate(direct = ifelse(`typeHMR` > 0, "Positive LFC", "Negative LFC"))

df_fig_type$taxon_id <- factor(df_fig_type$taxon_id, levels = df_fig_type$taxon_id)
df_fig_type$direct <- factor(df_fig_type$direct, levels = c("Positive LFC", "Negative LFC"))

p_type_l <- ggplot(df_fig_type, aes(x = taxon_id, y = `typeHMR`, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `typeHMR` - `typeHMRSE` , ymax = `typeHMR` + `typeHMRSE`),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log Fold Change", title = "HMR vs LMR Microbiome") +
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_type_l

#split top positive and negative LFC
top_pos <- df_fig_type %>% filter(typeHMR > 0) %>% slice_max(typeHMR, n = 5)
top_neg <- df_fig_type %>% filter(typeHMR < 0) %>% slice_min(typeHMR, n = 5)
top_pos_neg <- bind_rows(top_pos, top_neg)

#reorder for plotting
top_pos_neg <- top_pos_neg %>%
  arrange(typeHMR) %>% #sorts negative to positive
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(type_direction = ifelse(typeHMR > 0, "HMR", "LMR")) %>%
  arrange(typeHMR) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id)) %>%
  mutate(type_direction = factor(type_direction, levels = c("LMR", "HMR")))

#bars extending to the right (positive) - taxa enriched in HMR
#bars extending to the left (negative) - taxa enriched in LMR
levels(sample_data(pseq_genus)$type) #LMR is reference, HMR is comparison

colorsW <- c("Flavicella" = "#FFD1DF", #flavo
             "Marinicella" = "#CCCCCC", #other
             "uncultured_15" = "#CCCCCC", #other
             "uncultured_11" = "#CCCCCC", #other
             "Aestuariicella" = "#CCCCCC", #other
             "Pleionea" = "#CCCCCC", #other
             "Fulvivirga" = "#CCCCCC", #other
             "Psychrobium" = "#FFB347", #shew
             "Ekhidna" = "#CCCCCC", #other
             "Bermanella" = "#CCCCCC") #other

## → plot p_type_all_l ----
p_type_all_l <- ggplot(top_pos_neg, aes(x = taxon_id, y = typeHMR, fill = taxon_id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = typeHMR - typeHMRSE, ymax = typeHMR + typeHMRSE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip() 
p_type_all_l
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_type_all_l_plot.jpg", plot = p_type_all_l, scale = 1, width = 3, height = 3.5,
       units = "in", dpi = 300)
ggsave("p_type_all_l_plot.svg", plot = p_type_all_l, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)

# TEMPERATURE ----
# 3. seawater by temp ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "seawater") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

#testing for differential abundance of bacterial families across temp and type 
out <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "temperature + type", #main effects + interaction
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "temperature", #main comparison
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus)$temperature)
table(sample_data(pseq_genus)$type)

#extract results
res <- out$res
res_global <- out$res_global

#log fold changes table
tab_lfc <- res$lfc
colnames(tab_lfc)

col_name = c("Taxon", "Intercept", "temperature", "medium")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#standard errors
tab_se <- res$se
colnames(tab_se) <- col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#p-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#bias-corrected log abundances
samp_frac = out$samp_frac
#replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
#add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
#adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
#show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% datatable(caption = "Bias-corrected log observed abundances")

#visualization 
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

colnames(df_se)[4] <- "temp18SE"

df_fig_temp <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `temperature18`, `temp18SE`) %>% #LFC + SE
  dplyr::filter(`temperature18` != 0) %>%
  dplyr::arrange(desc(`temperature18`)) %>%
  dplyr::mutate(direct = ifelse(`temperature18` > 0, "Positive LFC", "Negative LFC"))

df_fig_temp$taxon_id <- factor(df_fig_temp$taxon_id, levels = df_fig_temp$taxon_id)
df_fig_temp$direct <- factor(df_fig_temp$direct, levels = c("Positive LFC", "Negative LFC"))

p_temp_sw <- ggplot(df_fig_temp, aes(x = taxon_id, y = `temperature18`, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `temperature18` - `temp18SE` , ymax = `temperature18` + `temp18SE`),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log Fold Change", title = "HMR vs LMR Microbiome") +
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_temp_sw

#split top positive and negative LFC
top_pos <- df_fig_temp %>% filter(temperature18 > 0) %>% slice_max(temperature18, n = 5)
top_neg <- df_fig_temp %>% filter(temperature18 < 0) %>% slice_min(temperature18, n = 5)
top_pos_neg <- bind_rows(top_pos, top_neg)

#reorder for plotting 
top_pos_neg <- top_pos_neg %>%
  arrange(temperature18) %>% #sorts negative to positive
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(temp_direction = ifelse(temperature18 > 0, "18˚C", "14˚C")) %>%
  arrange(temperature18) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#bars extending to the right (positive) - taxa enriched in 18
#bars extending to the left (negative) - taxa enriched in 14
levels(sample_data(pseq_genus)$temperature) #14 is reference, 18 is comparison

colorsW <- c("Celeribacter" = "#BEE3BA", #rhodo
             "Alkalimarinus" = "#6495ED", #altero
             "Streptococcus" = "#CCCCCC", #other
             "Neisseria" = "#CCCCCC", #other
             "Sedimentitalea" = "#BEE3BA", #rhodo
             "Owenweeksia" = "#CCCCCC", #other
             "Milano-WF1B-44" = "#CCCCCC", #other
             "Leisingera" = "#BEE3BA", #rhodo
             "Oceaniserpentilla" = "#CCCCCC", #other
             "Nautella" = "#BEE3BA") #rhodo

## → plot p_temp_all_sw ----
p_temp_all_sw <- ggplot(top_pos_neg, aes(x = taxon_id, y = temperature18, fill = taxon_id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = temperature18 - temp18SE, ymax = temperature18 + temp18SE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip() 
p_temp_all_sw
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_temp_all_sw_plot.jpg", plot = p_temp_all_sw, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
ggsave("p_temp_all_sw_plot.svg", plot = p_temp_all_sw, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)

# 4. larval by temp ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "larval") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

#testing for differential abundance of bacterial families across temp and type 
out <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "temperature + type", #main effects + interaction
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "temperature", #main comparison
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05,
  global = TRUE,
  n_cl = 1,
  verbose = TRUE)

table(sample_data(pseq_genus)$temperature)
table(sample_data(pseq_genus)$type)

#extract results
res <- out$res
res_global <- out$res_global

#log fold changes table
tab_lfc <- res$lfc
colnames(tab_lfc)

col_name = c("Taxon", "Intercept", "temperature", "medium")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#standard errors
tab_se <- res$se
colnames(tab_se) <- col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#p-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

#differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#bias-corrected log abundances
samp_frac = out$samp_frac
#replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
#add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
#adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
#show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% datatable(caption = "Bias-corrected log observed abundances")

#visualization
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

colnames(df_se)[4] <- "temp18SE"

df_fig_temp <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `temperature18`, `temp18SE`) %>% #LFC + SE
  dplyr::filter(`temperature18` != 0) %>%
  dplyr::arrange(desc(`temperature18`)) %>%
  dplyr::mutate(direct = ifelse(`temperature18` > 0, "Positive LFC", "Negative LFC"))

df_fig_temp$taxon_id <- factor(df_fig_temp$taxon_id, levels = df_fig_temp$taxon_id)
df_fig_temp$direct <- factor(df_fig_temp$direct, levels = c("Positive LFC", "Negative LFC"))

p_temp_l <- ggplot(df_fig_temp, aes(x = taxon_id, y = `temperature18`, fill = direct, color = direct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `temperature18` - `temp18SE` , ymax = `temperature18` + `temp18SE`),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log Fold Change", title = "HMR vs LMR Microbiome") +
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_temp_l

#split top positive and negative LFC
top_pos <- df_fig_temp %>% filter(temperature18 > 0) %>% slice_max(temperature18, n = 5)
top_neg <- df_fig_temp %>% filter(temperature18 < 0) %>% slice_min(temperature18, n = 5)
top_pos_neg <- bind_rows(top_pos, top_neg)

#reorder for plotting
top_pos_neg <- top_pos_neg %>%
  arrange(temperature18) %>% #sorts negative to positive
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

top_pos_neg <- bind_rows(top_pos, top_neg) %>%
  mutate(temp_direction = ifelse(temperature18 > 0, "18˚C", "14˚C")) %>%
  arrange(temperature18) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id))

#bars extending to the right (positive) - taxa enriched in 18
#bars extending to the left (negative) - taxa enriched in 14
levels(sample_data(pseq_genus)$temperature) #14 is reference, 18 is comparison

colorsW <- c("Neisseria" = "#CCCCCC", #other
             "Haemophilus" = "#CCCCCC", #other
             "Veillonella" = "#CCCCCC", #other
             "Gemella" = "#CCCCCC", #other
             "Porphyromonas" = "#CCCCCC", #other
             "uncultured_61" = "#CCCCCC", #other
             "uncultured_47" = "#CCCCCC", #other
             "uncultured_60" = "#CCCCCC", #other
             "Thalassospira" = "#CCCCCC", #other
             "Nautella" = "#BEE3BA") #rhodo

## → plot p_temp_all_l ----
p_temp_all_l <- ggplot(top_pos_neg, aes(x = taxon_id, y = temperature18, fill = taxon_id)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = temperature18 - temp18SE, ymax = temperature18 + temp18SE),
                width = 0.2, color = "black") +
  scale_fill_manual(values = colorsW) +
  labs(x = "Taxon", y = "Log fold change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black", face = "italic"),
        axis.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  coord_flip()
p_temp_all_l
setwd("~/Desktop/StephV3V4/plots")
ggsave("p_temp_all_l_plot.jpg", plot = p_temp_all_l, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)
ggsave("p_temp_all_l_plot.svg", plot = p_temp_all_l, scale = 1, width = 3, height = 2.5,
       units = "in", dpi = 300)

# 5. global test seawater by treat ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "seawater") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

sample_data(pseq_genus)$treat <- factor(
  sample_data(pseq_genus)$treat,
  levels = c("14LMR", "18LMR", "14HMR", "18HMR"))

colorsW <- c("14HMR" = "#663399", #purple
             "14LMR" = "#2A788EFF", #blue
             "18LMR" = "#D5AB09", #yellow
             "18HMR" = "#7AD151FF") #green
metadata$treat_c <- colorsW[metadata$treat]

out_global <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "treat",
  group = "treat",
  p_adj_method = "BH", #Adjusts p-values using the Benjamini-Hochberg FDR method.
  prv_cut = 0.10, #Removes taxa present in less than 10% of samples.
  lib_cut = 1000, #Removes samples with library size < 1000 reads.
  struc_zero = TRUE, #Detects structural zeros (taxa that are truly absent in a group).
  neg_lb = TRUE, #Ensures negative lower bounds for confidence intervals are allowed.
  alpha = 0.05, #Significance threshold for declaring taxa differentially abundant.
  global = TRUE, #Calculates a global test statistic across all treatment levels, not just pairwise comparisons.
  verbose = TRUE)
out_global

#extract results
res <- out_global$res
res_global <- out_global$res_global

#test statistics
tab_w_global = res_global[, c("taxon", "W")]
tab_w_global %>% datatable(caption = "Test Statistics from the Global Test Result") %>%
  formatRound(c("W"), digits = 2)

#p-values
tab_p_global = res_global[, c("taxon", "p_val")]
tab_p_global %>% datatable(caption = "P-values from the Global Test Result") %>%
  formatRound(c("p_val"), digits = 2)

#adjusted p-values
tab_q_global = res_global[, c("taxon", "q_val")]
tab_q_global %>% datatable(caption = "Adjusted p-values from the Global Test Result") %>%
  formatRound(c("q_val"), digits = 2)

#differentially abundant taxa
tab_diff_global = res_global[, c("taxon", "diff_abn")]
tab_diff_global %>% datatable(caption = "Differentially Abundant Taxa from the Global Test Result")

#visualization
sig_taxa = res_global %>% dplyr::filter(diff_abn == TRUE) %>% .$taxon

tab_lfc <- res$lfc
colnames(tab_lfc)

df_treat <- tab_lfc %>%
  dplyr::select(taxon, `treat14HMR`, `treat18LMR`, `treat18HMR`) %>%
  dplyr::filter(taxon %in% sig_taxa) %>%
  dplyr::rename(`14HMR vs 14LMR` = treat14HMR,
                `18LMR vs 14LMR` = treat18LMR,
                `18HMR vs 14LMR` = treat18HMR)

df_heat_global <- df_treat %>%
  pivot_longer(cols = -taxon, names_to = "treatment", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat_global$taxon <- factor(df_heat_global$taxon, levels = sort(sig_taxa))

lo  <- floor(min(df_heat_global$value))
up  <- ceiling(max(df_heat_global$value))
mid <- (lo + up) / 2

p_heat_global <- ggplot(df_heat_global, aes(x = treatment, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       na.value = "white", midpoint = mid, limits = c(lo, up), name = NULL) +
  geom_text(aes(treatment, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL,
       title = "Log fold changes for globally significant taxa (relative to 14LMR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_heat_global

#select top 20 taxa by max absolute value
top20_taxa <- df_heat_global %>%
  group_by(taxon) %>%
  summarise(max_abs = max(abs(value), na.rm = TRUE)) %>%
  slice_max(max_abs, n = 20) %>%
  pull(taxon)

df_heat_top20 <- df_heat_global %>% filter(taxon %in% top20_taxa)

#plotting
p_heat_global <- ggplot(df_heat_top20, aes(x = treatment, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    na.value = "white", midpoint = mid, limits = c(lo, up), name = NULL) +
  geom_text(aes(label = value), color = "black", size = 4) +
  labs(
    x = NULL, y = NULL,
    title = "LFC for globally sig taxa in seawater (relative to 14LMR)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1))
p_heat_global

# 6. global test larval by treat ----
pseq_genus <- ps.rarefied %>%
  subset_samples(medium == "larval") %>%
  subset_taxa(!is.na(Genus) & Genus != "")

sample_data(pseq_genus)$treat <- factor(
  sample_data(pseq_genus)$treat,
  levels = c("14LMR", "18LMR", "14HMR", "18HMR"))

out_global <- ancombc(
  data = pseq_genus,
  tax_level = "Genus",
  formula = "treat",
  group = "treat",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE,
  verbose = TRUE)
out_global

#extract results
res <- out_global$res
res_global <- out_global$res_global

#test statistics
tab_w_global = res_global[, c("taxon", "W")]
tab_w_global %>% datatable(caption = "Test Statistics from the Global Test Result") %>%
  formatRound(c("W"), digits = 2)

#p-values
tab_p_global = res_global[, c("taxon", "p_val")]
tab_p_global %>% datatable(caption = "P-values from the Global Test Result") %>%
  formatRound(c("p_val"), digits = 2)

#adjusted p-values
tab_q_global = res_global[, c("taxon", "q_val")]
tab_q_global %>% datatable(caption = "Adjusted p-values from the Global Test Result") %>%
  formatRound(c("q_val"), digits = 2)

#fifferentially abundant taxa
tab_diff_global = res_global[, c("taxon", "diff_abn")]
tab_diff_global %>% datatable(caption = "Differentially Abundant Taxa from the Global Test Result")

#visualization 
sig_taxa = res_global %>% dplyr::filter(diff_abn == TRUE) %>% .$taxon

tab_lfc <- res$lfc
colnames(tab_lfc)

df_treat <- tab_lfc %>%
  dplyr::select(taxon, `treat14HMR`, `treat18LMR`, `treat18HMR`) %>%
  dplyr::filter(taxon %in% sig_taxa) %>%
  dplyr::rename(`14HMR vs 14LMR` = treat14HMR,
                `18LMR vs 14LMR` = treat18LMR,
                `18HMR vs 14LMR` = treat18HMR)

df_heat_global <- df_treat %>%
  pivot_longer(cols = -taxon, names_to = "treatment", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat_global$taxon <- factor(df_heat_global$taxon, levels = sort(sig_taxa))

lo  <- floor(min(df_heat_global$value))
up  <- ceiling(max(df_heat_global$value))
mid <- (lo + up) / 2

p_heat_global <- ggplot(df_heat_global, aes(x = treatment, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       na.value = "white", midpoint = mid, limits = c(lo, up), name = NULL) +
  geom_text(aes(treatment, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL,
       title = "Log fold changes for globally significant taxa (relative to 14LMR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_heat_global

#select top 20 taxa by max absolute value
top20_taxa <- df_heat_global %>%
  group_by(taxon) %>%
  summarise(max_abs = max(abs(value), na.rm = TRUE)) %>%
  slice_max(max_abs, n = 20) %>%
  pull(taxon)

df_heat_top20 <- df_heat_global %>%
  filter(taxon %in% top20_taxa)

#plotting
p_heat_global <- ggplot(df_heat_top20, aes(x = treatment, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    na.value = "white", midpoint = mid, limits = c(lo, up), name = NULL) +
  geom_text(aes(label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL,
       title = "LFC for globally sig taxa in larvae (relative to 14LMR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_heat_global

############## fin
