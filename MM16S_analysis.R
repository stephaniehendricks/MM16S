#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al
#qiime2 https://github.com/jbisanz/qiime2R/tree/master
#https://micca.readthedocs.io/en/latest/phyloseq.html
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html

#clear working environment
rm(list=ls())

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R", force = TRUE)

#download packages
#data
library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
#comm eco
library(vegan) 
library(phyloseq)
library(qiime2R)
#plot aesthetics 
library(ggplot2)
library(ggvenn)
#stats 
library(car) 
library(btools) #faith's div
library(rstatix) 
library(lme4)
library(lmerTest)  
library(DHARMa)
library(sjPlot)
library(emmeans)
#set memory limit
library(usethis) 
usethis::edit_r_environ()

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

#remove chloroplast,  mitochondria, and eukaryota
clean <- subset_taxa(physeq, 
                        !Order %in% "Chloroplast" &
                        !Family %in% "Mitochondria" &
                        !grepl("Eukaryota", Kingdom) &
                        !is.na(Phylum) & Phylum != "")
clean #otu_table()   OTU Table:         [ 2832 taxa and 67 samples ]
physeq <- clean

#check cleaning worked
#get number of families 
length(unique(na.omit(tax_table(physeq)[, "Family"]))) 
#256 but 255 KNOWN - one NA (omitted above) and Unknown_Family (in 256 count)
#all families listed 
unique(na.omit(as.character(tax_table(physeq)[, "Family"])))

#get number of genera 
length(unique(na.omit(tax_table(physeq)[, "Genus"]))) 
#449 

#should return FALSE if removal worked
any(grepl("Chloroplast", tax_table(physeq)[, "Order"], ignore.case = TRUE), na.rm = TRUE)
any(grepl("Mitochondria", tax_table(physeq)[, "Family"], ignore.case = TRUE), na.rm = TRUE)
any(grepl("Eukaryota", tax_table(physeq)[, "Kingdom"], ignore.case = TRUE), na.rm = TRUE)

#prepare OTU matrix to plot rarefaction curves
#extract OTU table
otu_data <- otu_table(physeq)
#check dimensions
dim(otu_data) #2832  67
#convert to matrix
otu_matrix <- as(otu_data, "matrix")
#transpose matrix
otu_matrix_transposed <- t(otu_matrix)
#remove rows with 0 total counts
otu_matrix_transposed <- otu_matrix_transposed[rowSums(otu_matrix_transposed) > 0, ]
#remove any rows with NA 
otu_matrix_transposed <- na.omit(otu_matrix_transposed)

#plot rarefaction curves
#rarecurve(otu_matrix_transposed, step = 50, cex = 0.5)
#confirm value for min sample.size below  
min(sample_sums(physeq)) #[1] 82246
#rarefy samples without replacement to simulate even number of reads per sample
  #rarefaction depth chosen is the minimum sample depth in the dataset
ps.rarefied <- rarefy_even_depth(physeq, 
                                 rngseed = 1, 
                                 sample.size = 82246, 
                                 replace = FALSE)
#15OTUs were removed because they are no longer present in any sample after random subsampling

#number before rarefying
ntaxa(physeq) #2832
#number after rarefying
ntaxa(ps.rarefied) #2817

#get number of families 
length(unique(na.omit(tax_table(ps.rarefied)[, "Family"]))) 
#256 but 255 KNOWN - one NA (omitted above) and Unknown_Family (in 256 count)
#all families listed 
unique(na.omit(as.character(tax_table(ps.rarefied)[, "Family"])))

#removing NA family
not_na_family <- !is.na(tax_table(ps.rarefied)[, "Family"])
keep_taxa <- taxa_names(ps.rarefied)[not_na_family]
noNA <- prune_taxa(keep_taxa, ps.rarefied)
ntaxa(noNA) #2728

#replace noNA to be ps.rarefied
ps.rarefied <- noNA 
#check that NA is gone (don't use na.omit this time)
unique((as.character(tax_table(ps.rarefied)[, "Family"])))

# 3. calculate relative abundances ----
#agglomerate at Family level
fams <- tax_glom(ps.rarefied, taxrank = "Family", NArm = TRUE)

#transform to relative abundance
RA <- transform_sample_counts(fams, function(x) x / sum(x))

#get abundance matrix
abund_matrix <- otu_table(RA)
if (!taxa_are_rows(RA)) abund_matrix <- t(abund_matrix)

#extract taxonomy
tax <- as.data.frame(tax_table(fams))
tax$Family <- as.character(tax$Family)

#maximum relative abundance per family
family_max_ra <- apply(abund_matrix, 1, max)

#families reaching >=25% RA in at least one sample
top_family <- tax$Family[family_max_ra >= 0.25]

#relabel low abundance families as Other
tax$Family <- ifelse(
  tax$Family %in% top_family &
  tax$Family != "Unknown_Family",
  tax$Family, "Other")

#replace higher taxonomy ranks with Other
tax_ranks <- colnames(tax)
family_col_index <- which(tax_ranks == "Family")
for (i in seq_len(family_col_index - 1)) {
  tax[tax$Family == "Other", i] <- "Other"
}

#replace taxonomy
tax_table(fams) <- tax_table(as.matrix(tax))

#re-agglomerate
fams_filt <- tax_glom(fams, taxrank = "Family", NArm = TRUE)

table(tax_table(fams_filt)[, "Family"]) #only shows 10 families to plot

# plot giant bar plot for figs and supp ----
#get ready for plotting
colors <- c(
  Alteromonadaceae        = "#6495ED", #periwinkle
  Psychromonadaceae       = "#00CED1", #teal
  Crocinitomicaceae       = "#A569BD", #dark purple
  Flavobacteriaceae       = "#FFD1DF", #light pink
  Pseudoalteromonadaceae  = "#FFD700", #yellow
  Rhodobacteraceae        = "#BEE3BA", #light green
  Saccharospirillaceae    = "#F08080", #coral
  Shewanellaceae          = "#FFB347", #orange
  Vibrionaceae            = "#CDC0F8", #light purple
  Other                   = "#CCCCCC" #grey
)

#make Other last in legend 
#extract plot data from phyloseq object
p <- plot_bar(fams_filt, fill = "Family")
plot_data <- p$data
family_levels <- c(sort(unique(plot_data$Family[plot_data$Family != "Other"])), "Other")
plot_data$Family <- factor(plot_data$Family, levels = family_levels)

#add word "day" after numbers
plot_data$day <- paste0(plot_data$day, " dpf")
#make seawater samples first then larval
plot_data$medium <- factor(plot_data$medium, levels = c("seawater", "larval"))

#plot barplot
quartz()
family_barplot <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  facet_wrap(~ medium + treatC + day, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 11, colour = "black"),
    legend.text = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 11, colour = "black"),
    plot.title = element_text(size = 11, colour = "black"),
    legend.position = "right") +
  labs(x = NULL, y = "Relative abundance") +
  scale_y_continuous(expand = c(0, 0))
family_barplot
setwd("~/Desktop/StephV3V4/plots")
ggsave("family_plot.jpg", plot = family_barplot, width = 24, height = 4, units = "in", dpi = 300)
#ggsave("family_barplot.svg", plot = family_barplot, device = "svg", width = 24, height = 4)     

# SUMMARIES ----
# 1. seawater summary stats ----
RA_sw <- subset_samples(RA, medium == "seawater")
table(sample_data(RA_sw)$medium) #34 samples
RA_sw <- psmelt(RA_sw)

family_of_interest <- "Flavobacteriaceae" #change this to family of interest

summary_sw <- RA_sw %>%
  dplyr::select(Sample, treat, day, type, temperature, Family, Abundance) %>%
  tidyr::complete(Sample, Family = family_of_interest, fill = list(Abundance = 0)) %>%
  tidyr::fill(treat, day, type, temperature, .direction = "downup") %>%
  dplyr::filter(Family == family_of_interest) %>%
  dplyr::group_by(treat, day, type, temperature) %>%
  dplyr::summarise(min_RA  = min(Abundance),
    max_RA  = max(Abundance),
    mean_RA = mean(Abundance),
    .groups = "drop") %>%
  dplyr::mutate(across(c(min_RA, max_RA, mean_RA), ~ .x * 100))
summary_sw

overall_mean_RA_sw <- mean(summary_sw$mean_RA)
overall_mean_RA_sw #average treatment-level abundance

#overall_mean_RA_sw <- mean(RA_sw$Abundance[RA_sw$Family == family_of_interest]) * 100
#overall mean abundance across all samples

overall_mean_RA_filt_sw <- summary_sw %>%
  filter(day == "6") %>%
  #filter(type == "LMR") %>%
  summarise(overall_mean_RA_sw = mean(mean_RA, na.rm = TRUE)) %>%
  pull(overall_mean_RA_sw)
overall_mean_RA_filt_sw

# 2. larvae summary stats ---- 
RA_l <- subset_samples(RA, medium == "larval")
table(sample_data(RA_l)$medium) #33 samples
RA_l <- psmelt(RA_l)

family_of_interest <- "Shewanellaceae" #change this to family of interest

summary_l <- RA_l %>%
  dplyr::select(Sample, treat, day, type, temperature, Family, Abundance) %>%
  tidyr::complete(Sample, Family = family_of_interest, fill = list(Abundance = 0)) %>%
  tidyr::fill(treat, day, type, temperature, .direction = "downup") %>%
  dplyr::filter(Family == family_of_interest) %>%
  dplyr::group_by(treat, day, type, temperature) %>%
  dplyr::summarise(min_RA  = min(Abundance),
    max_RA  = max(Abundance),
    mean_RA = mean(Abundance),
    .groups = "drop") %>%
  dplyr::mutate(across(c(min_RA, max_RA, mean_RA), ~ .x * 100))
summary_l

overall_mean_RA_l <- mean(summary_l$mean_RA)
overall_mean_RA_l #average treatment-level abundance

#overall_mean_RA_l <- mean(RA_l$Abundance[RA_l$Family == family_of_interest]) * 100
#overall mean abundance across all samples

overall_mean_RA_filt_l <- summary_l %>%
  filter(day == "6") %>%
  filter(treat == "18LMR") %>%
  summarise(overall_mean_RA_l = mean(mean_RA)) %>%
  pull(overall_mean_RA_l) 
overall_mean_RA_filt_l

# plot families over time ----
#define colors for each treat 
treatColors <- c("14HMR" = "#663399", #purple
                 "14LMR" = "#2A788EFF", #blue
                 "18LMR" = "#D5AB09", #yellow
                 "18HMR" = "#7AD151FF") #green

#family of interest - numbers coincide with y axis limits
#Vibrionaceae 60,000
#Alteromonadaceae 40,000
#Rhodobacteraceae 30,000
#Saccharospirillaceae 80,000
#Flavobacteriaceae
#Pseudoalteromonadaceae
target_family <- "Saccharospirillaceae" #change this to family of interest

#counts - not relative abundance 
family_counts <- tax_glom(ps.rarefied, taxrank = "Family")
family_counts <- psmelt(family_counts)

#larval samples only
family_counts_filt <- family_counts %>%
  dplyr::filter(medium == "larval") %>%
  dplyr::select(Sample, day, treat, Family, Abundance)

family_counts_filt$day <- factor(family_counts_filt$day, levels = c("2", "4", "6"))

#summary stats
family_counts_filt <- family_counts_filt %>%
  tidyr::complete(Sample, Family = target_family,
                  fill = list(Abundance = 0)) %>%
  tidyr::fill(day, treat, .direction = "downup") %>%
  dplyr::filter(Family == target_family)

family_means <- family_counts_filt %>%
  dplyr::group_by(day, treat) %>%
  dplyr::summarise(mean_abundance = mean(Abundance),
                   se = sd(Abundance) / sqrt(n()),
                   .groups = "drop")

#plot family over time 
quartz()
family_counts_plot <- ggplot(family_counts_filt, aes(x = day, y = Abundance, color = treat)) +
  geom_point(size = 3, alpha = 0.5, position = position_jitter(width = 0.15)) +
  geom_line(data = family_means, aes(y = mean_abundance, group = treat), linewidth = 1) +
  geom_errorbar(data = family_means, aes(y = mean_abundance,
                    ymin = mean_abundance - se,
                    ymax = mean_abundance + se),
                width = 0.1, linewidth = 0.5) +
  labs(title = target_family,
       x = "Days post-fertilization",
       y = "Abundance") +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black", face = "italic"),
        legend.position = "none") +
  scale_color_manual(values = treatColors) +
  scale_y_continuous(limits = c(0, max(family_counts_filt$Abundance) * 1.2),
                     labels = function(x) paste0(x / 1000, "k"),
                     expand = c(0, 0))
family_counts_plot
#ggsave(filename = "sacc.svg", plot = family_counts_plot, width = 2.75, height = 4.25, units = "in", dpi = 300)

# plot venns ----
#aggregate to family
ps_family <- tax_glom(ps.rarefied, taxrank = "Family", NArm = TRUE)

#convert to long format
df <- psmelt(ps_family) %>% dplyr::select(Sample, day, medium, type, Family, Abundance)

#presence and absence
df_pa <- df %>% mutate(present = Abundance > 0)

#families per group
fam_presence <- df_pa %>%
  filter(present) %>%
  group_by(day, medium, type) %>%
  summarise(families = list(unique(Family)), .groups = "drop")
print(fam_presence)

#plot venns LMR and HMR at 6 dpf
plot_venn <- function(day_val, type_val) {
  sub <- fam_presence %>%
    filter(day == day_val, type == type_val)
  venn_list <- list(
    seawater = unique(unlist(sub$families[sub$medium == "seawater"])),
    larval   = unique(unlist(sub$families[sub$medium == "larval"])))
  ggvenn(venn_list,
    fill_color = c("skyblue", "#D8BFD8"),
    fill_alpha = 0.5,
    stroke_size = 1,
    text_size = 8,
    set_name_size = 8) +
    ggtitle(paste0(type_val, " (", day_val, " dpf)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24, margin = margin(b = 15)))
}

quartz()
plot_venn(6, "LMR")
ggsave("venn_LMR.jpg", plot = plot_venn(6, "LMR"))

# ALPHA DIVERSITY ----
rich <- estimate_richness(ps.rarefied)
fpd <- estimate_pd(ps.rarefied)

rich$Faith_PD <- fpd$PD
rich$SampleID <- rownames(rich)

metadata$SampleID <- rownames(metadata)

rich_combined <- rich %>%
  mutate(Faith_PD = fpd$PD,
         SampleID = rownames(rich)) %>%
  left_join(metadata, by = "SampleID")
colnames(rich_combined)

#define colors for points
treatColors2 <- c("14LMR" = "#014D4E",
                  "14HMR" = "#341539",
                  "18LMR" = "#D5AB09",
                  "18HMR" = "#228C22")

rich_combined$treatC <- factor(rich_combined$treatC, levels = c("14°C LMR", "18°C LMR", "14°C HMR", "18°C HMR"))

# 1. seawater summary and plot ----
seawater <- rich_combined %>% filter(medium == "seawater")

rich_summary_sw <- seawater %>%
  group_by(day, treat) %>%
  rstatix::get_summary_stats(Faith_PD, type = "mean_sd")
rich_summary_sw

#plot boxplot
quartz()
alphadiv_sw_plot <- ggplot(seawater, aes(x = day, y = Faith_PD, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +
  geom_point( aes(color = treat), size = 3, alpha = 0.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  theme_classic() +
  labs(x = "Days post-fertilization", y = "Faith's phylogenetic diversity") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values = treatColors) +
  scale_color_manual(values = treatColors2) +
  coord_cartesian(ylim = c(min(seawater$Faith_PD) - 1,
                           max(seawater$Faith_PD) + 2)) +
  facet_grid(. ~ treatC)
alphadiv_sw_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("FPD_alphadiv_day_water_boxplot.jpg", plot = alphadiv_sw_plot, scale = 1, width = 4, height = 4.5, units = "in", dpi = 300)
#ggsave("FPD_alphadiv_day_water_boxplot.svg", plot = alphadiv_sw_plot, scale = 1, width = 4, height = 4, units = "in", dpi = 300)

#alpha diversity testing for normality and homogeneity of variances
# → seawater normality ----
rich_sw <- subset(rich_combined, medium == "seawater")
rich_sw$day <- factor(rich_sw$day, levels = c("2", "4", "6"))

#SW
shapiro.test(rich_sw$Faith_PD) #W = 0.92828, p-value = 0.02788
by(rich_sw$Faith_PD, rich_sw$type, shapiro.test) #LMR is sig
#LMR W = 0.69424, p-value = 0.0001459
#HMR W = 0.91862, p-value = 0.1222
by(rich_sw$Faith_PD, rich_sw$temperature, shapiro.test) 
#14 W = 0.92085, p-value = 0.1526
#18 W = 0.91954, p-value = 0.1451
by(rich_sw$Faith_PD, rich_sw$treat, shapiro.test) #14LMR is sig
#14LMR W = 0.71386, p-value = 0.003229
#18LMR W = 0.9264, p-value = 0.4839
#14HMR W = 0.87863, p-value = 0.1519
#18HMR W = 0.91269, p-value = 0.3352
by(rich_sw$Faith_PD, rich_sw$day, shapiro.test) 
#2 W = 0.89193, p-value = 0.147
#4 W = 0.8718, p-value = 0.06888
#6 W = 0.89241, p-value = 0.149
by(rich_sw$Faith_PD, rich_sw$chamber, shapiro.test) #harry is sig
#bill W = 0.87826, p-value = 0.09889
#harry W = 0.85289, p-value = 0.03986
#tom W = 0.94768, p-value = 0.6145

#Levene
leveneTest(Faith_PD ~ type, data = rich_sw) #group  1  2.2752 0.1413
leveneTest(Faith_PD ~ temperature, data = rich_sw) #group  1  0.7525 0.3922
leveneTest(Faith_PD ~ treat, data = rich_sw) #group  3   0.901 0.4522
leveneTest(Faith_PD ~ day, data = rich_sw) #group  2  3.8928 0.03102 *
leveneTest(Faith_PD ~ chamber, data = rich_sw) #group  2  1.8953 0.1673

# → seawater stats ----
set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

rich_sw$temperature <- factor(rich_sw$temperature)
rich_sw$type <- factor(rich_sw$type)
rich_sw$day <- factor(rich_sw$day)

lmm_sw <- lmer(Faith_PD ~ temperature * type * day + (1|chamber), data = rich_sw)
anova(lmm_sw)
#                     Sum Sq Mean Sq NumDF   DenDF F value       Pr(>F)    
#Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#temperature           32.74   32.74     1 20.225  2.4881   0.13023    
#type                 837.15  837.15     1 19.996 63.6210 1.222e-07 ***
#  day                   53.00   26.50     2 20.158  2.0139   0.15945    
#temperature:type       0.86    0.86     1 20.225  0.0655   0.80052    
#temperature:day       13.74    6.87     2 20.054  0.5221   0.60113    
#type:day             151.56   75.78     2 20.158  5.7591   0.01051 *  
#  temperature:type:day  34.01   17.00     2 20.054  1.2922   0.29658 

#testing with no random effect (chamber) to see if genotype is sig
lm_sw <- lm(Faith_PD ~ temperature * type * day, data = rich_sw)
anova(lm_sw)
#                     Df Sum Sq Mean Sq F value    Pr(>F)    
#temperature           1  29.20   29.20  2.1643   0.15541    
#type                  1 895.20  895.20 66.3541 4.368e-08 ***
#  day                   2  63.15   31.58  2.3405   0.11979    
#temperature:type      1   1.59    1.59  0.1180   0.73448    
#temperature:day       2   7.50    3.75  0.2779   0.75996    
#type:day              2 148.04   74.02  5.4865   0.01167 *  
#  temperature:type:day  2  35.14   17.57  1.3023   0.29206 

res <- residuals(lmm_sw)
fit <- fitted(lmm_sw)
#check linearity and homogeneity (variance)
plot(fit, res)
abline(h = 0, lty = 2)

#check normality of residuals
qqnorm(res)
qqline(res)
hist(res, breaks = 20)
shapiro.test(res) #W = 0.97414, p-value = 0.5842

#check residuals vs predictors
plot(rich_sw$temperature, res)
plot(rich_sw$type, res)
plot(rich_sw$day, res)
plot(lmm_sw)

#make sure no sig problems detected
sim_res <- simulateResiduals(lmm_sw)
plot(sim_res)

testUniformity(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

#post hoc
#does Faith's PD differ between LMR and HMR, averaging across temperature?
emm_type_day <- emmeans(lmm_sw, ~ type | day, lmer.df = "satterthwaite")
pairs(emm_type_day, adjust = "tukey")
#HMR has greater Faith's PD than LMR at every sampling day
#difference is largest at 4 dpf
#consistent with the significant type × day interaction

#day = 2:
#  contrast  estimate   SE   df t.ratio p.value
#LMR - HMR    -6.45 2.22 20.2  -2.902  0.0087

#day = 4:
#  contrast  estimate   SE   df t.ratio p.value
#LMR - HMR   -15.92 2.09 19.9  -7.603  <.0001

#day = 6:
#  contrast  estimate   SE   df t.ratio p.value
#LMR - HMR    -7.75 2.22 20.2  -3.485  0.0023

#within each microbial type, does Faith’s PD change over time (from 2 → 4 → 6 dpf)?
emm_day_type <- emmeans(lmm_sw, ~ day | type, lmer.df = "satterthwaite")
pairs(emm_day_type, adjust = "tukey")
#Faith's PD did not significantly change across days in LMR
#Faith's PD was significantly higher at day 4 than day 2 in HMR
#Faith's PD was significantly higher at day 4 than day 6 in HMR

#type = LMR:
#  contrast    estimate   SE   df t.ratio p.value
#day2 - day4    2.722 2.22 20.2   1.224  0.4532
#day2 - day6    1.622 2.35 20.7   0.690  0.7716
#day4 - day6   -1.100 2.22 20.2  -0.495  0.8746

#type = HMR:
#  contrast    estimate   SE   df t.ratio p.value
#day2 - day4   -6.747 2.09 19.9  -3.222  0.0115
#day2 - day6    0.326 2.09 19.9   0.156  0.9868
#day4 - day6    7.073 2.09 19.9   3.377  0.0081

#how much variance does random effect explain (chamber/genotype)
VarCorr(lmm_sw)
vc <- as.data.frame(VarCorr(lmm_sw))
chamber_var <- vc$vcov[vc$grp == "chamber"]
residual_var <- vc$vcov[vc$grp == "Residual"]
ICC <- chamber_var / (chamber_var + residual_var)
ICC #[1] 0.02485713

#using a package to get ICC 
tab_model(lmm_sw) #0.02

# 2. larval summary and plot ----
larval <- rich_combined %>% filter(medium == "larval")

rich_summary_l <- larval %>%
  group_by(day, treat) %>%
  rstatix::get_summary_stats(Faith_PD, type = "mean_sd")
rich_summary_l

#plot boxplot
alphadiv_l_plot <- ggplot(larval, aes(x = day, y = Faith_PD, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +
  geom_point(aes(color = treat), size = 3, alpha = 0.5,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  theme_classic() +
  labs(x = "Days post-fertilization",
       y = "Faith's phylogenetic diversity") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values = treatColors) +
  scale_color_manual(values = treatColors2) +
  coord_cartesian(ylim = c(min(larval$Faith_PD) - 1, max(larval$Faith_PD) + 2)) +
  facet_grid(. ~ treatC)
alphadiv_l_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("FPD_alphadiv_day_larval_boxplot.jpg", plot = alphadiv_l_plot, scale = 1, width = 4, height = 4.5, units = "in", dpi = 300)
#ggsave("FPD_alphadiv_day_larval_boxplot.svg", plot = alphadiv_l_plot, scale = 1, width = 4, height = 4, units = "in", dpi = 300)

# → larval normality ----
rich_l <- subset(rich_combined, medium == "larval")
rich_l$day <- factor(rich_l$day, levels = c("2", "4", "6"))

#SW
shapiro.test(rich_l$Faith_PD) #W = 0.89077, p-value = 0.003095
by(rich_l$Faith_PD, rich_l$type, shapiro.test) 
#LMR W = 0.9294, p-value = 0.2673
#HMR W = 0.94658, p-value = 0.3739
by(rich_l$Faith_PD, rich_l$temperature, shapiro.test) #18 is sig
#14 W = 0.90471, p-value = 0.09566
#18 W = 0.85367, p-value = 0.01221
by(rich_l$Faith_PD, rich_l$treat, shapiro.test) 
#14LMR W = 0.93364, p-value = 0.5823
#18LMR W = 0.96454, p-value = 0.8521
#14HMR W = 0.94391, p-value = 0.6236
#18HMR W = 0.8437, p-value = 0.06348
by(rich_l$Faith_PD, rich_l$day, shapiro.test) #day 6 is sig
#2 W = 0.86873, p-value = 0.09662
#4 W = 0.90492, p-value = 0.1836
#6 W = 0.83681, p-value = 0.02868
by(rich_l$Faith_PD, rich_l$chamber, shapiro.test) #bill and tom are sig
#bill W = 0.84915, p-value = 0.04161
#harry W = 0.90176, p-value = 0.1942
#tom W = 0.82299, p-value = 0.01888

#Levene
leveneTest(Faith_PD ~ type, data = rich_l) #group  1  2.7028 0.1103
leveneTest(Faith_PD ~ temperature, data = rich_l) #group  1  0.0516 0.8218
leveneTest(Faith_PD ~ treat, data = rich_l) #group  3  1.3692 0.2718
leveneTest(Faith_PD ~ day, data = rich_l) #group  2  0.2544  0.777
leveneTest(Faith_PD ~ chamber, data = rich_l) #group  2  0.4278 0.6559

# → larval stats ----
set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

rich_l$temperature <- factor(rich_l$temperature)
rich_l$type <- factor(rich_l$type)
rich_l$day <- factor(rich_l$day)

lmm_l <- lmer(Faith_PD ~ temperature * type * day + (1|chamber), data = rich_l)
anova(lmm_l)
#                      Sum Sq Mean Sq NumDF   DenDF  F value          Pr(>F)    
#Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
#temperature             0.21    0.21     1 19.162   0.0184    0.8936    
#type                 2586.98 2586.98     1 19.000 224.6458 5.581e-12 ***
#  day                    22.25   11.13     2 19.115   0.9662    0.3983    
#temperature:type        2.82    2.82     1 19.162   0.2450    0.6263    
#temperature:day        13.29    6.65     2 19.154   0.5772    0.5709    
#type:day               41.54   20.77     2 19.115   1.8037    0.1916    
#temperature:type:day   11.39    5.69     2 19.154   0.4943    0.6175  

#testing with no random effect (chamber) to see if genotype is sig
lm_l <- lm(Faith_PD ~ temperature * type * day, data = rich_l)
anova(lm_l)
#                     Df Sum Sq Mean Sq F value    Pr(>F)    
#temperature           1    8.38    8.38   0.6005    0.4470    
#type                  1 2616.71 2616.71 187.5042 6.155e-12 ***
#  day                   2   28.30   14.15   1.0141    0.3798    
#temperature:type      1    0.97    0.97   0.0697    0.7944    
#temperature:day       2    7.70    3.85   0.2757    0.7617    
#type:day              2   38.71   19.35   1.3868    0.2718    
#temperature:type:day  2   10.23    5.11   0.3664    0.6975  

res <- residuals(lmm_l)
fit <- fitted(lmm_l)
#check linearity and homogeneity (variance)
plot(fit, res)
abline(h = 0, lty = 2)

#check normality of residuals
qqnorm(res)
qqline(res)
hist(res, breaks = 20)
shapiro.test(res) #W = 0.94849, p-value = 0.1202

#check residuals vs predictors
plot(rich_l$temperature, res)
plot(rich_l$type, res)
plot(rich_l$day, res)
plot(lmm_l)

#make sure no sig problems detected
sim_res <- simulateResiduals(lmm_l)
plot(sim_res)

testUniformity(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

#how much variance does random effect explain (chamber/genotype)
VarCorr(lmm_l)
vc <- as.data.frame(VarCorr(lmm_l))
chamber_var <- vc$vcov[vc$grp == "chamber"]
residual_var <- vc$vcov[vc$grp == "Residual"]
ICC <- chamber_var / (chamber_var + residual_var)
ICC #[1] 0.1748165

#using a package to get ICC 
tab_model(lmm_l) #0.17

# BETA DIVERSITY ----
set.seed(123)

sample_data(ps.rarefied)$day <- factor(sample_data(ps.rarefied)$day, levels = c("2", "4", "6"))
sample_data(ps.rarefied)$treat <- factor(sample_data(ps.rarefied)$treat, levels = c("14LMR", "14HMR", "18LMR", "18HMR"))

#all samples
bc_dist <- phyloseq::distance(ps.rarefied, method = "bray")
ordination <- ordinate(ps.rarefied, method = "PCoA", distance = bc_dist)

#plot pcoa for all samples
plot_ordination(ps.rarefied, ordination, color = "treat") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps.rarefied)), size = 3, vjust = -0.5) +
  scale_color_manual(values = treatColors) +
  facet_wrap(~ medium) +
  ggtitle("PCoA on Bray-Curtis Distance") +
  theme_classic()

# 1. seawater bray and pcoa ----
beta_sw <- subset_samples(ps.rarefied, medium == "seawater")
meta_sw <- data.frame(sample_data(beta_sw))

#bray and pcoa
bc_dist_sw <- phyloseq::distance(beta_sw, method = "bray")
ordination_sw <- ordinate(beta_sw, method = "PCoA", distance = bc_dist_sw)

#variance explained
var_exp_sw <- ordination_sw$values$Relative_eig
pc1_sw <- round(var_exp_sw[1] * 100, 1)
pc2_sw <- round(var_exp_sw[2] * 100, 1)

#plot pcoa for seawater
betadiv_pcoa_sw <- plot_ordination(beta_sw, ordination_sw, color = "treat", shape = "day")
betadiv_pcoa_sw$layers <- NULL

betadiv_pcoa_sw_plot <- betadiv_pcoa_sw +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(aes(group = treat), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = treatColors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = paste0("PCo1 (", pc1_sw, "%)"),
       y = paste0("PCo2 (", pc2_sw, "%)")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        legend.position = "none")
betadiv_pcoa_sw_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("betadiv_pcoa_water.jpg", plot = betadiv_pcoa_sw_plot, scale = 1, width = 4, height = 4.5, units = "in", dpi = 300)
#ggsave("betadiv_pcoa_water.svg", plot = betadiv_pcoa_sw_plot, scale = 1, width = 4, height = 4.25, units = "in", dpi = 300)

# → seawater stats ---- 
#PERMANOVA
adonis2(bc_dist_sw ~ temperature * type + chamber, data = meta_sw, permutations = 999, by = "terms",)
#temperature       1   0.2332 0.02611  1.2518  0.244    
#type              1   2.5978 0.29082 13.9427  0.001 ***
#  chamber           2   0.5613 0.06284  1.5064  0.119    
#temperature:type  1   0.3232 0.03618  1.7345  0.129     

#test homogeneity of dispersion 
disp_sw <- betadisper(bc_dist_sw, meta_sw$treat)
anova(disp_sw) #Groups     3 0.03030 0.010100  0.5016 0.6841
permutest(disp_sw) #Groups     3 0.03030 0.010100 0.5016    999   0.69

# 2. larval bray and pcoa ----
beta_l <- subset_samples(ps.rarefied, medium == "larval")
meta_l <- data.frame(sample_data(beta_l))

#bray and pcoa
bc_dist_l <- phyloseq::distance(beta_l, method = "bray")
ordination_l <- ordinate(beta_l, method = "PCoA", distance = bc_dist_l)

#variance explained
var_exp_l <- ordination_l$values$Relative_eig
pc1_l <- round(var_exp_l[1] * 100, 1)
pc2_l <- round(var_exp_l[2] * 100, 1)

#plot pcoa for larval
betadiv_pcoa_l <- plot_ordination(beta_l, ordination_l, color = "treat", shape = "day")
betadiv_pcoa_l$layers <- NULL

quartz()
betadiv_pcoa_day_l_plot <- betadiv_pcoa_l +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(aes(group = treat), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = treatColors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_continuous(limits = c(-0.5, 1),
                     breaks = c(-0.5, 0.0, 0.5)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c(-0.5, 0.0, 0.5)) +
  labs(x = paste0("PCo1 (", sprintf("%.1f", 39.0), "%)"),
       y = paste0("PCo2 (", sprintf("%.1f", 16.0), "%)")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        legend.position = "none")
betadiv_pcoa_day_l_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("betadiv_pcoa_day_larval.jpg", plot = betadiv_pcoa_day_l_plot, scale = 1, width = 4, height = 4.5, units = "in", dpi = 300)
#ggsave("betadiv_pcoa_day_larval.svg", plot = betadiv_pcoa_day_l_plot, scale = 1, width = 4, height = 4.25, units = "in", dpi = 300)

# → larval stats ----
#PERMANOVA
adonis2(bc_dist_l ~ temperature * type + chamber, data = meta_l, permutations = 999, by = "terms")
#temperature       1   0.2994 0.03435  1.8835  0.073 .  
#type              1   3.1119 0.35710 19.5778  0.001 ***
#  chamber           2   0.4972 0.05706  1.5641  0.114    
#temperature:type  1   0.5143 0.05902  3.2358  0.006 **  
  
#test homogeneity of dispersion 
disp_l <- betadisper(bc_dist_l, meta_l$treat)
anova(disp_l) #Groups     3 0.00735 0.0024512  0.1918 0.9011
permutest(disp_l) #Groups     3 0.00735 0.0024512 0.1918    999  0.888

# 3. day 6 bray and pcoa ----
beta_day6 <- subset_samples(ps.rarefied, day == "6")
meta_day6 <- data.frame(sample_data(beta_day6))

#bray and pcoa
bc_dist_day6 <- phyloseq::distance(beta_day6, method = "bray")
ordination_day6 <- ordinate(beta_day6, method = "PCoA", distance = bc_dist_day6)

#variance explained
var_exp_day6 <- ordination_day6$values$Relative_eig
pc1_day6 <- round(var_exp_day6[1] * 100, 1)
pc2_day6 <- round(var_exp_day6[2] * 100, 1)

#plot pcoa for all samples at 6 dpf
quartz()
betadiv_pcoa_day6_plot <- plot_ordination(beta_day6, ordination_day6, color = "treat", shape = "medium") +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(aes(group = treat), type = "norm", level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = treatColors) +
  scale_shape_manual(values = c(8, 3)) +
  scale_x_continuous(limits = c(-0.75, 1),
                     breaks = c(-0.5, 0.0, 0.5)) +
  scale_y_continuous(limits = c(-1.25, 1.5),
                     breaks = c(-0.5, 0.0, 0.5)) +
  labs(x = paste0("PCo1 (", pc1_day6, "%)"),
       y = paste0("PCo2 (", pc2_day6, "%)")) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        legend.position = "none")
betadiv_pcoa_day6_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("betadiv_pcoaplot_day6.jpg", plot = betadiv_pcoa_day6_plot, scale = 1, width = 4, height = 3.2, units = "in", dpi = 300)
#ggsave("betadiv_pcoaplot_day6.svg", plot = betadiv_pcoa_day6_plot, scale = 1, width = 4, height = 4.25,  units = "in", dpi = 300)

# → day 6 stats ----
#PERMANOVA
adonis2(bc_dist_day6 ~ temperature * type + chamber, data = meta_day6, permutations = 999, by = "terms",)
#temperature       1   0.3012 0.04772 1.5119  0.148    
#type              1   1.9524 0.30928 9.7988  0.001 ***
#  chamber           2   0.4147 0.06569 1.0407  0.400    
#temperature:type  1   0.4564 0.07230 2.2906  0.034 *  
  
#test homogeneity of dispersion 
disp_day6 <- betadisper(bc_dist_day6, meta_day6$treat)
anova(disp_day6) #Groups     3 0.00449 0.0014968  0.1436 0.9324
permutest(disp_day6) #Groups     3 0.00449 0.0014968 0.1436    999  0.939

# FAPROTAX ----
#code to run FAPROTAX
#./collapse_table.py -i combined.txt -o functional.txt -g FAPROTAX.txt -d "taxonomy" -c "#" -v --force
setwd("~/Desktop/StephV3V4/FAPROTAX_1.2.11") 

functional <- read_tsv("functional.txt") 
head(functional)

#remove samples
samples_to_remove <- c("BS14_6dpf_L", "TS14_2dpf_L")
functional <- functional[, !(colnames(functional) %in% samples_to_remove)]

functional_long <- functional %>%
  filter(rowSums(across(-group)) > 0) %>%
  pivot_longer(cols = -group,
  names_to = "sample",
  values_to = "abundance") %>%
  group_by(group) %>%
  filter(sum(abundance) > 0) %>%
  ungroup()

#update metadata columns
functional_long <- functional_long %>%
  mutate(type = case_when(
      str_detect(sample, "F") ~ "HMR",
      str_detect(sample, "S") ~ "LMR",
      TRUE ~ "Unknown"),
    temperature = case_when(
      str_detect(sample, "14") ~ "14",
      str_detect(sample, "18") ~ "18",
      TRUE ~ "Unknown"),
    medium = case_when(
      str_detect(sample, "L") ~ "larval",
      str_detect(sample, "W") ~ "seawater",
      TRUE ~ "Unknown"),
    day = case_when(
      str_detect(sample, "2dpf") ~ "2",
      str_detect(sample, "4dpf") ~ "4",
      str_detect(sample, "6dpf") ~ "6",
      TRUE ~ "Unknown"),
    treat = paste(temperature, type, sep = ""))
print(functional_long)

#filter functional groups with total abundance >25,000 across all samples
functional_filt <- functional_long %>%
  group_by(group) %>%
  filter(sum(abundance) > 25000) %>%
  ungroup()

#reorder type
functional_filt$type <- factor(functional_filt$type, levels = c("LMR", "HMR"))

#labels for plot
labels <- c("aerobic_chemoheterotrophy",
             "chemoheterotrophy",
             "fermentation",
             "nitrate_reduction")

plot_df <- functional_filt %>%
  filter(group %in% labels) %>%
  mutate(group = factor(group, levels = rev(sort(unique(group)))))

plot_df$medium <- factor(plot_df$medium, levels = c("seawater", "larval"))

#plot heatmap
functional_plot <- ggplot(plot_df, aes(x = sample, y = group, fill = abundance)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "#F3E5F5", high = "#4A148C") +  
  labs(x = "",
       y = "Functional group") +
  facet_grid(~ medium + type, scales = "free_x", space = "free_x", switch = "y") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(),
        strip.text.x = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11, colour = "black"),
        panel.spacing.x = unit(0.3, "lines"))
functional_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("functional_plot.jpg", plot = functional_plot, scale = 1, width = 8, height = 3, units = "in", dpi = 300)
#ggsave("functional_plot.svg", plot = functional_plot, scale = 1, width = 8, height = 2.5, units = "in", dpi = 300)

# ADDITIONAL PLOTS (legends and supp) ---- 
#for treat 
svg("~/Desktop/StephV3V4/plots/legend_treat_square.svg", width = 10, height = 7)                
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1) #blank plot
plotSym <- c(15, 15, 15, 15) #square

legend("right", 
       legend = c("14\u00b0C LMR", 
                  "18\u00b0C LMR",
                  "14\u00b0C HMR",
                  "18\u00b0C HMR"),
       pch = plotSym, 
       pt.cex = 5.1, #point size
       col = c("#2A788EFF", "#D5AB09", "#663399ff", "#7AD151FF"), #point color
       cex = 2.5, 
       bty = 'n', 
       title = "Treatment",        
       title.adj = 0.1, 
       title.cex = 3)
dev.off()

#for dpf
svg("~/Desktop/StephV3V4/plots/legend_dpf.svg", width = 10, height = 7)  
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1) #blank plot
plotSym <- c(16, 17, 15) #circle, triangle, square
legend("right", 
       legend = c("2", 
                  "4",
                  "6"),
       pch = plotSym, pt.cex = 5.1, #size of points
       #pt.bg = c("black", "black", "black"), #fill color
       col = "black", #border color
       cex = 2.5, 
       bty = 'n', 
       title = "dpf",        
       title.adj = 0.1, #move title
       title.cex = 3)
dev.off()

#for larval vs seawater
svg("~/Desktop/StephV3V4/plots/legend_medium.svg", width = 10, height = 7)                
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1) #blank plot
plotSym <- c(3, 8) #plus and star
legend("right", 
       legend = c("seawater", 
                  "larval"),
       pch = plotSym, pt.cex = 5.1, #size of points
       #pt.bg = c("black", "black", "black"), #fill color
       col = "black", #border color
       cex = 2.5, 
       bty = 'n', 
       #title = "sample",        
       title.adj = 0.1, #move title
       title.cex = 3)
dev.off()

#boxplots of alpha diversity - for supp fig
rarefiedPlot <- plot_richness(ps.rarefied, x = "type", measures = c("Observed", "Shannon", "Chao1", "Simpson")) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Water condition",
       y = "Alpha diversity measure")
plot(rarefiedPlot)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rarefiedPlot, filename = "rarefiedPlotSupp-type.jpg", scale = 1, width = 7.2, height = 6, units = c("in"), dpi = 300)

rareTemp <- plot_richness(ps.rarefied, x = "temperature", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Temperature (˚C)",
       y = "Alpha diversity measure")
plot(rareTemp)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rareTemp, filename = "rarefiedPlotSupp-temp.jpg", scale = 1, width = 7.2, height = 6, units = c("in"), dpi = 300)

rareDay <- plot_richness(ps.rarefied, x = "day", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Days post-fertilization",
       y = "Alpha diversity measure")
plot(rareDay)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rareDay, filename = "rarefiedPlotSupp-day.jpg", scale = 1, width = 7.2, height = 6.5, units = c("in"), dpi = 300)

rareMedium <- plot_richness(ps.rarefied, x = "medium", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Sample type",
       y = "Alpha diversity measure")
plot(rareMedium)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rareMedium, filename = "rarefiedPlotSupp-medium.jpg", scale = 1, width = 7.2, height = 6.5, units = c("in"), dpi = 300)

rareChamber <- plot_richness(ps.rarefied, x = "chamber", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Chamber",
       y = "Alpha diversity measure")
plot(rareChamber)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rareChamber, filename = "rarefiedPlotSupp-chamber.jpg", scale = 1, width = 7.2, height = 6.5, units = c("in"), dpi = 300)

plot_richness(ps.rarefied, x = "medium", measures = c("Shannon")) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Sample type",
       y = "Alpha diversity measure") +
  facet_wrap(~ day)  

rareDay <- plot_richness(ps.rarefied, x = "treat", measures = c("Shannon")) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Treatment",
       y = "Alpha diversity measure") +
  facet_wrap(~ day + medium)  
plot(rareDay)

#reorder treatments
ps.rarefied@sam_data$treat <- factor(ps.rarefied@sam_data$treat, 
                                     levels = c("14LMR", "18LMR", "14HMR", "18HMR"))
rareTreat <- plot_richness(ps.rarefied, x = "treat", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Treatment",
       y = "Alpha diversity measure")
plot(rareTreat)
setwd("~/Desktop/StephV3V4/plots")
ggsave(plot = rareTreat, filename = "rarefiedPlotSupp-treat.jpg", scale = 1, width = 7.2, height = 6.5, units = c("in"), dpi = 300)

#extract BEFORE grouping into Other - for supp figure
p_all <- plot_bar(RA, fill = "Family")
plot_data_all <- p_all$data

plot_data_all$day <- paste0(plot_data_all$day, " dpf")
plot_data_all$medium <- factor(plot_data_all$medium,
                               levels = c("seawater", "larval"))

#order families alphabetically
plot_data_all$Family <- factor(plot_data_all$Family,
                               levels = sort(unique(plot_data_all$Family)))

quartz()
family_barplot_all <- ggplot(plot_data_all, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  facet_wrap(~ medium + treatC + day, scales = "free_x", nrow = 1) +
  scale_fill_viridis_d(guide = guide_legend(nrow = 60)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 11, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11, colour = "black"),
        legend.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 5, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none",
        legend.key.size = unit(0.1, "cm")) +
  labs(x = NULL, y = "Relative abundance") +
  scale_y_continuous(expand = c(0, 0))
family_barplot_all
#setwd("~/Desktop/StephV3V4/plots")
#ggsave("family_barplot_all.jpg", plot = family_barplot_all, width = 24, height = 4, units = "in", dpi = 300)

#larval and seawater plot
alphadiv_combined_facet_plot <- ggplot(rich_combined, aes(x = day, y = Faith_PD, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +
  geom_point(aes(color = treat), size = 3, alpha = 0.5,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  theme_classic() +
  labs(x = "Days post-fertilization", y = "Faith's phylogenetic diversity") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values = treatColors) +
  scale_color_manual(values = treatColors2) +
  coord_cartesian(ylim = c(min(rich_combined$Faith_PD) - 1,
                           max(rich_combined$Faith_PD) + 2)) +
  facet_grid(medium ~ treatC)
alphadiv_combined_facet_plot
setwd("~/Desktop/StephV3V4/plots")
ggsave("alphadiv_combined_facet_plot.jpg", plot = alphadiv_combined_facet_plot, scale = 1, width = 5, height = 3.2, units = "in", dpi = 300)

############## fin
