#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al

## This script is to calculate the bacterial load from the qPCR of Amy's 
## samples using the standard curve from the controls.

# 1. Load required packages ----
library(tidyverse) # for data wrangling
library(ggplot2) # for plotting
library(ggrepel) # for labels on plots
library(car) # for levene test
library(ARTool) # for non parametric data with interactions 

## Optional: turn off scientific notation for this session.
options(scipen = 999)

# 2. Load qPCR data ----
## Import the CSV file I made from the excel output from the run. I removed all
## of the unnecessary columns for this analysis, so it just has "Well", "Sample",
## and "Cq".
setwd("~/Desktop/StephV3V4")
dat <- read_csv("Taqman_pan_bacterial_20251031_140007_795BR05383_STRADER_LAB -  Quantification Cq Results.csv")

# Keep only certain columns
df <- dat %>% select(Well, Sample, Cq)

# 3. Standard Curve ----
standard_conc <- c(S1 = 75, S2 = 7.5, S3 = 0.75, S4 = 0.075, S5 = 0.0075, S6 = 0.00075)

## Mutate the data frame to contain the known ng for the standards
df <- df %>%
  mutate(is_standard = Sample %in% names(standard_conc),
    quantity_ng = if_else(is_standard, standard_conc[Sample], NA_real_))

## → Fit standard curve line ----
## Filter data frame to only have the standards, and calculate the log10 of the
## known ng values. The log10 is used because qPCR is exponential, and the Cq is
## linearly related to the log of the starting DNA amount. 
std_curve <- df %>%
  filter(is_standard, !is.na(Cq)) %>%
  mutate(log_quantity = log10(quantity_ng))

## Fit a linear model where Cq is the dependent variable and log10(ng) is the 
## independent variable. 
model <- lm(Cq ~ log_quantity, data = std_curve)

## Look at the summary of the model (it looks really good R2 > 99%)
summary(model)
#Residual standard error: 0.6707 on 16 degrees of freedom
#Multiple R-squared:  0.9934,	Adjusted R-squared:  0.993 
#F-statistic:  2400 on 1 and 16 DF,  p-value: < 0.00000000000000022

## → Visualize the standard curve ----

## Extract the line equation and r2 from the model info
coefs <- coef(model)
intercept <- round(coefs[1], 3)
slope <- round(coefs[2], 3)
r2 <- round(summary(model)$r.squared, 3)

## Create label for plot
eqn_text <- paste0("Cq = ", intercept, " + ", slope, " × log10(ng)\nR² = ", r2)

## Create the plot
ggplot(std_curve, aes(x = log_quantity, y = Cq)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(aes(label = Sample),
                  size = 4,
                  box.padding = 0.4,       # Add more space around labels
                  point.padding = 0.5,     # Add more space around points
                  max.overlaps = Inf,      # Show all labels
                  segment.color = "gray50",
                  segment.size = 0.4,
                  force = 2,               # Stronger repulsion
                  force_pull = 0.5,        # Gentle pull back to origin
                  nudge_x = 0.3,           
                  nudge_y = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Standard Curve: 16S qPCR",
       x = expression(log[10]*"(ng of 16S DNA)"),
       y = "Cq Value") +
  annotate("text",
           x = max(std_curve$log_quantity),  # Far right
           y = max(std_curve$Cq),            # Top
           label = eqn_text,
           hjust = 1,                        # Right-aligned
           size = 4,
           fontface = "italic",
           color = "red") +
  theme_classic(base_size = 14)

# 4. Calculate sample concentrations ----
## Mutate the data frame to have new columns with the calculated concentration 
## (ng). Note that the log_predicted_quantity is the log10 value of the predicted
## DNA amount based on the Cq, while the predicted_quantity is the actual amount
## of DNA present in the well based on the linear model equation. 
df <- df %>%
  mutate(log_predicted_quantity = (Cq - coef(model)[1]) / coef(model)[2],
    predicted_quantity = 10^log_predicted_quantity)

# 5. Filter data by experiment ----

# Define sample sets
mamma_samples <- c("MMWS18", "MMWF18", "MMWS14", "MMWF14")

## → Mamma Mia ----
mamma_mia <- df %>%
  filter(Sample %in% mamma_samples, !is.na(predicted_quantity)) %>%
  mutate(temp = str_extract(Sample, "14|18"),
    type = case_when(
      Sample == "MMWS14" ~ "LMR",
      Sample == "MMWS18" ~ "LMR",
      Sample == "MMWF14" ~ "HMR",
      Sample == "MMWF18" ~ "HMR",
      TRUE ~ NA_character_),
    predicted_quantity = round(predicted_quantity, 3))

mamma_mia$dilution_factor <- "2"
mamma_mia$dilution_factor[mamma_mia$type == "HMR"] <- "4"
mamma_mia$dilution_factor <- as.numeric(mamma_mia$dilution_factor)
mamma_mia <- mamma_mia %>% 
  mutate(quantTotal_ng = predicted_quantity*dilution_factor)

## Stats
#create a new treat column in mamma_mia by combining column temperature and type in metadata
mamma_mia$treat <- paste(mamma_mia$temp, mamma_mia$type, sep = "")

shapiro.test(mamma_mia$quantTotal_ng) #W = 0.71161, p-value = 0.001093
by(mamma_mia$quantTotal_ng, mamma_mia$type, shapiro.test) #both pass
leveneTest(quantTotal_ng ~ type, data = mamma_mia) #group  1  5.9711 0.03464 *

by(mamma_mia$quantTotal_ng, mamma_mia$temp, shapiro.test) #both fail
leveneTest(quantTotal_ng ~ temp, data = mamma_mia) #group  1  6.4543 0.02934 *

mamma_mia$type <- factor(mamma_mia$type)
mamma_mia$temp <- factor(mamma_mia$temp)

#using ANOVA ART for non parametric data with interactions 
model <- art(quantTotal_ng ~ type * temp, data = mamma_mia)
summary(model)
anova(model)
#            Df Df.res F value   Pr(>F)   
#1 type       1      8 24.9231 0.001063 **
#  2 temp       1      8  8.9109 0.017464  *
#  3 type:temp  1      8  4.0563 0.078784  .

art.con(model, "type:temp")
# contrast        estimate   SE df t.ratio p.value
#HMR,14 - HMR,18    3.000 1.33  8   2.250  0.1895
#HMR,14 - LMR,14    7.333 1.33  8   5.500  0.0026
#HMR,14 - LMR,18    7.667 1.33  8   5.750  0.0019
#HMR,18 - LMR,14    4.333 1.33  8   3.250  0.0468
#HMR,18 - LMR,18    4.667 1.33  8   3.500  0.0330
#LMR,14 - LMR,18    0.333 1.33  8   0.250  0.9941

qqnorm(residuals(model))
qqline(residuals(model))

# 6. Average technical replicates ----

## → Mamma Mia ----
mamma_summary <- mamma_mia %>%
  group_by(Sample, temp, type) %>%
  summarise(mean_ng = mean(quantTotal_ng, na.rm = TRUE),
    sd_ng = sd(quantTotal_ng, na.rm = TRUE),
    .groups = "drop")
mamma_summary

# 7. Plot bacterial load ----

# Now add Treatment using case_when (safer and more flexible than recode)
mamma_summary <- mamma_summary %>%
  mutate(treat = case_when(
    Sample == "MMWF14" ~ "14˚C HMR",
    Sample == "MMWF18" ~ "18˚C HMR",
    Sample == "MMWS14" ~ "14˚C LMR",
    Sample == "MMWS18" ~ "18˚C LMR",
    TRUE ~ NA_character_))

mamma_summary$treat <- factor(mamma_summary$treat,
  levels = c("14˚C LMR", "18˚C LMR", "14˚C HMR", "18˚C HMR"))

treatColors <- c("14˚C LMR" = "#2A788EFF", "14˚C HMR" = "#663399", 
                 "18˚C LMR" = "#D5AB09", "18˚C HMR" = "#7AD151FF")

# Generate the plot
load <- ggplot(mamma_summary, aes(x = treat, y = mean_ng, fill = treat)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  geom_errorbar(aes(ymin = mean_ng - sd_ng, ymax = mean_ng + sd_ng),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = NULL,
       y = "Mean absolute bacterial load (ng)") +
  scale_fill_manual(values = treatColors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # <--- Here
  theme_classic() +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.ticks = element_blank(),  
    legend.text = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 11, colour = "black"),
    plot.title = element_text(size = 11, colour = "black"),
    legend.position = "none")
load
setwd("~/Desktop/StephV3V4/plots")
ggsave("load_plot.jpg", plot = load, scale = 1, width = 2, height = 4.5,
       units = "in", dpi = 300)
ggsave("load_plot.svg", plot = load, 
       scale = 1, width = 1.5, height = 3.5, units = "in", dpi = 300)

############## fin
