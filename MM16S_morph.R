#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al

#clear working environment
rm(list=ls())

#download packages
library(dplyr)
library(rstatix)
library(ggplot2) 

#set working directory
setwd("~/Desktop/StephV3V4/morphology")
# 1. length analysis - on clean data (outliers removed) ----
data <- read.csv("MM16S_length_clean.csv") #6 dpf data only

data <- data %>% rename(type = water)
data <- data %>% rename(temperature = temp)

sum_length <- data %>% 
  group_by(treat) %>% 
  get_summary_stats(Length_mean, type = "mean_sd") 
sum_length
#  treat  variable        n  mean    sd
#1 HMR_14 Length_mean    32  231.  21.3
#2 HMR_18 Length_mean    30  228.  21.4
#3 LMR_14 Length_mean    30  232.  25.0
#4 LMR_18 Length_mean    29  209.  25.0

# 2. stats ----
#Shapiro-Wilk normality: p>0.05 is NOT sig diff from normal
shapiro.test(data$Length_mean) #p = 0.5609 NORMAL 

data$temperature <- as.factor(data$temperature)
#Levene's variance: p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(data, Length_mean ~ temperature) #p = 0.995 YES
levene_test(data, Length_mean ~ type) #p = 0.169 YES

#length data meets assumptions of normality and homogeneity

#6dpf / 0hpi length analysis
#meets criteria for ANOVA
length_aov <- aov(data = data, Length_mean ~ temperature * type)
summary(length_aov)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#temperature        1   5007    5007   9.301 0.00283 **
#  type               1   2268    2268   4.214 0.04233 * 
#  temperature:type   1   2818    2818   5.234 0.02395 * 

# 3. plotting ---- 
data$treat <- paste(data$temperature, data$type, sep = "")
data$tempC <- paste0(data$temperature, "\u00B0C")  #\u00B0 = degree symbol

#define colors for boxes
treatColors <- c(
  "14LMR" = "#2A788EFF", #blue
  "14HMR" = "#663399", #purple
  "18LMR" = "#D5AB09", #yellow
  "18HMR" = "#7AD151FF" #green
)

#define colors for points
treatColors2 <- c("14LMR" = "#014D4E",
                  "14HMR" = "#341539",
                  "18LMR" = "#D5AB09",
                  "18HMR" = "#228C22")

data$treat <- paste0(gsub("°C", "", data$tempC), data$type)
data$treat <- factor(data$treat, levels = c("14LMR", "14HMR", "18LMR", "18HMR"))
data$temp_plot <- factor(data$tempC, levels = c("14°C", "18°C"))

lengthBoxplot_temp_water <- ggplot(data, aes(x = temp_plot, y = Length_mean, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +
  geom_point(aes(color = treat), size = 3, alpha = 0.5, position = position_jitterdodge(
      jitter.width = 0.1, dodge.width = 0.7)) +
  scale_fill_manual(values = treatColors) +
  scale_color_manual(values = treatColors2) +
  theme_classic() +
  labs(x = NULL, y = "Mean body length (µm)") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") 
print(lengthBoxplot_temp_water)
setwd("~/Desktop/StephV3V4/plots")
ggsave("lengthBoxplot_temp_water.svg", plot = lengthBoxplot_temp_water, scale = 1, 
       width = 5.3, height = 3.5, units = "in", dpi = 300)

############## fin
