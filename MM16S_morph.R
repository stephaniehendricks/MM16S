#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for mamma mia project Stephanie Hendricks et al

#clear working environment
rm(list=ls())

#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github("chankinonn/GroupStruct", force = TRUE)

#download packages
library(dplyr)
library(stringr)
library(rstatix)
library(ggplot2) 
library(tidyr)
library(GroupStruct)
library(Rmisc)
library(data.table)
library(lme4)

#set working directory
setwd("~/Desktop/StephV3V4/morphology")
# 1. length analysis - on clean data (outliers removed) ----
data <- read.csv("MM16S_length_clean.csv") #6dpf data only

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
shapiro.test(data$Length_mean) #p = 0.5609 NORMAL across timepoints

data$temp <- as.factor(data$temp)
#Levene's variance: p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(data, Length_mean ~ temp) #p = 995 YES
levene_test(data, Length_mean ~ water) #p = 0.169 YES

#length data meets assumptions of normality and homogeneity

#6dpf / 0hpi length analysis
#meets criteria for ANOVA
length_aov <- aov(data = data, Length_mean ~ water * temp)
summary(length_aov)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#water         1   2320    2320   4.310 0.04009 * 
#temp          1   4955    4955   9.205 0.00297 **
#water:temp    1   2818    2818   5.234 0.02395 *
res_length_aov <- TukeyHSD(length_aov)
res_length_aov
# sig diff: 
# LMR-HMR p = 0.0400861 SIG
# 18-14 p = 0.0029755 SIG
#                     diff       lwr       upr     p adj
#LMR:14-HMR:14   0.7588957 -14.60924 16.127029 0.9992343 Not
#HMR:18-HMR:14  -3.3818447 -18.74998 11.986289 0.9397951 Not
#LMR:18-HMR:14 -21.9371503 -37.44144 -6.432862 0.0019321 SIG
#HMR:18-LMR:14  -4.1407405 -19.75478 11.473299 0.9002830 Not
#LMR:18-LMR:14 -22.6960460 -38.44411 -6.947978 0.0015237 SIG
#LMR:18-HMR:18 -18.5553056 -34.30337 -2.807238 0.0139394 SIG

# 3. plotting ---- 
data$treat <- paste(data$temp, data$water, sep = "")

#add ˚C to temp and treat
data$tempC <- paste0(data$temp, "\u00B0C")  #\u00B0 = degree symbol
data$treatC <- paste(data$tempC, data$water, sep = " ")

orderList = c("14LMR", "18LMR", "14HMR", "18HMR")
names(orderList) = c("14ºC LMR", "18ºC LMR", "14ºC HMR", "18ºC HMR")

#define colors for each treat group
treatColors <- c("14LMR" = "#2A788EFF", "14HMR" = "#663399", 
                 "18LMR" = "#D5AB09", "18HMR" = "#7AD151FF")
treatColors2 <- c("14LMR" = "#014D4E", "14HMR" = "#341539",
                  "18LMR" = "#BA8E23", "18HMR" = "#228C22") 

lengthBoxplot <- ggplot(data, aes(x = treat, y = Length_mean, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.7)) +
  geom_point(aes(color = treat),
             size = 3, alpha = 0.5,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  scale_fill_manual(values = treatColors) + 
  scale_color_manual(values = treatColors2) + 
  scale_x_discrete(limits = orderList, labels = names(orderList)) +
  theme_classic() +
  labs(x = NULL,
       y = "Mean body length (µm)") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.title = element_text(size = 11, colour = "black"),
        legend.position = "none") 
print(lengthBoxplot)
setwd("~/Desktop/StephV3V4/plots")
ggsave("lengthBoxplot.jpg", plot = lengthBoxplot, scale = 1, 
       width = 4, height = 5, units = "in", dpi = 300)

waterColors <- c("LMR" = "#FFB347", "HMR" = "#FFD1DF")
waterColors2 <- c("LMR" = "#CC5500", "HMR" = "#E75480")

lengthBoxplot_type <- ggplot(data, aes(x = water, y = Length_mean, fill = water)) +
  geom_boxplot(width = 0.6) +
  geom_point(aes(color = water), size = 3, alpha = 0.5, position = position_jitter(width = 0.15)) +
  scale_fill_manual(values = waterColors) +
  scale_color_manual(values = waterColors2) +
  scale_x_discrete(limits = c("LMR", "HMR")) +
  theme_classic() +
  labs(x = NULL, y = "Mean body length (µm)") +
  theme(axis.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none")
print(lengthBoxplot_type)
setwd("~/Desktop/StephV3V4/plots")
ggsave("lengthBoxplot_type.jpg", plot = lengthBoxplot_type, scale = 1, 
       width = 2, height = 4.5, units = "in", dpi = 300)
ggsave("lengthBoxplot_type.svg", plot = lengthBoxplot_type, scale = 1, 
       width = 2.5, height = 4.5, units = "in", dpi = 300)

tempColors <- c("14°C" = "#BEE3BA", "18°C" = "#6495ED")
tempColors2 <- c("14°C" = "#95BB72", "18°C" = "#00008B")

lengthBoxplot_temp <- ggplot(data, aes(x = tempC, y = Length_mean, fill = tempC)) +
  geom_boxplot(width = 0.6) +
  geom_point(aes(color = tempC), size = 3, alpha = 0.5, position = position_jitter(width = 0.15)) +
  scale_fill_manual(values = tempColors) +
  scale_color_manual(values = tempColors2) +
  scale_x_discrete(limits = c("14°C", "18°C")) +
  theme_classic() +
  labs(x = NULL, y = "Mean body length (µm)") +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "none")
print(lengthBoxplot_temp)
setwd("~/Desktop/StephV3V4/plots")
ggsave("lengthBoxplot_temp.jpg", plot = lengthBoxplot_temp, scale = 1, 
       width = 2, height = 4.5, units = "in", dpi = 300)
ggsave("lengthBoxplot_temp.svg", plot = lengthBoxplot_temp, scale = 1, 
       width = 2.5, height = 4.5, units = "in", dpi = 300)

############## fin
