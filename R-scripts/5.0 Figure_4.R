rm(list=ls())

#load packages
library(ggplot2)
library(Rmisc)
library(cowplot)
library(visreg)
library(lmerTest)
library(ggpubr)
library(ggbeeswarm)
library(ggpattern)

#load data
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
data <- read.csv(file= "data/datafish1-114.csv", header = TRUE)

#organise data structure
str(data)
data$ox <- as.factor(data$ox)
data$ploidy <- as.factor(data$ploidy)
data$sex <- as.factor(data$sex)
data$batch <- as.factor(data$batch)
head(data)

# weeks to days
data$age_days <- data$age * 7 
data

# create and reorder treatments
data$treat <- as.factor(paste(data$r_temp, data$a_temp, data$ploidy, data$ox))
levels(data$treat)
data$treat <- factor(data$treat,levels = c("23 23 2n hyp", "23 23 2n norm", "23 28 2n hyp", "23 28 2n norm",
                                           "23 23 3n hyp", "23 23 3n norm", "23 28 3n hyp", "23 28 3n norm",
                                           "28 28 2n hyp", "28 28 2n norm", "28 23 2n hyp", "28 23 2n norm",
                                           "28 28 3n hyp", "28 28 3n norm", "28 23 3n hyp", "28 23 3n norm")) # change the order of the boxplot
#-------------------------------------------------------------------------------
#create figure 4
p3 <- ggplot() + 
  geom_rect_pattern(aes(xmin=0,xmax=2.5,ymin=-Inf,ymax=Inf), pattern_fill= '#49adf5',
      pattern = 'gradient',
      pattern_fill2 = '#49adf5', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=2.5,xmax=4.5,ymin=-Inf,ymax=Inf), pattern_fill= '#49adf5',
      pattern = 'gradient',
      pattern_fill2 = '#eb4634', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=4.5,xmax=6.5,ymin=-Inf,ymax=Inf), pattern_fill= '#49adf5',
      pattern = 'gradient',
      pattern_fill2 = '#49adf5', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=6.5,xmax=8.5,ymin=-Inf,ymax=Inf), pattern_fill= '#49adf5',
      pattern = 'gradient',
      pattern_fill2 = '#eb4634', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=8.5,xmax=10.5,ymin=-Inf,ymax=Inf), pattern_fill= '#eb4634',
      pattern = 'gradient',
      pattern_fill2 = '#eb4634', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=10.5,xmax=12.5,ymin=-Inf,ymax=Inf), pattern_fill= '#eb4634',
      pattern = 'gradient',
      pattern_fill2 = '#49adf5', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=12.5,xmax=14.5,ymin=-Inf,ymax=Inf), pattern_fill= '#eb4634',
      pattern = 'gradient',
      pattern_fill2 = '#eb4634', alpha = 0.5) +
  geom_rect_pattern(aes(xmin=14.5,xmax=Inf,ymin=-Inf,ymax=Inf), pattern_fill= '#eb4634',
      pattern = 'gradient',
      pattern_fill2 = '#49adf5', alpha = 0.5) +
  scale_x_discrete() +
  geom_boxplot(data = data, aes(x = treat, y = length, fill = ox), width = 0.6) +
  scale_fill_manual(values=c("grey99", "grey65"), labels = c("Hypoxia", "Normoxia")) + #change the labels of the fill variable
  geom_beeswarm(data = data, aes(x = treat, y = length, color = ploidy), size = 0.5) +
  scale_color_manual(values=c("#295E11", "#f4dc01")) +
  geom_vline(xintercept=c(4.5, 12.5), linetype = 2) +
  geom_vline(xintercept=8.5, linetype = 1) +
  xlab("Treatment (RT, AT, ploidy, oxygen)") +
  ylab("Length (mm)") +
  guides(color = "none") +
  labs(fill="Oxygen") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p3 

#create labels for figure 4
p4<-p3 + 
  annotate("text", x=2.5, y=29, label = "2n", fontface = "bold", size = 4.4) +
  geom_label(aes(x=4.5, y=29.3, label = "RT 23°C"), fontface = "bold", size = 4.5, fill = "white") +
  annotate("text", x=6.5, y=29, label = "3n", fontface = "bold", size = 4.4) +
  annotate("text", x=10.5, y=29, label = "2n", fontface = "bold", size = 4.4) +
  geom_label(aes(x=12.5, y=29.3, label = "RT 28°C"), fontface = "bold", size = 4.5, fill = "white") +
  annotate("text", x=14.5, y=29, label = "3n", fontface = "bold", size = 4.4) +
  annotate("point", x = 1, y = 23.91, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 2, y = 23.91, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 3, y = 23.91, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 4, y = 23.91, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 5, y = 21.57, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 6, y = 25.21, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 7, y = 21.57, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 8, y = 21.57, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 9, y = 18.71, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 10, y = 18.71, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 11, y = 18.71, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 12, y = 18.71, shape = 23, fill = "#295E11", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 13, y = 20.31, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 14, y = 20.31, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 15, y = 20.31, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5) +
  annotate("point", x = 16, y = 20.31, shape = 23, fill = "#f4dc01", color = "#5d5d5d", alpha = 0.7, size = 3.5)

p4  
if(T){
  png("figures/Figure 4.png",width = 10,height = 5,units = "in",res = 600)
 
  p4
  
  dev.off()  
}

################################################################################
# Calculate average length of fish of a specific rearing T and ploidy level at the start of acclimation #
#load data
getwd()
dat <- read.csv(file= "data/Groeidata R.csv", header = TRUE, sep = ";", dec = "." )

#Check and organize data structure
str(dat)
names(dat)[1] <- "tank"
names(dat)[2] <- "condition"
dat$condition<-as.factor(dat$condition)
names(dat)[3] <- "temp"
names(dat)[4] <- "ploidy"
names(dat)[5] <- "fdate"
names(dat)[6] <- "fish"
names(dat)[7] <- "age" #this is age in days
names(dat)[8] <- "tl"
dat$cond <- paste(dat$ploidy, dat$temp)
dat$cond <- as.factor(dat$cond)

head(dat)

# subset data for calcualtions
# Lenght at maturity 
data2n28 <- subset(dat, dat$age=="98"&dat$cond=="2n 28")
data3n28 <- subset(dat, dat$age=="98"&dat$cond=="3n 28")
data2n23 <- subset(dat, dat$age=="126"&dat$cond=="2n 23")
data3n23 <- subset(dat, dat$age=="126"&dat$cond=="3n 23")
data3n23.2 <- subset(dat, dat$age=="119"&dat$cond=="3n 23"&dat$fdate=="28-1-2021") #B6 reared as extra batch, but not measured at 126 days
mean(data2n28$tl)
log10(mean(data2n28$tl))
mean(data3n28$tl, na.rm = TRUE)
log10(mean(data3n28$tl, na.rm = TRUE))
mean(data2n23$tl, na.rm = TRUE)
log10(mean(data2n23$tl, na.rm = TRUE))
mean(data3n23$tl, na.rm = TRUE)
log10(mean(data3n23$tl, na.rm = TRUE))
mean(data3n23.2$tl, na.rm = TRUE) #length at maturity differs slightly. Length at start acclimation differs more (larger difference in age)
log10(mean(data3n23.2$tl, na.rm = TRUE))

#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/5.0_figure_4_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################

