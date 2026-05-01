############Final length and weight analysis############
rm(list=ls())

#load packages
library(ggplot2)
library(Rmisc)
library(cowplot)
library(visreg)
library(lmerTest)
library(ggpubr)
library(ggimage)
library(rsvg)

#Load dataset 1
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
data <- read.csv(file= "data/datafish1-114.csv", header = TRUE)

#check and reorganize data
str(data)
data$ox <- as.factor(data$ox)
data$ploidy <- as.factor(data$ploidy)
data$sex <- as.factor(data$sex)
data$batch <- as.factor(data$batch)

head(data)

# weeks to days
data$age_days <- data$age * 7 
data$treat <- as.factor(paste(data$r_temp, data$a_temp, data$ploidy, data$ox))
levels(data$treat)

# Plotting weight for visual check #
plot(weight~length, data=data)
plot(weight~length, col=data$sex, data=data)# pch=c(data$ploidy),
legend(28,100,unique(data$sex),col=1:length(data$sex),pch=1)

# Calculate average length of fish of a specific rearing T and ploidy level at the start of acclimation #
#load dataset 2
getwd()
dat <- read.csv(file= "data/Groeidata R.csv", header = TRUE, sep = ";", dec = "." )

#check data and reorganize
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

#calculate length
# Lenght at maturity #
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

# Length at start acclimation #
mean(data2n28$tl) #length at maturity same as length at start acclimation for batch 1
data3n28sa <- subset(dat, dat$age=="105"&dat$cond=="3n 28") # acclimation started one week later for batch 3
mean(data3n28sa$tl, na.rm = TRUE)
data2n23sa <- subset(dat, dat$age=="133"&dat$cond=="2n 23") # acclimation started one week later for batch 4
mean(data2n23sa$tl)
data3n23sa <- subset(dat, dat$age=="133"&dat$cond=="3n 23") #acclimation started one week later for batch 2
mean(data3n23sa$tl, na.rm = TRUE)
mean(data3n23.2$tl, na.rm = TRUE) #acclimation started at 17 weeks instead of 18 for bacth 6. Fish were classified as adults this week

data$l1 <- 0 # average length at start acclimation
data$l1[which(data$batch=="b1")] <- paste(mean(data2n28$tl))
data$l1[which(data$batch=="b3")] <- paste(mean(data3n28sa$tl, na.rm = TRUE))
data$l1[which(data$batch=="b4")] <- paste(mean(data2n23sa$tl))
data$l1[which(data$batch=="b2")] <- paste(mean(data3n23sa$tl, na.rm = TRUE))
data$l1[which(data$batch=="b6")] <- paste(mean(data3n23.2$tl, na.rm = TRUE))

#organize data structure
str(data)
data$ploidy <-as.factor(data$ploidy)
data$ox <-as.factor(data$ox)
data$sex <-as.factor(data$sex)
data$batch <-as.factor(data$batch)
data$l1 <- as.numeric(data$l1)

# split male and female and log transform length and weight#
data$ploisex <- paste(data$ploidy, data$sex)
data$ploisex[6:8] <- "NA"
data$ploisex <- as.factor(data$ploisex)

newdata <- data[-c(6:8),]

newdata$log_length <- paste(log10(newdata$length))
newdata$log_l1 <- paste(log10(newdata$l1))
newdata$log_length <- as.numeric(newdata$log_length)
newdata$log_l1 <- as.numeric(newdata$log_l1)

newdata$log_weight <- paste(log10(newdata$weight))
newdata$log_length <- paste(log10(newdata$length))
newdata$log_length <- as.numeric(newdata$log_length)
newdata$log_weight <- as.numeric(newdata$log_weight)

newdata$a_temp <- as.numeric(newdata$a_temp)

#drop level ploisex, because ploisex only has 3 levels
levels(droplevels(newdata$ploisex))

# Best model for length from 1.0 model selection length
m13 <- lmer(log_length~ploisex+ox+a_temp+r_temp+log_l1+(1|batch), na.action = na.exclude, data=newdata) # - all 2-way interactions
summary(m13)
anova(m13)
AIC(m13)

# best model for Weight from 2.0 model selection weight

m3 <- lmer(log10(weight)~ploisex*(log_length+a_temp+ox)+r_temp+(1|batch), na.action = na.exclude, data=newdata) 

################################################################################
######################         Create figures 5, 6 and 7   #####################
################################################################################
p1a <- visreg(m13, "a_temp",by = "ploisex",overlay = TRUE,band = FALSE,gg = TRUE,points = list(position = position_jitter(width = .18))) +
  aes(linetype = ploisex,shape = ploisex) + 
  xlab("Acclimation temperature (°C)") +
  ylab(expression(Log[10]~"final length (mm)")) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw()
  theme(legend.position="none")
#annotate("text", x=-Inf, y=Inf, label=expression(paste("*", italic(p), " = 0.024")) , vjust=1.5, hjust=-0.1)
#female=circle, male=triangle, 3n=diamand
p1a

p4a <- visreg(m3,"a_temp",by="ploisex", overlay=TRUE, band = FALSE, gg=TRUE,points = list(position = position_jitter(width = .18))) + 
  aes(linetype = ploisex,shape = ploisex) + 
  xlab("Acclimation temperature (°C)") +
  ylab(expression(Log[10]~"mass (mg)")) + #ylim(4.5,6.5) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw() +
  theme(legend.position="none")
  #annotate("text", x=-Inf, y=Inf, label=expression(paste("** int. ploidy/sex, ", italic(p), " = 0.01")) , vjust=1.5, hjust=-0.1) 

p4a

p2a <- visreg(m13,"ox",by="ploisex", overlay=TRUE, band = FALSE, gg=TRUE) +
  aes(linetype = ploisex,shape = ploisex) + 
  xlab("Oxygen condition") +
  ylab(expression(Log[10]~"final length (mm)")) + #ylim(2.8,3.4) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw() +
  theme(legend.position="none")
  #annotate("text", x=-Inf, y=Inf, label=expression(paste("**", italic(p), " = 0.001")) , vjust=1.5, hjust=-0.1)

p2a

p5a <- visreg(m3,"ox",by="ploisex", overlay=TRUE, band = FALSE, gg=TRUE) + 
  aes(linetype = ploisex,shape = ploisex) +
  xlab("Oxygen condition") +
  ylab(expression(Log[10]~"mass (mg)")) + #ylim(4.5,6.5) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw() +
  theme(legend.position="none")
  #annotate("text", x=-Inf, y=Inf, label=expression(paste("***", italic(p), " < 0.001")) , vjust=1.5, hjust=-0.1)

p5a

p3a <- visreg(m13,"log_l1",by="ploisex", overlay=TRUE, jitter = TRUE, band = FALSE, gg=TRUE) +
  aes(linetype = ploisex,shape = ploisex) +
  geom_abline(slope = 1, linetype = 2, color = "gray") +
  #geom_image(df_img, aes(x=0, y = length), image = image) + #doesn't work to add an icon of a fish
  xlab(expression(Log[10]~"average initial length (mm)")) +
  ylab(expression(Log[10]~"final length (mm)")) + #ylim(2.8,3.4) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")))
  #annotate("text", x=-Inf, y=Inf, label=expression(paste("***", italic(p), " < 0.001")) , vjust=1.5, hjust=-0.1)

p3a

p6a <- visreg(m3,"log_length",by="ploisex", overlay=TRUE, jitter = TRUE, band = FALSE, gg=TRUE) +
  aes(linetype = ploisex,shape = ploisex) +
  geom_abline(slope = 1, linetype = 2, intercept = 4.5, color = "gray") +
  xlab(expression(Log[10]~"final length (mm)")) + 
  ylab(expression(Log[10]~"mass (mg)")) + #ylim(4.5,6.5) +
  scale_fill_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11", "2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_colour_manual(name = "Ploidy level/sex",values = c("2n f" = "#295E11","2n m" = "#295E11","3n m" = "#f4dc01")) +
  scale_linetype_manual(name = "Ploidy level/sex",values = c("2n f" = "dashed","2n m" = "solid","3n m" = "solid")) +
  scale_shape_manual(name = "Ploidy level/sex",values = c("2n f" = 16,"2n m" = 17,"3n m" = 18)) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")))
  #annotate("text", x=-Inf, y=Inf, label=expression(paste("***", italic(p), " < 0.001")) , vjust=1.5, hjust=-0.1)

p6a
#-------------------------------------------------------------------------------
#Combine plots to make the final figures
require(grid)   # for the textGrob() function
#figure 5
figure3 <- ggarrange(p3a, p6a,
                     labels = "AUTO",
                     ncol = 2, nrow = 1,
                     common.legend = TRUE, legend = "right")
figure3 #figure 5 #female=circle, male=triangle, 3n=diamond

if(T){
  png("figures/Figure 5.png",width = 10,height = 5,units = "in",res = 600)
  
  figure3
  
  dev.off()  
}

#figure 6
figure1 <- ggarrange(p1a + rremove("xlab"), p4a + rremove("xlab"), # remove x-axis labels from plots
                    labels = "AUTO",
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "right")
                    #align = "hv", 
                    #font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
figure1 #figure 6 #female=circle, male=triangle, 3n=diamond
annotate_figure(figure1, bottom = textGrob("Acclimation temperature (°C)                 "))

if(T){
  png("figures/Figure 6.png",width = 10,height = 5,units = "in",res = 600)
  
  figure1
  annotate_figure(figure1, bottom = textGrob("Acclimation temperature (°C)                 "))
  
  dev.off()  
}

#figure 7
figure2 <- ggarrange(p2a + rremove("xlab"), p5a + rremove("xlab"), # remove x-axis labels from plots
                     labels = "AUTO",
                     ncol = 2, nrow = 1,
                     common.legend = TRUE, legend = "right")

figure2 #figure 7 #female=circle, male=triangle, 3n=diamond
annotate_figure(figure2, bottom = textGrob("Oxygen condition                  "))

if(T){
  png("figures/Figure 7.png",width = 10,height = 5,units = "in",res = 600)
  
  figure2 #figure 7 #female=circle, male=triangle, 3n=diamond
  annotate_figure(figure2, bottom = textGrob("Oxygen condition                  "))
  
  dev.off()  
}
#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/6.0 Figure_5_6_7_printing_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################