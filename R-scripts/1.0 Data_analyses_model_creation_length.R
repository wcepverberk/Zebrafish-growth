rm(list=ls())

#load packages
library(ggplot2)
library(Rmisc)
library(cowplot)
library(visreg)
library(lmerTest)
library(ggpubr)
library(MuMIn)
library(AICcmodavg)

#load dataset 1
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
data <- read.csv(file= "data/datafish1-114.csv", header = TRUE)

#Check and organize data structure
str(data)
data$ox <- as.factor(data$ox)
data$ploidy <- as.factor(data$ploidy)
data$sex <- as.factor(data$sex)
data$batch <- as.factor(data$batch)

head(data)

# convert weeks to days
data$age_days <- data$age * 7 
data
data$treat <- as.factor(paste(data$r_temp, data$a_temp, data$ploidy, data$ox))
levels(data$treat)

# Calculate average length of fish of a specific rearing T and ploidy level at the start of acclimation #
#load dataset 2
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


# Calculate Lenght at maturity
data2n28 <- subset(dat, dat$age=="98"&dat$cond=="2n 28")
data3n28 <- subset(dat, dat$age=="98"&dat$cond=="3n 28")
data2n23 <- subset(dat, dat$age=="126"&dat$cond=="2n 23")
data3n23 <- subset(dat, dat$age=="126"&dat$cond=="3n 23")
data3n23.2 <- subset(dat, dat$age=="119"&dat$cond=="3n 23"&dat$fdate=="28-1-2021") #B6 reared as extra batch, but not measured at 126 days
mean(data2n28$tl)
mean(data3n28$tl, na.rm = TRUE)
mean(data2n23$tl, na.rm = TRUE)
mean(data3n23$tl, na.rm = TRUE)
mean(data3n23.2$tl, na.rm = TRUE) #length at maturity differs slightly. Length at start acclimation differs more (larger difference in age)

#length at start acclimation
mean(data2n28$tl) #length at maturity same as length at start acclimation for batch 1
data3n28sa <- subset(dat, dat$age=="105"&dat$cond=="3n 28") # acclimation started one week later for batch 3
mean(data3n28sa$tl, na.rm = TRUE)
data2n23sa <- subset(dat, dat$age=="133"&dat$cond=="2n 23") # acclimation started one week later for batch 4
mean(data2n23sa$tl)
data3n23sa <- subset(dat, dat$age=="133"&dat$cond=="3n 23") #acclimation started one week later for batch 2
mean(data3n23sa$tl, na.rm = TRUE)
mean(data3n23.2$tl, na.rm = TRUE) #acclimation started at 17 weeks instead of 18 for bacth 6. Fish were classified as adults this week

# calcuate average length at start acclimation
data$l1 <- 0 #
data$l1[which(data$batch=="b1")] <- paste(mean(data2n28$tl))
data$l1[which(data$batch=="b3")] <- paste(mean(data3n28sa$tl, na.rm = TRUE))
data$l1[which(data$batch=="b4")] <- paste(mean(data2n23sa$tl))
data$l1[which(data$batch=="b2")] <- paste(mean(data3n23sa$tl, na.rm = TRUE))
data$l1[which(data$batch=="b6")] <- paste(mean(data3n23.2$tl, na.rm = TRUE))

str(data)
data$ploidy <-as.factor(data$ploidy)
data$ox <-as.factor(data$ox)
data$sex <-as.factor(data$sex)
data$batch <-as.factor(data$batch)
data$l1 <- as.numeric(data$l1)

################################################################################
# Model selection #

m0 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp+ox+log10(l1))+(1|batch), na.action = na.exclude, data=data) 
summary(m0)
anova(m0)
AIC(m0)

m1 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp+ox)+ploidy*log10(l1)+(1|batch), na.action = na.exclude, data=data) # - l1 3-way interaction
summary(m1)
anova(m1)
AIC(m1) # same AIC as m1
r.squaredGLMM(m1)

m2 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp+log10(l1))+ploidy*ox+(1|batch), na.action = na.exclude, data=data) # - ox 3-way interaction
summary(m2)
anova(m2)
AIC(m2)
r.squaredGLMM(m2)

m3 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(r_temp+log10(l1)+ox)+ploidy*a_temp+(1|batch), na.action = na.exclude, data=data) # -a_temp 3-way interaction
summary(m3)
anova(m3)
AIC(m3)
r.squaredGLMM(m3)

m4 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(a_temp+log10(l1)+ox)+ploidy*r_temp+(1|batch), na.action = na.exclude, data=data) # -r_temp 3-way interaction
summary(m4)
anova(m4)
AIC(m4)
r.squaredGLMM(m4)

m5 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp)+ploidy*(log10(l1)+ox)+(1|batch), na.action = na.exclude, data=data) # - l1 - ox 3-way interaction
summary(m5)
anova(m5)
AIC(m5)
r.squaredGLMM(m5)

m6 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*r_temp+ploidy*(a_temp+log10(l1)+ox)+(1|batch), na.action = na.exclude, data=data) # - l1 - ox - a_temp 3-way interaction
summary(m6)
anova(m6)
AIC(m6)
r.squaredGLMM(m6)

m7 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(r_temp+a_temp+log10(l1)+ox)+(1|batch), na.action = na.exclude, data=data) # - l1 - ox - a_temp - r_temp 3-way interaction
summary(m7)
anova(m7)
AIC(m7)
r.squaredGLMM(m7)

m8 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(r_temp+a_temp+log10(l1))+ox+(1|batch), na.action = na.exclude, data=data) # - ox 2-way interaction
summary(m8)
anova(m8)
AIC(m8)
r.squaredGLMM(m8)

m9 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(r_temp+a_temp+ox)+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - l1 2-way interaction (same as m11)
summary(m9)
anova(m9)
AIC(m9)
r.squaredGLMM(m9)

m10 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(r_temp+ox)+a_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - a_temp 2-way interaction
summary(m10)
anova(m10)
AIC(m10)
r.squaredGLMM(m10)

m11 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(ox+a_temp)+r_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - r_temp 2-way interaction
summary(m11)
anova(m11)
AIC(m11)
r.squaredGLMM(m11)

m12 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy:ox+a_temp+r_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - r_temp - a_temp 2-way interaction
summary(m12)
anova(m12)
AIC(m12)
r.squaredGLMM(m12)

m13 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ox+a_temp+r_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - all 2-way interactions
summary(m13)
anova(m13)
AIC(m13)
r.squaredGLMM(m13)

m13a <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ox+a_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - all 2-way interactions
summary(m13a)
anova(m13a)
AIC(m13a)
r.squaredGLMM(m13a)


m14 <-lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy:ox+a_temp+r_temp+a_temp:r_temp+log10(l1)+(1|batch), na.action = na.exclude, data=data) # - r_temp - a_temp 2-way interaction, + a_temp:r_temp interaction
summary(m14)
anova(m14)
AIC(m14)
r.squaredGLMM(m14)

m15 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")+ploidy*(r_temp+a_temp+log10(l1)+ox)+a_temp:r_temp+(1|batch), na.action = na.exclude, data=data) # - l1 - ox - a_temp - r_temp 3-way interaction
summary(m15)
anova(m15)
AIC(m15)
r.squaredGLMM(m15)

m16 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*ox+a_temp+r_temp+log10(l1)+ploidy:ox:r_temp+(1|batch), na.action = na.exclude, data=data) # + 3-way interaction ploidy:a_temp:ox
summary(m16)
anova(m16)
AIC(m16)
r.squaredGLMM(m16)

m17 <- lmer(log10(length)~ploidy/relevel(sex, ref="m")*ox+a_temp+r_temp+log10(l1)+ploidy:ox:a_temp+(1|batch), na.action = na.exclude, data=data) # + 3-way interaction ploidy:a_temp:ox
summary(m17)
anova(m17)
AIC(m17)
r.squaredGLMM(m17)

#Table of created models
if(T){
  
  #make model list
  fit.list.table.1.m <- list(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17)
  
  # give model names
  fit.names.table.1.m <-c(
    "m0, ploidy/relevel(sex=male)*(r_temp+a_temp+ox+log10(l1))+(1|batch)",
    "m1, ploidy/relevel(sex=male)*(r_temp+a_temp+ox)+ploidy*log10(l1)+(1|batch)",
    "m2, ploidy/relevel(sex=male)*(r_temp+a_temp+log10(l1))+ploidy*ox+(1|batch)",
    "m3, ploidy/relevel(sex=male)*(r_temp+log10(l1)+ox)+ploidy*a_temp+(1|batch)",
    "m4, ploidy/relevel(sex=male)*(a_temp+log10(l1)+ox)+ploidy*r_temp+(1|batch)",
    "m5, ploidy/relevel(sex=male)*(r_temp+a_temp)+ploidy*(log10(l1)+ox)+(1|batch)",
    "m6, ploidy/relevel(sex=male)*r_temp+ploidy*(a_temp+log10(l1)+ox)+(1|batch)",
    "m7, ploidy/relevel(sex=male)+ploidy*(r_temp+a_temp+log10(l1)+ox)+(1|batch)",
    "m8, ploidy/relevel(sex=male)+ploidy*(r_temp+a_temp+log10(l1))+ox+(1|batch)",
    "m9, ploidy/relevel(sex=male)+ploidy*(r_temp+a_temp+ox)+log10(l1)+(1|batch)",
    "m10, ploidy/relevel(sex=male)+ploidy*(r_temp+ox)+a_temp+log10(l1)+(1|batch)",
    "m11, ploidy/relevel(sex=male)+ploidy*(ox+a_temp)+r_temp+log10(l1)+(1|batch)",
    "m12, ploidy/relevel(sex=male)+ploidy:ox+a_temp+r_temp+log10(l1)+(1|batch)",
    "m13, ploidy/relevel(sex=male)+ox+a_temp+r_temp+log10(l1)+(1|batch)",
    "m14, ploidy/relevel(sex=male)+ploidy:ox+a_temp+r_temp+a_temp:r_temp+log10(l1)+(1|batch)",
    "m15, ploidy/relevel(sex=male)+ploidy*(r_temp+a_temp+log10(l1)+ox)+a_temp:r_temp+(1|batch)",
    "m16, ploidy/relevel(sex=male)*ox+a_temp+r_temp+log10(l1)+ploidy:ox:r_temp+(1|batch)", 
    "m17, ploidy/relevel(sex=male)*ox+a_temp+r_temp+log10(l1)+ploidy:ox:a_temp+(1|batch)"
      )
  
  ###compare by using AIC
  fit.table.1.m<-aictab(fit.list.table.1.m,fit.names.table.1.m, second.ord = F,sort = TRUE, digits = 3, LL=TRUE)
  fit.table.1.m
  
  r2<-lapply(fit.list.table.1.m, r.squaredGLMM)
  rsq_table <- do.call(rbind, r2)
  
  # Add model names if needed
  rsq_table2 <- data.frame(Model = fit.names.table.1.m, rsq_table)
  rsq_table2
  
  # Ensure the output is a proper data frame
  r2_n <- data.frame(
    Model = fit.names.table.1.m,
    R2m = sapply(r2, function(x) x[1]),  
    R2c = sapply(r2, function(x) x[2])   
  )
  
  #rename column for merging tables
  fit.table.1.m <- fit.table.1.m %>%
    dplyr::rename(Model = Modnames)
  
  # Merge the tables by "Model"
  fit.table.1.m <- merge(fit.table.1.m, r2_n, by = "Model")
  
  # Reorder by AIC
  fit.table.1.m <- fit.table.1.m[order(fit.table.1.m$AIC), ]
  fit.table.1.m
  
  ## Export tables to .csv ----
  write.csv(fit.table.1.m,"outputs/1.0 Table_s1_Model_comparison_for_length.csv",row.names = FALSE)
}

# make table of the selected model
# Table S1-S2 ----
Table_S1 <- capture.output(summary(m13)) 
Table_S2 <- capture.output(anova(m13))

write.table(Table_S1,"outputs/1.0 Table_summary of the model with the highest support for length.txt",row.names = F,quote = F)
write.table(Table_S2,"outputs/1.0 Table_anova of the model with the highest support for length.txt",row.names = F,quote = F)

#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/1.0 Data_analyses_model_creation_length_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################