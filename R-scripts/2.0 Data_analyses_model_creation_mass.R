rm(list=ls())

#load packages
library(ggplot2)
library(Rmisc)
library(cowplot)
library(visreg)
library(lmerTest)
library(MuMIn)
library(AICcmodavg)
library(dplyr)

#Load dataset 1
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
data <- read.csv(file= "Data/datafish1-114.csv", header = TRUE)

#check and reorganize data structure
str(data)
data$ox <- as.factor(data$ox)
data$ploidy <- as.factor(data$ploidy)
data$sex <- as.factor(data$sex)
data$batch <- as.factor(data$batch)
head(data)

# weeks to days
data$age_days <- data$age * 7 # weeks to days
data
data$treat <- as.factor(paste(data$r_temp, data$a_temp, data$ploidy, data$ox))
levels(data$treat)

#load dataset 2
getwd()
dat <- read.csv(file= "Data/Groeidata R.csv", header = TRUE, sep = ";", dec = "." )

#check and reorganize data structure
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

# Lenght at maturity
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
mean(data2n28$tl) #length at maturity same as length at start acclimation for bacth 1
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

str(data)
data$l1 <- as.numeric(data$l1)
head(data)

###############    model creating and selection     ############################
m0 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+r_temp+a_temp+ox)+(1|batch), na.action = na.exclude, data=data)
m1 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+r_temp+a_temp)+ox+(1|batch), na.action = na.exclude, data=data) 
m2 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+r_temp+ox)+a_temp+(1|batch), na.action = na.exclude, data=data) 
m3 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+a_temp+ox)+r_temp+(1|batch), na.action = na.exclude, data=data) 
m4 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp)+log10(length)+ox+(1|batch), na.action = na.exclude, data=data) 
m5 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+r_temp)+ox+a_temp+(1|batch), na.action = na.exclude, data=data) 
m6 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+a_temp)+r_temp+ox+(1|batch), na.action = na.exclude, data=data) 
m7 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(a_temp+ox)+log10(length)+r_temp+(1|batch), na.action = na.exclude, data=data) 
m8 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*ox+(log10(length)+r_temp+a_temp)+(1|batch), na.action = na.exclude, data=data) 
m9 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*a_temp+(log10(length)+r_temp+ox)+(1|batch), na.action = na.exclude, data=data) 
m10 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*r_temp+(log10(length)+a_temp+ox)+(1|batch), na.action = na.exclude, data=data)
m11 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*ox+(r_temp+a_temp)+log10(length)+(1|batch), na.action = na.exclude, data=data)
m12 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")+ox+(r_temp+a_temp)+log10(length)+(1|batch), na.action = na.exclude, data=data)
m13 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp+log10(length)+ox+ox:r_temp)+(1|batch), na.action = na.exclude, data=data)
m14 <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(r_temp+a_temp+log10(length)+ox+ox:a_temp)+(1|batch), na.action = na.exclude, data=data)

AIC(m0)
AIC(m1)
AIC(m2)
AIC(m3) #-209.0237 more complex
AIC(m4)
AIC(m5)
AIC(m6) #-210.203  Best model
AIC(m7)
AIC(m8)
AIC(m9)
AIC(m10)
AIC(m11)
AIC(m12)
AIC(m13)
AIC(m14)

r.squaredGLMM(m3) #R2m 0.9112511, R2c 0.9112511
r.squaredGLMM(m6) #R2m 0.9002051, R2c 0.9002051

#model refinement
m3  <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(log10(length)+a_temp+ox)+r_temp+(1|batch), na.action = na.exclude, data=data) 
m3a <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(a_temp+ox)+log10(length)+r_temp+(1|batch), na.action = na.exclude, data=data) 
m3b <- lmer(log10(weight)~ploidy/relevel(sex, ref="m")*(a_temp)+ox+log10(length)+r_temp+(1|batch), na.action = na.exclude, data=data) 
summary(m3a)
summary(m3b)
AIC(m3a)   #-195.1608 Lesser model
AIC(m3b) # -194.033 Lesser model


#Table of model selection
if(T){
  
 # make model list
  fit.list.table.1.m <- list(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
  
  # make model names list
  fit.names.table.1.m <-c(
                          "m0, ploidy/relevel(sex=male)*(log10(length)+r_temp+a_temp+ox)+(1|batch)",
                          "m1, ploidy/relevel(sex=male)*(log10(length)+r_temp+a_temp)+ox+(1|batch)",
                          "m2, ploidy/relevel(sex=male)*(log10(length)+r_temp+ox)+a_temp+(1|batch)",
                          "m3, ploidy/relevel(sex=male)*(log10(length)+a_temp+ox)+r_temp+(1|batch)",
                          "m4, ploidy/relevel(sex=male)*(r_temp+a_temp)+log10(length)+ox+(1|batch)",
                          "m5, ploidy/relevel(sex=male)*(log10(length)+r_temp)+ox+a_temp+(1|batch)",
                          "m6, ploidy/relevel(sex=male)*(log10(length)+a_temp)+r_temp+ox+(1|batch)",
                          "m7, ploidy/relevel(sex=male)*(a_temp+ox)+log10(length)+r_temp+(1|batch)",
                          "m8, ploidy/relevel(sex=male)*ox+(log10(length)+r_temp+a_temp)+(1|batch)",
                          "m9, ploidy/relevel(sex=male)*a_temp+(log10(length)+r_temp+ox)+(1|batch)",
                          "m10, ploidy/relevel(sex=male)*r_temp+(log10(length)+a_temp+ox)+(1|batch)",
                          "m11, ploidy/relevel(sex=male)*ox+(r_temp+a_temp)+log10(length)+(1|batch)",
                          "m12, ploidy/relevel(sex=male)+ox+(r_temp+a_temp)+log10(length)+(1|batch)",
                          "m13, ploidy/relevel(sex=male)*(r_temp+a_temp+log10(length)+ox+ox:r_temp)+(1|batch)",
                          "m14, ploidy/relevel(sex=male)*(r_temp+a_temp+log10(length)+ox+ox:a_temp)+(1|batch)"
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
  write.csv(fit.table.1.m,"outputs/2.0 Table_s2_Model_comparison_for_Mass.csv",row.names = FALSE)
}
 
# Model m3 is the selected model
# Table S1-S2 ----
Table_S1 <- capture.output(summary(m3))
Table_S2 <- capture.output(anova(m3))

write.table(Table_S1,"outputs/2.0 Table_summary of the model with the highest support for mass.txt",row.names = F,quote = F)
write.table(Table_S2,"outputs/2.0 Table_anova of the model with the highest support for mass.txt",row.names = F,quote = F)

#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/2.0 Data_analyses_model_creation_mass_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################