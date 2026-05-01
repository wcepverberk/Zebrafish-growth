rm(list=ls())

#load packages
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(lme4)
library(lmerTest)
library(nlme)

#load dataset 1
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
getwd()
dat <- read.csv(file= "data/Groeidata R.csv", header = TRUE, sep = ";", dec = "." )

#check and organize data 
str(dat)

names(dat)[1] <- "tank"
names(dat)[2] <- "condition"
dat$condition<-as.factor(dat$condition)
names(dat)[3] <- "temp"
names(dat)[4] <- "ploidy"
names(dat)[5] <- "fdate"
names(dat)[6] <- "fish"
names(dat)[7] <- "age"
names(dat)[8] <- "tl"
dat$cond <- paste(dat$ploidy, dat$temp)
dat$cond <- as.factor(dat$cond)

head(dat)  

# visualize data
plot(tl~age,
     col = factor(dat$cond),
     xlab = "Days post fertilization (dpf)", 
     ylab = "Length (mm)",
     data=dat)
legend("topleft",
       legend = levels(factor(dat$cond)),
       pch = 1,
       col = factor(levels(factor(dat$cond))))

#load dataset 1
data <- read.csv(file= "data/datafish1-114.csv", header = TRUE)

#check and organize data 
str(data)
data$ox <- as.factor(data$ox)
data$ploidy <- as.factor(data$ploidy)
data$sex <- as.factor(data$sex)
data$batch <- as.factor(data$batch)

head(data)

# convert weeks to days
data$age_days <- data$age * 7  

# Create new data frames and combine them #
newdata <- data.frame(data$fish_id[which(data$ox=="norm")], data$r_temp[which(data$ox=="norm")],
                      data$ploidy[which(data$ox=="norm")], data$length[which(data$ox=="norm")],
                      data$f_date[which(data$ox=="norm")], data$age_days[which(data$ox=="norm")])
names(newdata) <- c("fish", "temp", "ploidy", "tl", "fdate", "age")  
newdata$tank <- "NA"
head(newdata)
newdata2 <- data.frame(dat$fish, dat$temp, dat$ploidy, dat$tl, dat$fdate, dat$age, dat$tank)
names(newdata2) <- c("fish", "temp", "ploidy", "tl", "fdate", "age", "tank")
newdata2

newdata3 <- rbind(newdata, newdata2)  
newdata3$cond <- paste(newdata3$ploidy, newdata3$temp)
newdata3$cond <- as.factor(newdata3$cond)
head(newdata3)

str(newdata3)
newdata3$age<-as.integer(newdata3$age)
newdata3<-subset(newdata3, newdata3$tl!="NA")

# visualize data in new plot with extra length data #
plot(tl~age,
     col = factor(newdata3$cond),
     xlab = "Days post fertilization (dpf)", 
     ylab = "Length (mm)",
     data=newdata3)
legend("topleft",
       legend = levels(factor(newdata3$cond)),
       pch = 1,
       col = factor(levels(factor(newdata3$cond))))

# size at maturity #
data2n28 <- subset(newdata3, newdata3$age=="98"&newdata3$cond=="2n 28")
data3n28 <- subset(newdata3, newdata3$age=="98"&newdata3$cond=="3n 28")
data2n23 <- subset(newdata3, newdata3$age=="126"&newdata3$cond=="2n 23")
data3n23 <- subset(newdata3, newdata3$age=="126"&newdata3$cond=="3n 23")

sizemat <- rbind(data2n28, data3n28, data2n23, data3n23)
sizemat$temp <-as.factor(sizemat$temp)

sizemats <- summarySE(sizemat, measurevar = "tl", groupvars = c("ploidy", "temp"), na.rm = TRUE)
sizemats

# visualize data in a Bar plot
ggplot(sizemats, aes(x= temp, y= tl, fill = ploidy)) + #set of aesthetic mappings, fill is the grouping variable to fill the bars with different colours
  geom_col(position=position_dodge(), colour="black", size = .3, width=.5, show.legend = TRUE) + #colour and size ar about the outline of the bars, position is to not put them stacked
  geom_errorbar(aes(ymin=tl-se, ymax=tl+se), position=position_dodge(.5), size = .3, width=.2) + #with position_dodge() the errorbars are moved horizontally
  xlab("Temperature") +
  ylab("Length") +
  scale_fill_manual(name="Ploidy level", values = c("2n" = "dodgerblue3", "3n" = "indianred3")) + #scale_fill_hue can be used for automatic fill colour per category
  theme_bw()

# visualize data in a Box plots
ggplot(sizemat, aes(x= temp, y= tl, fill = ploidy)) + #set of aesthetic mappings, fill is the grouping variable to fill the bars with different colours
  geom_boxplot() + 
  xlab("Temperature") +
  ylab("Length") +
  scale_fill_manual(name="Ploidy level", values = c("2n" = "dodgerblue3", "3n" = "indianred3")) + #scale_fill_hue can be used for automatic fill colour per category
  theme_bw() 

p_facet <- ggplot(sizemat, aes(x= temp, y= tl, fill = ploidy)) + #set of aesthetic mappings, fill is the grouping variable to fill the bars with different colours
  geom_boxplot() + 
  xlab("Temperature") +
  ylab("Length") +
  scale_fill_manual(name="Ploidy level", values = c("2n" = "dodgerblue3", "3n" = "indianred3")) + #scale_fill_hue can be used for automatic fill colour per category
  theme_bw() +
  facet_wrap(~ploidy)
p_facet

p2n <- ggplot(subset(sizemat, ploidy %in% "2n"), aes(x= temp, y= tl, fill = ploidy)) + #set of aesthetic mappings, fill is the grouping variable to fill the bars with different colours
  geom_boxplot() + 
  xlab("Temperature") +
  ylab("Length") +
  scale_fill_manual(name="Ploidy level", values = c("2n" = "dodgerblue3", "3n" = "indianred3")) + #scale_fill_hue can be used for automatic fill colour per category
  theme_bw() 

p2n

p3n <- ggplot(subset(sizemat, ploidy %in% "3n"), aes(x= temp, y= tl, fill = ploidy)) + #set of aesthetic mappings, fill is the grouping variable to fill the bars with different colours
  geom_boxplot() + 
  xlab("Temperature") +
  ylab("Length") +
  scale_fill_manual(name="Ploidy level", values = c("2n" = "dodgerblue3", "3n" = "indianred3")) + #scale_fill_hue can be used for automatic fill colour per category
  theme_bw() 

p3n

# load stored data of Von_Bertalanffy-Pütter_growth_curves
dpe <- read.csv(file= "data/parameter_estimates_def.csv", header = TRUE)
dpeci <- read.csv(file= "data/parameter_estimates_ci_def.csv", header = TRUE)

dpe
dpeci

#combine datasets 
linf <- cbind(dpe, dpeci)
linf
linf <- linf[,-4]
linf
linf$temp <- c("23", "28", "23", "28")
linf$ploidy <- c("2n", "2n", "3n", "3n")
linf

write.table(linf, file = "Outputs/Linf.txt", sep = ",", quote = FALSE, row.names = T)

################################################################################
#create figure 4
### Change ploidy colors ###

p_2n <- ggplot() +
  geom_boxplot(data = subset(sizemat, ploidy %in% "2n"), aes(x = temp, y = tl, fill = cond)) + ylim(13,50) +
  geom_ribbon(data = subset(linf, cond %in% c("2n 23", "2n 28")), aes(x = temp,ymin = Linf.LCI,ymax = Linf.UCI, group = 1), fill = "#295E11", alpha = 0.2) + 
  geom_point(data = subset(linf, cond %in% c("2n 23", "2n 28")), aes(x = temp, y = Linf, color = cond), size = 2) +
  geom_line(data = subset(linf, cond %in% c("2n 23", "2n 28")), aes(x = temp, y = Linf, group = 1, linetype = "cond"), show.legend = FALSE) + #group consists of one observation: group = 1
  xlab("Temperature (°C)") +
  ylab("Length (mm)") +
  scale_fill_manual(name = "Condition", values = c("2n 28" = "#295E11", "2n 23" = "#295E11", "3n 28" = "#f4dc01", "3n 23" = "#f4dc01")) + #scale_fill_hue can be used for automatic fill colour per category
  scale_color_manual(name = "Condition", values = c("2n 28" = "#295E11", "2n 23" = "#295E11", "3n 28" = "#f4dc01", "3n 23" = "#f4dc01")) + #otherwise point color not changed (one dimentional geom = color, two dimensional = color + fill)
  #adding all the legend values here as we combine the plots later and want one shared legend
  theme_bw()

p_2n

p_3n <- ggplot() +
  geom_boxplot(data = subset(sizemat, ploidy %in% "3n"), aes(x = temp, y = tl, fill = cond)) + ylim(13,50) +
  geom_ribbon(data = subset(linf, cond %in% c("3n 23", "3n 28")), aes(x = temp,ymin = Linf.LCI,ymax = Linf.UCI, group = 1), fill = "#f4dc01", alpha = 0.2) +
  geom_point(data = subset(linf, cond %in% c("3n 23", "3n 28")), aes(x = temp, y = Linf, color = cond), size = 2) +
  geom_line(data = subset(linf, cond %in% c("3n 23", "3n 28")), aes(x = temp, y = Linf, group = 1, linetype = "cond"), show.legend = FALSE) +
  xlab("Temperature (°C)") +
  ylab("  ") +
  scale_fill_manual(name = "Condition", values = c("2n 28" = "#295E11", "2n 23" = "#295E11", "3n 28" = "#f4dc01", "3n 23" = "#f4dc01")) + #scale_fill_hue can be used for automatic fill colour per category
  scale_colour_manual(name = "Condition", values = c("2n 28" = "#295E11", "2n 23" = "#295E11", "3n 28" = "#f4dc01", "3n 23" = "#f4dc01")) + #otherwise point color not changed
  theme_bw()

p_3n

# combine p_2n and p_3n into Figure 4
require(grid)
figure <- ggarrange(p_2n + rremove("xlab"), p_3n + rremove("xlab"), # remove x-axis labels from plots
                    labels = "AUTO",
                    ncol = 2, nrow = 1,
                    common.legend = FALSE, legend = "right")

figure# figure 3

annotate_figure(figure, bottom = textGrob("Temperature (°C)        "))

if(T){
  png("figures/Figure 3.png",width = 10,height = 5,units = "in",res = 600)
  figure# figure 3
  
  annotate_figure(figure, bottom = textGrob("Temperature (°C)        "))
  
  dev.off()  
}
################################################################################
#model selection
#check and organize data
str(sizemat)
sizemat$tank <- as.factor(sizemat$tank)
sizemat$tl_log<-log10(sizemat$tl)

#Create models
lm<-gls(tl ~ 1,data=sizemat, method = "REML", na.action = na.omit)
lme<-lme(tl ~ 1,random=~1|fish,data=sizemat,method = "REML", na.action = na.omit)
anova(lm,lme) #mixed model is not better

m1 <- lm(tl~ploidy+temp+fish, data = sizemat) #additive
summary(m1)
anova(m1)
AIC(m1)

m2 <- lm(tl~ploidy*temp+fish, data = sizemat) # interaction
summary(m2)
anova(m2)
AIC(m2)

m3 <- lm(tl~ploidy+temp, data = sizemat) # additive
summary(m3)
anova(m3)
AIC(m3)

m4 <- lm(tl~ploidy*temp, data = sizemat) # interaction
summary(m4)
anova(m4)
AIC(m4)

m5 <- lm(tl_log~ploidy*temp+fish, data = sizemat) # interaction
summary(m5)
anova(m5)
AIC(m5)

m6 <- lm(tl_log~ploidy+temp, data = sizemat) # interaction
summary(m6)
anova(m6)
AIC(m6)

m7 <- lm(tl_log~ploidy+temp+fish, data = sizemat) # interaction
summary(m7)
anova(m7)
AIC(m7)

#Table of model selection
if(T){
  #create model list
  fit.list.table.1.m <- list(m1,m2,m3,m4)
  
  #create model names list
  fit.names.table.1.m <-c("m1, ploidy+temp+fish",
                          "m2, ploidy*temp+fish",
                          "m3, ploidy+temp",
                          "m4, ploidy*temp")
  
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
  write.csv(fit.table.1.m,"outputs/4.0 Table_Model_comparison_for_figure 3 statistics.csv",row.names = FALSE)
}

# Output summary and anova statistics in Table S1-S2 ----
Table_S1 <- capture.output(summary(m4))
Table_S2 <- capture.output(anova(m4))

write.table(Table_S1,"outputs/4.0 Table_summary of the model for figure 3.txt",row.names = F,quote = F)
write.table(Table_S2,"outputs/4.0 Table_anova of the model with the highest support for figure 3.txt",row.names = F,quote = F)
#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/4.0 Figure_3_printing_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################

