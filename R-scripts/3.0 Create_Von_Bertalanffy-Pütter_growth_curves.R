#create Von Bertalanffy-Pütter growth curves 
rm(list=ls())

#load packages
library(FSAdata) # for data
library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(car)     # for Boot()
library(plyr)
library(dplyr)   # for filter(), mutate()
library(ggplot2)
library(Rmisc)

#load dataset
setwd("/Users/Jboerrigter/Nextcloud/Iris/Paper 3 - Growth paper/Revision/Revision for BiO/Final files/JB")
load(file = "data/newdata3.Rda") #data for all ages

# Multiple VB growth functions in one plot #
vb <- vbFuns()
predict2 <- function(x) predict(x,data.frame(age=ages))

agesum <- group_by(newdata3,cond) %>%
  summarize(minage=min(age),maxage=max(age))
agesum

( conds <- levels(newdata3$cond) ) # To simplify coding below, the levels and number of “groups” are saved into objects.
( nconds <- length(conds) )

cfs <- cis <- preds1 <- preds2 <- NULL # in the loop parameter estimates are saved, first set to NULL

for (i in 1:nconds) {
  ## Loop notification (for peace of mind)
  cat(conds[i],"Loop\n")
  ## Isolate cond's data
  tmp1 <- filter(newdata3,cond==conds[i])
  ## Fit von B to that condition
  sv1 <- vbStarts(tl~age,data=tmp1)
  ## DHO ... do this to remove the t0 starting value (not needed below) ########
  sv1 <- sv1[-3]
  ## DHO ... note that t0 is set to your constant value here ###################
  fit1 <- nls(tl~vb(age,Linf,K,t0=-1.37),data=tmp1,start=sv1)
  ## Extract and store parameter estimates and CIs
  cfs <- rbind(cfs,coef(fit1))
  boot1 <- Boot(fit1)
  tmp2 <-  confint(boot1)
  ## DHO ... I removed the "t0" portion here ###################################
  cis <- rbind(cis,c(tmp2["Linf",],tmp2["K",]))
  ## Predict mean lengths-at-age with CIs
  ##   preds1 -> across all ages
  ##   preds2 -> across observed ages only
  ages <- seq(0,190,0.2)
  boot2 <- Boot(fit1,f=predict2)
  tmp2 <- data.frame(cond=conds[i],age=ages,
                     predict(fit1,data.frame(age=ages)),
                     confint(boot2))
  preds1 <- rbind(preds1,tmp2)
  tmp2 <- filter(tmp2,age>=agesum$minage[i],age<=agesum$maxage[i])
  preds2 <- rbind(preds2,tmp2)
}

rownames(cfs) <- rownames(cis) <- conds
## DHO ... I removed the "t0" portion here #####################################
# colnames(cis) <- paste(rep(c("Linf","K"),each=2),
#                        rep(c("LCI","UCI"),times=2),sep=".")
colnames(preds1) <- colnames(preds2) <- c("cond","age","fit","LCI","UCI")

preds1 # for all ages (observed and predicted)
preds2 # only for observed ages
cfs # parameter estimates
cis # confidence intervals of parameter estimates

#Store outcome
write.table(cfs, file = "Outputs/parameter_estimates.txt", sep = ",", quote = FALSE, row.names = T) # these are the parameter estimates
write.table(cis, file = "Outputs/parameter_estimates_ci.txt", sep = ",", quote = FALSE, row.names = T) # these are the confidence intervals of parameter estimates

save(preds1, file = "data/preds1.Rda")
save(preds2, file = "data/preds2.Rda")

load(file = "data/preds1.Rda") #data for all ages (observed and predicted)
load(file = "data/preds2.Rda") #data for observed ages

################################################################################
#create Von Bertalanffy-Pütter growth curves 
# Change plot (low T's dashed lines)
vbFitPlot2 <- ggplot() + 
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI,fill=cond),alpha=0.2) +
  geom_point(data=newdata3,aes(y=tl,x=age,color=cond),alpha=0.25,size=2,
             position=position_dodge(width=2.5)) + # to prevent overlapping points
  geom_line(data=preds1[preds1$cond %in% c("2n 23", "3n 23"),], aes(y=fit,x=age,color=cond, linetype=cond), linewidth=1) + #put linetype in the aes here, specify in scale_linetype_manual
  geom_line(data=preds1[preds1$cond %in% c("2n 28", "3n 28"),], aes(y=fit,x=age,color=cond, linetype=cond),linewidth=1) +
  scale_color_manual(values=c('dodgerblue1', 'dodgerblue3', 'indianred1', 'indianred3'),
                     aesthetics=c("fill","color")) +
  scale_y_continuous(name="Standard length (mm)",limits=c(0,30),expand=c(0,0)) +
  scale_x_continuous(name="Days post fertilization (dpf)",expand=c(0,0),
                     limits=c(0,190),breaks=seq(0,190,50)) +
  scale_linetype_manual(name = "Condition", values = c("2n 23" = "dotted", "3n 23" = "dotted", "2n 28" = "solid", "3n 28" = "solid")) +
  labs(color = "Condition") +
  guides(fill = "none") + #removes the legend for fill (used for geom_ribbon)
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position=c(0.9,0.12))
vbFitPlot2

  
# Add vertical lines for size at maturity
vbFitPlot2 +
  geom_vline(xintercept=c(98, 126), linetype='dashed', color='gray85') +
  annotate("text", x=98, y=1, label= "Maturity 28°C") +
  annotate("text", x=126, y=2.5, label= "Maturity 23°C")

# Change colors of the plot
vbFitPlot2 <- ggplot() + 
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI,fill=cond),alpha=0.2) +
  geom_point(data=newdata3,aes(y=tl,x=age,color=cond, shape=as.factor(newdata3$temp)),alpha=0.25,size=2,
             position=position_dodge(width=2.5)) + # to prevent overlapping points
  geom_line(data=preds1[preds1$cond %in% c("2n 23", "3n 23"),], aes(y=fit,x=age,color=cond, linetype=cond), linewidth=1) + #put linetype in the aes here, specify in scale_linetype_manual
  geom_line(data=preds1[preds1$cond %in% c("2n 28", "3n 28"),], aes(y=fit,x=age,color=cond, linetype=cond),linewidth=1) +
  scale_color_manual(values=c('#295E11', '#295E11', '#d9bc00', '#d9bc00'), 
                     aesthetics=c("fill","color")) +
  scale_y_continuous(name="Standard length (mm)",limits=c(0,30),expand=c(0,0)) +
  scale_x_continuous(name="Days post fertilization (dpf)",expand=c(0,0),
                     limits=c(0,190),breaks=seq(0,190,50)) +
  scale_linetype_manual(name = "Condition", values = c("2n 23" = "dotted", "3n 23" = "dotted", "2n 28" = "solid", "3n 28" = "solid")) +
  scale_shape_manual(values = c(17,16)) +
  labs(color = "Condition") +
  guides(fill = "none") + #removes the legend for fill (used for geom_ribbon)
  guides(shape = "none") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position=c(0.9,0.12))
vbFitPlot2

#levels(as.factor(newdata3$temp))

# Add vertical lines for size at maturity
if(T){
png("figures/Figure 2.png",width = 10,height = 5,units = "in",res = 600)
vbFitPlot2 +
  geom_vline(xintercept=c(98, 126), linetype='dashed', color='gray85') +
  annotate("text", x=98, y=1, label= "Maturity 28°C") +
  annotate("text", x=126, y=2.5, label= "Maturity 23°C")

dev.off()  
}
#------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("outputs/3.0 Create_Von_Bertalanffy-Pütter_growth_curves_session.txt")
sessionInfo()
sink() 
################################################################################
###########################           End          #############################
################################################################################