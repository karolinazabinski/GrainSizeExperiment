###########
##SETTING UP WORKSPACE##
###########
#remove objects from the global environment and set wd
rm(list=ls())
setwd("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment")

#libraries needed for data manipulation
library(dplyr) #manipulate df
library(tidyverse)
#loading packages needed for data analysis
library(nlme) 
library(lme4)
library(doBy) #creates summary tables
library(MuMIn)
library(visreg)
library(car) #has Anova() function
library(emmeans)
#and visualization
library(ggplot2)
library(gridExtra)

##############################
##MANIPULATING METADATA
##############################
meta <-  read.csv("GSEData.csv", header=T, dec=".")
str(meta) #most things in class they should be
#convert tank to factor
meta$TankID <- as.factor(meta$TankID)
str(meta)
#convert Dead:LeafGrowth5_mm from int to numeric
cols <- names(meta)[5:12]
meta[cols] <- lapply(meta[cols], as.numeric)
str(meta)
#convert DaysGrowing:No_RootBundles to numeric
cols2 <- names(meta)[15:23]
meta[cols2] <- lapply(meta[cols2], as.numeric)
str(meta)

##CALCULATIONS
#sums of biomass
meta <- meta %>% mutate(term.g = Term_FoilPost_g -Term_FoilPre_g) %>%
  mutate(side.g = Side_FoilPost_g - Side_FoilPre_g) %>% mutate(rhiz.g =Rhiz_FoilPost_g - Rhiz_FoilPre_g) %>%
  mutate(root.g = Root_FoilPost_g-Root_FoilPre_g)
str(meta)
#above, below, A:B
meta <- meta %>% mutate(sumAbove = rowSums(select(meta, term.g, side.g), na.rm=TRUE))
meta$sumAbove[meta$sumAbove == 0] <- NA

meta <- meta %>% mutate(sumBel = rowSums(select(meta, rhiz.g, root.g), na.rm=TRUE))
meta$sumBel[meta$sumBel == 0] <- NA

meta <- meta %>% mutate(ab = sumAbove/sumBel)

meta <- meta %>% mutate(totalBio = sumAbove + sumBel)

#average root hair length
meta <- meta %>% mutate(rhLength = (root.hair1.mm + root.hair2.mm + root.hair3.mm)/3)

#total growth and total growth rate and productivity
meta$totalGrowth <-as.numeric(rowSums(meta[, c("LeafGrowth1_mm", "LeafGrowth2_mm", "LeafGrowth3_mm", "LeafGrowth4_mm", "LeafGrowth5_mm")], na.rm = TRUE))
meta$totalGrowth[meta$totalGrowth == 0] <- NA

meta <- meta %>% mutate(growthRate = totalGrowth/DaysGrowing)
meta$growthRate[meta$growthRate == 0] <- NA

meta$productivity <- as.numeric(meta$totalGrowth*meta$ShootWidth_mm/meta$DaysGrowing/100) #all NAs remained NAs

#No SS
meta <- meta %>% mutate(NoSS = ShootCount-1)
meta$NoSS[meta$NoSS == 0] <- NA

#average ss side
meta <- meta %>% mutate(avgSsSize = side.g/NoSS)

##creating col for pop-treatment
for(i in 1:nrow(meta)){
  if(meta$PlantOrigin[i]=="MILL"&meta$SedOrigin[i]=="MILL"){
    meta$poptreat[i] <- "MILL-S"}
  if(meta$PlantOrigin[i]=="MILL"&meta$SedOrigin[i]=="BLAK"){
    meta$poptreat[i] <- "MILL-O"}
  if(meta$PlantOrigin[i]=="BLAK"&meta$SedOrigin[i]=="BLAK"){
    meta$poptreat[i] <- "BLAK-S"}
  if(meta$PlantOrigin[i]=="BLAK"&meta$SedOrigin[i]=="MILL"){
    meta$poptreat[i] <- "BLAK-O"}
}

##################
##FUNCTION NEEDED FOR GRAPHS
#################
#set function to calculate the mean, standard deviation, length, and standard error of variables
fun <- function(x,...) {c(m=mean(x, na.rm=T), sd=sd(x,na.rm=T), n = length(x), se=sd(x,na.rm=T)/sqrt(length(x)))}


###########################
##BOXPLOT OF PRODUCTIVITY##
##########################
prod <- meta %>% select (SampleID:TankID, productivity, poptreat)

#summary stats
sum.prod <- summaryBy(productivity ~ PlantOrigin + SedOrigin, data=prod, FUN=fun)

#plot
prod.plot <- ggplot(sum.prod, aes(x=SedOrigin, y=productivity.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin)) +
  geom_errorbar(aes(ymin=productivity.m-productivity.se, ymax =productivity.m+productivity.se), width=.1) +
  theme_classic() +
  ylab(expression(paste( "Productivity (cm^2/day)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("darkolivegreen", "brown2"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) +
  ylim(0,5)
prod.plot

#boxplot
#set.seed(5)
#prod.plot <- ggplot(prod, aes(x=PlantOrigin, y=productivity, group=PlantOrigin)) + 
  #geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  #theme_bw() + theme(panel.grid.major=element_blank(),
                     #panel.grid.minor=element_blank())
#prod.plot

#ANALYSIS
hist(prod$productivity)
prod$logProd <- log(prod$productivity)^2
prod.2way  <- lm(logProd ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin, data=prod, na.action = na.exclude)
prod.1way <- lm(logProd ~ PlantOrigin  + SedOrigin, data=prod, na.action = na.exclude )
#see if 2way can drop
anova(prod.2way, prod.1way)  #drop

prod.full <- lm(productivity ~ PlantOrigin  + SedOrigin, data=prod, na.action = na.exclude )
plot(prod.full)
qqnorm(resid(prod.full))
qqline(resid(prod.full))
shapiro.test(resid(prod.full))# passes
summary(prod.full) 
Anova(prod.full, type = 2) #planOrigin

###########################
##BOXPLOT OF AB##
##########################
aboveBelow <- meta %>% select (SampleID:TankID, ab, poptreat)
#setting seed so that code is run from same random path
set.seed(5)
ab.plot <- ggplot(aboveBelow, aes(x=PlantOrigin, y=ab, group=PlantOrigin)) + 
  geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  theme_bw()
ab.plot

#ANALYSIS
ab.2way  <- lmer(ab ~ PlantOrigin*SedOrigin + (1|TankID), data=aboveBelow, na.action = na.exclude)
ab.1way <- lmer(ab ~ PlantOrigin  + SedOrigin + (1|TankID), data=aboveBelow, na.action = na.exclude )
model.sel(ab.2way, ab.1way) 
#ab.1way is best fit

ab.full <- lmer(ab ~ PlantOrigin  + SedOrigin + (1|TankID), data=aboveBelow, na.action = na.exclude )
plot(ab.full)
qqnorm(resid(ab.full))
qqline(resid(ab.full))
shapiro.test(resid(ab.full)) #p=0.9384
summary(ab.full) 
Anova(ab.full)

#looking at contrasts for 2-way interaction
ab.emm <- emmeans(ab.full, ~ PlantOrigin * SedOrigin)
contrast(ab.emm, "consec", simple = "each", combine = FALSE, adjust = "mvt")



###########################
##BOXPLOT OF Avg SS Size##
##########################
avgSs <- meta %>% select (SampleID:TankID, avgSsSize, poptreat)


#summary stats
sum.avgSsSize <- summaryBy(avgSsSize ~ PlantOrigin + SedOrigin, data=avgSs, FUN=fun)

#plot
avgSs.plot <- ggplot(sum.avgSsSize, aes(x=SedOrigin, y=avgSsSize.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=avgSsSize.m-avgSsSize.se, ymax =avgSsSize.m+avgSsSize.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Average Side Shoot Size (g/no.shoots)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("BL", "MP")) + ylab("Average sideshoot size (g)")
avgSs.plot

#ANALYSIS
hist(avgSs$avgSsSize) #right tail looking; one outlier
avgSs$sqrtSs <- sqrt(avgSs$avgSsSize)
hist(avgSs$sqrtSs) #more normal with an outlier
avgSs.2way  <- lmer(sqrtSs ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=avgSs, na.action = na.exclude)
avgSs.1way <- lmer(sqrtSs ~  PlantOrigin + SedOrigin + (1|TankID), data=avgSs, na.action = na.exclude)
#see if 2way can drop
anova(avgSs.2way, avgSs.1way) #drop 
avgSs.full <-lmer( sqrtSs ~ PlantOrigin + SedOrigin + (1|TankID), data=avgSs, na.action = na.exclude)
plot(avgSs.2way)
qqnorm(resid(avgSs.2way))
qqline(resid(avgSs.2way)) #has 2 outliers
summary(avgSs.2way) 
Anova(avgSs.2way, type = 2) #sed origin

###########################
##BOXPLOT OF SS MASS##
##########################
ss <- meta %>% select (SampleID:TankID, side.g)
#summary stats
sum.ss <- summaryBy(side.g ~ PlantOrigin + SedOrigin, data=ss, FUN=fun)

#plot
ss.plot <- ggplot(sum.ss, aes(x=SedOrigin, y=side.g.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=side.g.m-side.g.se, ymax =side.g.m+side.g.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Side Shoot Mass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=30))+
  theme(legend.position = "none" , plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
ss.plot

#ANALYSIS
#not rooted
ss.2way  <- lmer(side.g ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=ss, na.action = na.exclude)
ss.1way <- lmer(side.g ~ PlantOrigin + SedOrigin + (1|TankID), data=ss, na.action = na.exclude)
#see if 2way can drop
anova(ss.2way, ss.1way) #drop 
ss.full <- lmer(side.g ~ PlantOrigin + SedOrigin + (1|TankID), data=ss, na.action = na.exclude)
#not rooted
plot(ss.full)
qqnorm(resid(ss.full))
qqline(resid(ss.full)) #normal
summary(ss.full) 
Anova( ss.full) 


###########################
##BOXPLOT OF SS COUNT##
##########################
ssCount <- meta %>% select (SampleID:TankID, NoSS)
#summary stats
sum.ssCount <- summaryBy(NoSS ~ PlantOrigin + SedOrigin, data=ssCount, FUN=fun)

#plot
ssCount.plot <- ggplot(sum.ssCount, aes(x=SedOrigin, y=NoSS.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=NoSS.m-NoSS.se, ymax =NoSS.m+NoSS.se), width=.1, size =1) +
  theme_classic() +
  ylab(expression(paste( "Number of Side Shoots"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
ssCount.plot

#ANALYSIS
#not rooted
ssCount.2way  <- glm(NoSS ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin, data=ssCount, na.action = na.exclude, family= "poisson")
ssCount.1way <- glm(NoSS ~ PlantOrigin + SedOrigin, data=ssCount, na.action = na.exclude, family= "poisson")
#see if can drop 2way
anova(ssCount.2way, ssCount.1way, test = "Chisq") #drop
ssCount.full <- glm(NoSS ~ PlantOrigin + SedOrigin, data=ssCount, na.action = na.exclude, family= "poisson")
summary(ssCount.full)
Anova(ssCount.full, type=2) 
#no effects
###########################
##BOXPLOT OF TERMINAL BIOMASSe##
##########################
term <- meta %>% select (SampleID:TankID, term.g, poptreat)
#summary stats
sum.term <- summaryBy(term.g ~ PlantOrigin + SedOrigin, data=term, FUN=fun)

#plot
term.plot <- ggplot(sum.term, aes(x=SedOrigin, y=term.g.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=term.g.m-term.g.se, ymax =term.g.m+term.g.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Terminal Shoot (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = "none", plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_x_discrete(NULL, labels = c("BL", "MP"))
term.plot

#ANALYSIS
str(term)
hist(term$term.g)
term.2way  <- lmer(term.g ~ PlantOrigin*SedOrigin +PlantOrigin +SedOrigin + (1|TankID), data=term, na.action = na.exclude)
term.1way <- lmer(term.g ~ PlantOrigin +SedOrigin+ (1|TankID) , data=term, na.action = na.exclude)
#determine if 2way can drop
anova(term.2way, term.1way) #drop
term.full <- lmer(term.g ~ PlantOrigin +SedOrigin + (1|TankID), data=term, na.action = na.exclude)
plot(term.full)
qqnorm(resid(term.full))
qqline(resid(term.full))
summary(term.full) 
Anova(term.full, type = 2) #plant origin

###########################
##BOXPLOT OF ROOT BUNDLES COUNT##
##########################
rootBund <- meta %>% select (SampleID:TankID, No_RootBundles)
#summary stats
sum.rootBund <- summaryBy(No_RootBundles ~ PlantOrigin + SedOrigin, data=rootBund, FUN=fun)

#plot
rootBund.plot <- ggplot(sum.rootBund, aes(x=SedOrigin, y=No_RootBundles.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=No_RootBundles.m-No_RootBundles.se, ymax =No_RootBundles.m+No_RootBundles.se), width=.1,size=1) +
  theme_classic() +
  ylab(expression(paste( "No. Root Bundles"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = "none", plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
rootBund.plot

#ANALYSIS
#not rooted
rootBund.2way  <- glmer(No_RootBundles ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID), data=rootBund, na.action = na.exclude, family= "poisson")
rootBund.1way <- glmer(No_RootBundles ~ PlantOrigin + SedOrigin+ (1|TankID), data=rootBund, na.action = na.exclude, family= "poisson")
#see if 2way can drop
anova(rootBund.2way, rootBund.1way) #drop
rootBund.full <- glmer(No_RootBundles ~ PlantOrigin + SedOrigin+ (1|TankID), data=rootBund, na.action = na.exclude, family= "poisson")
summary(rootBund.full) 
Anova(rootBund.full) #plant
#############
## TOTAL GROWTH RATE
#############
growthRate.data <- meta %>% select (SampleID:TankID, growthRate, poptreat)

#summary stats
sum.growth <- summaryBy(growthRate ~ PlantOrigin + SedOrigin, data=growthRate.data, FUN=fun)

#plot
growth.plot <- ggplot(sum.growth, aes(x=SedOrigin, y=growthRate.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=growthRate.m-growthRate.se, ymax =growthRate.m+growthRate.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Growth Rate (cm/day)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) +
  ylim(35, 60)
growth.plot

#analysis
hist(growthRate.data$growthRate)
growth.2way <- lmer(growthRate ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID) , data=growthRate.data, na.action = na.exclude)
growth.1way <- lmer(growthRate ~ PlantOrigin + SedOrigin+ (1|TankID), data=growthRate.data, na.action = na.exclude)
#see if 2way can drop
anova(growth.2way, growth.1way) #drop

growth.full <- lmer(growthRate ~ PlantOrigin + SedOrigin+ (1|TankID), data=growthRate.data, na.action = na.exclude)
plot(growth.full)
qqnorm(resid(growth.full))
qqline(resid(growth.full))
shapiro.test(resid(growth.full)) #passes
summary(growth.full) 
Anova(growth.full) #no differences

############
##TOTAL GROWTH
#############
totalGrowth.data <- meta %>% select (SampleID:TankID, totalGrowth, poptreat)

#summary stats
sum.Totalgrowth <- summaryBy(totalGrowth ~ PlantOrigin + SedOrigin, data=totalGrowth.data, FUN=fun)

#plot
Totalgrowth.plot <- ggplot(sum.Totalgrowth, aes(x=SedOrigin, y=totalGrowth.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin)) +
  geom_errorbar(aes(ymin=totalGrowth.m-totalGrowth.se, ymax =totalGrowth.m+totalGrowth.se), width=.1) +
  theme_classic() +
  ylab(expression(paste( "Total Growth (cm)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("darkolivegreen", "brown2"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
Totalgrowth.plot

#analysis
hist(totalGrowth.data$totalGrowth)
Totalgrowth.2way <- lm(totalGrowth ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin  , data=totalGrowth.data, na.action = na.exclude)
Totalgrowth.1way <- lm(totalGrowth ~ PlantOrigin + SedOrigin , data=totalGrowth.data, na.action = na.exclude)
#see if 2way can drop
anova(Totalgrowth.2way, Totalgrowth.1way) #drop

Totalgrowth.full <-  lm(totalGrowth ~ PlantOrigin + SedOrigin , data=totalGrowth.data, na.action = na.exclude)
plot(Totalgrowth.full)
qqnorm(resid(Totalgrowth.full))
qqline(resid(Totalgrowth.full))
shapiro.test(resid(Totalgrowth.full)) #passes
summary(Totalgrowth.full) 
Anova(Totalgrowth.full, type=2) #no differences


#####################
##BOXPLOT OF LONGEST ROOT_MM##
#####################
longRt <- meta %>% select(SampleID:TankID, poptreat, LongestRoot_mm)
#summary stats
sum.longestRoot <- summaryBy(LongestRoot_mm ~ PlantOrigin + SedOrigin, data=longRt, FUN=fun)

#plot
longestRoot.plot <- ggplot(sum.longestRoot, aes(x=SedOrigin, y=LongestRoot_mm.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=LongestRoot_mm.m-LongestRoot_mm.se, ymax =LongestRoot_mm.m+LongestRoot_mm.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Longest Root (mm)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=30))+
  theme(legend.position = "none", plot.margin = margin(1,1,1.5,1.2, "cm") )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
longestRoot.plot

#ANALYSIS
hist(longRt$LongestRoot_mm)
longRt.2way  <- lmer(LongestRoot_mm ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+  (1|TankID), data=longRt, na.action = na.exclude)
longRt.1way <- lmer(LongestRoot_mm ~ PlantOrigin  + SedOrigin + (1|TankID), data=longRt, na.action = na.exclude )
#seeding if 2way can be dropped
anova(longRt.2way, longRt.1way) #drop
longRt.full <- lmer(LongestRoot_mm ~ PlantOrigin  + SedOrigin + (1|TankID), data=longRt, na.action = na.exclude )
plot(longRt.full)
qqnorm(resid(longRt.full))
qqline(resid(longRt.full))
shapiro.test(resid(longRt.full)) #passes
summary(longRt.full) 
Anova(longRt.full) 

###################
##BOXPLOT OF ROOT BIOMASS##
###################
rtBiomass <- meta %>% select(SampleID:TankID, poptreat, root.g)
#summary stats
sum.roots <- summaryBy(root.g ~ PlantOrigin + SedOrigin, data=rtBiomass, FUN=fun)

#plot
roots.plot <- ggplot(sum.roots, aes(x=SedOrigin, y=root.g.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin)) +
  geom_errorbar(aes(ymin=root.g.m-root.g.se, ymax =root.g.m+root.g.se), width=.1) +
  theme_classic() +
  ylab(expression(paste( "Roots Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("darkolivegreen", "brown2"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
roots.plot


#ANALYSIS
hist(rtBiomass$root.g)
rtBiomass$sqRoot <-(rtBiomass$root.g)^(1/4)
hist(rtBiomass$sqRoot)
rtBiomass.2way  <- lmer(sqRoot ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID), data=rtBiomass, na.action = na.exclude)
rtBiomass.1way <- lmer(sqRoot ~ PlantOrigin  + SedOrigin + (1|TankID), data=rtBiomass, na.action = na.exclude )
#see if 2-way can be dropped
anova(rtBiomass.2way, rtBiomass.1way)  #drop

rootBio.full <-lmer(sqRoot ~ PlantOrigin  + SedOrigin + (1|TankID), data=rtBiomass, na.action = na.exclude )
#dropped random effect bc there was a singular boundary fit error
plot(rootBio.full)
qqnorm(resid(rootBio.full))
qqline(resid(rootBio.full))
shapiro.test(resid(rootBio.full)) #p<0.05
summary(rootBio.full) 
Anova(rootBio.full, type = 2) #plant

###################
##BOXPLOT OF BELOW BIOMASS##
###################
belBiomass <- meta %>% select(SampleID:TankID, poptreat, sumBel)
#summary stats
sum.below <- summaryBy(sumBel ~ PlantOrigin + SedOrigin, data=belBiomass, FUN=fun)

#plot
below.plot <- ggplot(sum.below, aes(x=SedOrigin, y=sumBel.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=sumBel.m-sumBel.se, ymax =sumBel.m+sumBel.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Total Below Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = "none", plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
below.plot

#ANALYSIS
belBiomass$sqrtBel <- sqrt(belBiomass$sumBel)
hist(belBiomass$sqrtBel)
belBiomass.2way  <- lmer(sqrtBel ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID), data=belBiomass, na.action = na.exclude)
belBiomass.1way <- lmer(sqrtBel ~ PlantOrigin  + SedOrigin + (1|TankID), data=belBiomass, na.action = na.exclude )
#see if 2way can be dropped
anova(belBiomass.2way, belBiomass.1way) #drop

belBiomass.full <-lmer(sqrtBel ~ PlantOrigin  + SedOrigin + (1|TankID), data=belBiomass, na.action = na.exclude )
plot(belBiomass.full)
qqnorm(resid(belBiomass.full))
qqline(resid(belBiomass.full))
summary(belBiomass.full) 
Anova(belBiomass.full, type=2) #plant

###################
##BOXPLOT OF RHIZOME BIOMASS##
###################
rhizBiomass <- meta %>% select(SampleID:TankID, poptreat, rhiz.g)
#summary stats
sum.rhiz <- summaryBy(rhiz.g ~ PlantOrigin + SedOrigin, data=rhizBiomass, FUN=fun)

#plot
rhiz.plot <- ggplot(sum.rhiz, aes(x=SedOrigin, y=rhiz.g.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin)) +
  geom_errorbar(aes(ymin=rhiz.g.m-rhiz.g.se, ymax =rhiz.g.m+rhiz.g.se), width=.1) +
  theme_classic() +
  ylab(expression(paste( "Rhizome Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("darkolivegreen", "brown2"))+
  theme(text=element_text(size=20))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
rhiz.plot

#ANALYSIS

rhizBiomass.2way  <- lmer(rhiz.g ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID), data=rhizBiomass, na.action = na.exclude)
rhizBiomass.1way <- lmer(rhiz.g ~ PlantOrigin  + SedOrigin + (1|TankID), data=rhizBiomass, na.action = na.exclude )
anova(rhizBiomass.2way, rhizBiomass.1way) #drop

rhizBiomass.full <- lmer(rhiz.g ~ PlantOrigin  + SedOrigin + (1|TankID), data=rhizBiomass, na.action = na.exclude )
plot(rhizBiomass.full)
qqnorm(resid(rhizBiomass.full))
qqline(resid(rhizBiomass.full))
shapiro.test(resid(rhizBiomass.full)) #passes
summary(rhizBiomass.full) 
Anova(rhizBiomass.full, type=2) #plant

###################
##BOXPLOT OF TOTAL BIOMASS##
###################
totalBiomass <- meta %>% select(SampleID:TankID, poptreat, totalBio)
#summary stats
sum.totBio <- summaryBy(totalBio ~ PlantOrigin + SedOrigin, data=totalBiomass, FUN=fun)

#plot
totBio.plot <- ggplot(sum.totBio, aes(x=SedOrigin, y=totalBio.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=2)+ geom_line(aes(linetype=PlantOrigin),size=1) +
  geom_errorbar(aes(ymin=totalBio.m-totalBio.se, ymax =totalBio.m+totalBio.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Total Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=30))+
  theme(legend.position = "none", plot.margin = margin(1,1,1.5,1.2, "cm")  )+
  scale_x_discrete(NULL, labels = c("BL", "MP"))
totBio.plot

#ANALYSIS
hist(totalBiomass$totalBio)
totalBiomass$sqrtTotBio <- sqrt(totalBiomass$totalBio)
totBiomass.2way  <- lmer(sqTotalBio ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin+ (1|TankID), data=totalBiomass, na.action = na.exclude)
totBiomass.1way <- lmer(sqTotalBio ~ PlantOrigin  + SedOrigin + (1|TankID), data=totalBiomass, na.action = na.exclude )
#see if 2way can drop
anova(totBiomass.2way, totBiomass.1way) #drop
totBiomass.full <-  lmer(sqrtTotBio ~ PlantOrigin  + SedOrigin + (1|TankID), data=totalBiomass, na.action = na.exclude )
plot(totBiomass.full)
qqnorm(resid(totBiomass.full))
qqline(resid(totBiomass.full))
shapiro.test(resid(totBiomass.full))
summary(totBiomass.full) 
Anova(totBiomass.full) #3 outliers; plant origin

###############
##TOTAL ABOVE
###############
totalAbove <- meta %>% select (SampleID:TankID, sumAbove)

#summary stats
sum.above <- summaryBy(sumAbove ~ PlantOrigin + SedOrigin, data=totalAbove, FUN=fun)

#plot
above.plot <- ggplot(sum.above, aes(x=SedOrigin, y=sumAbove.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin),size=1) +
  geom_errorbar(aes(ymin=sumAbove.m-sumAbove.se, ymax =sumAbove.m+sumAbove.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Total Above Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=20))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
above.plot

#analysis
hist(totalAbove$sumAbove)
totalAbove$sqrtAbove <- sqrt(totalAbove$sumAbove)
hist(totalAbove$sqrtAbove)
above.2way  <- lmer(sqrtAbove ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=totalAbove, na.action = na.exclude)
above.1way <- lmer(sqrtAbove ~  PlantOrigin + SedOrigin + (1|TankID), data=totalAbove, na.action = na.exclude)
#see if can drop 2way
anova(above.2way, above.1way) #drop
above.full <- lmer(sqrtAbove ~  PlantOrigin + SedOrigin + (1|TankID), data=totalAbove, na.action = na.exclude)
#not rooted
plot(above.full)
qqnorm(resid(above.full))
qqline(resid(above.full))
summary(above.full) 
Anova(above.full, type=2) #three outlier residuals; sed and plant sig

################
##RHIZOME LENGTH
################
rhizLength <- meta %>% select (SampleID:TankID, Rhiz_Length_mm)

#summary stats
sum.rhizLength <- summaryBy(Rhiz_Length_mm ~ PlantOrigin + SedOrigin, data=rhizLength, FUN=fun)

#plot
rhizLength.plot <- ggplot(sum.rhizLength, aes(x=SedOrigin, y=Rhiz_Length_mm.m, group=PlantOrigin, color=PlantOrigin)) +
  geom_point(size=1)+ geom_line(aes(linetype=PlantOrigin), size=1) +
  geom_errorbar(aes(ymin=Rhiz_Length_mm.m-Rhiz_Length_mm.se, ymax =Rhiz_Length_mm.m+Rhiz_Length_mm.se), width=.1, size=1) +
  theme_classic() +
  ylab(expression(paste( "Rhizome Length (mm)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("#0C7BDC", "#FFC20A"))+
  theme(text=element_text(size=30))+
  theme(legend.position = )+
  scale_x_discrete(NULL, labels = c("BL", "MP")) 
rhizLength.plot

#analysis
rhizLength.2way  <- lmer(Rhiz_Length_mm ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=rhizLength, na.action = na.exclude)
rhizLength.1way <- lmer(Rhiz_Length_mm ~  PlantOrigin + SedOrigin+ (1|TankID), data=rhizLength, na.action = na.exclude)
#see if can drop 2way
anova(rhizLength.2way, rhizLength.1way) #keep
rhizLength.full <- lmer(Rhiz_Length_mm ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=rhizLength, na.action = na.exclude)
plot(rhizLength.full)
qqnorm(resid(rhizLength.full))
qqline(resid(rhizLength.full))
shapiro.test(resid(rhizLength.full)) #passes
summary(rhizLength.full) 
Anova(rhizLength.full, type=2)  #plant, sed, and interaction sig
Anova(lmer(Rhiz_Length_mm ~ PlantOrigin*SedOrigin + PlantOrigin + SedOrigin + (1|TankID), data=rhizLength, na.action = na.exclude, 
           contrasts=list(PlantOrigin=contr.sum, SedOrigin=contr.sum)), type=3)
