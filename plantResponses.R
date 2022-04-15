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
cols <- names(meta)[5:13]
meta[cols] <- lapply(meta[cols], as.numeric)
str(meta)
#convert DaysGrowing:No_RootBundles to numeric
cols2 <- names(meta)[16:19]
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

#########################
##REGRESSING LONGEST ROOT TO CHEMISTRY##
##########################
#Selecting f2size and severity
rlSul <- meta %>% select (SampleID:TankID, LongestRoot_mm, sulfide.um, poptreat)
hist(rlSul$LongestRoot_mm) #flat nor curve
hist(rlSul$sulfide.um) #very not normal; has those outliers

#with outliers
rlSul.regclass <- ggplot(rlSul, aes(x= LongestRoot_mm, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Sulfide (um) vs Longest Root (mm)"))) +
  geom_smooth(method=lm) 
rlSul.regclass 

##WITHOUT OUTLIERS
#finding outlier
boxplot.stats(rlSul$sulfide.um)$out
outSul.tot <- boxplot.stats(rlSul$sulfide.um)$out
out_indSul.tot <- which(rlSul$sulfide.um %in% c(outSul.tot))
out_indSul.tot
#outliers are: 2-10 (BLAK-S), 3-1 (MILL-O), 3-2 (MILL-O), 3-11 (MILL-O), 3-12 (MILL-O)

#removing outliers
rlSul.out <- rlSul
rlSul.out <- rlSul.out[rlSul.out$SampleID != "2-10",]
rlSul.out <- rlSul.out[rlSul.out$SampleID != "3-1",]
rlSul.out <- rlSul.out[rlSul.out$SampleID != "3-2",]
rlSul.out <- rlSul.out[rlSul.out$SampleID != "3-11",]
rlSul.out <- rlSul.out[rlSul.out$SampleID != "3-12",]
str(rlSul.out)
droplevels(rlSul.out)
hist(rlSul.out$sulfide.um)

#without outliers
rlSulOUT.regclass <- ggplot(rlSul.out, aes(x= LongestRoot_mm, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Sulfide (um) vs Longest Root (mm)"))) +
  geom_smooth(method=lm) 
rlSulOUT.regclass 

##ARRANGING TWO TGP PLOTS FOR A:B
grid.arrange(rlSul.regclass, rlSulOUT.regclass, ncol=2)

#########################
##REGRESSING AVG ROOT HAIR LENGTH TO CHEMISTRY##
##########################
#Selecting f2size and severity
rhSul <- meta %>% select (SampleID:TankID, rhLength, sulfide.um, poptreat)
hist(rhSul$rhLength) #right tailed
hist(rhSul$sulfide.um) #very not normal; has those outliers

#with outliers
rhSul.regclass <- ggplot(rhSul, aes(x= rhLength, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Sulfide (um) vs Average Root Hair Length (mm)"))) +
  geom_smooth(method=lm) 
rhSul.regclass 

#############################
##REG PROPORTION ROOT HAIR AREA TO SULF##
#############################
#Selecting f2size and severity
proprhSul <- meta %>% select (SampleID:TankID, propotion.root.hair.area.mm, sulfide.um, poptreat)
hist(proprhSul$propotion.root.hair.area.mm) #pretty norm
hist(proprhSul$sulfide.um) #very not normal; has those outliers

#with outliers
proprhSul.regclass <- ggplot(proprhSul, aes(x= propotion.root.hair.area.mm, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Sulfide (um) vs Proportion of Root Hair Area"))) +
  geom_smooth(method=lm) 
proprhSul.regclass 

#############################
##AB VS Sulfide##
#############################
#Selecting f2size and severity
abSul <- meta %>% select (SampleID:TankID, ab, sulfide.um, poptreat)
hist(abSul$ab) #pretty norm
hist(abSul$sulfide.um) #very not normal; has those outliers

#with outliers
abSul.regclass <- ggplot(abSul, aes(x= ab, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Sulfide (um) vs Above:Below"))) +
  geom_smooth(method=lm) 
abSul.regclass 


###########################
##Avg SS Size vs Below##
##########################
avgSSBel <- meta %>% select (SampleID:TankID, avgSsSize, rhiz.g, sumBel, poptreat)
hist(avgSSBel$avgSsSize) #pretty norm right outlier
hist(avgSSBel$rhiz.g) #broad right tail dist

avgSSBel.regclass <- ggplot(avgSSBel, aes(x= sumBel, y= avgSsSize, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[1], "Avg Ss Size (g) vs Below (g)"))) +
  geom_smooth(method=lm) 
avgSSBel.regclass

###########################
##productivity vs Avg SS Size##
##########################
prodAvgSS <- meta %>% select (SampleID:TankID, avgSsSize, productivity, poptreat)
hist(prodAvgSS$avgSsSize) #pretty norm right outlier
hist(prodAvgSS$productivity) #norm

prodAvgSS.regclass <- ggplot(prodAvgSS, aes(x= productivity, y= avgSsSize, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[1], "Avg SS Size (g) vs", F[0], "Productivity (cm^2)/day"))) +
  geom_smooth(method=lm) 
prodAvgSS.regclass

###########################
##longest root vs avg root hair length##
##########################
rootRh <- meta %>% select (SampleID:TankID, LongestRoot_mm, rhLength, poptreat)
hist(rootRh$LongestRoot_mm) #pretty norm 
hist(rootRh$rhLength) #right tailed

rootRh.regclass <- ggplot(rootRh, aes(x= LongestRoot_mm, y= rhLength, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "Longest root (mm) vs Average root hair lenght(cm)"))) +
  geom_smooth(method=lm) 
rootRh.regclass

###########################
##root vs proportion root hair area##
##########################
rhizR <- meta %>% select (SampleID:TankID, propotion.root.hair.area.mm, root.g, poptreat)
hist(rhizR$propotion.root.hair.area.mm) #pretty norm 
hist(rhizR$root.g) #right tailed

rhizR.regclass <- ggplot(rhizR, aes(x= root.g, y= propotion.root.hair.area.mm, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "root (g) vs propotion of root hair area"))) +
  geom_smooth(method=lm) 
rhizR.regclass

###########################
##root.g vs sulfide##
##########################
belSulf <- meta %>% select (SampleID:TankID, sumBel, sulfide.um, root.g, poptreat)
hist(belSulf$root.g) # right tailed
hist(belSulf$sulfide.um) #right tailed

rhizR.regclass <- ggplot(belSulf, aes(x= root.g, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "below (g) vs sulfide (um)"))) +
  geom_smooth(method=lm) 
rhizR.regclass

###########################
##prop root hair vs sulfide##
##########################
propRHSulf <- meta %>% select (SampleID:TankID, propotion.root.hair.area.mm, sulfide.um, poptreat)
hist(propRHSulf$propotion.root.hair.area.mm) #norm
hist(propRHSulf$sulfide.um) #right tailed

propRHSulf.regclass <- ggplot(propRHSulf, aes(x= propotion.root.hair.area.mm, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "proportion of root hairs vs sulfide (um)"))) +
  geom_smooth(method=lm) 
propRHSulf.regclass

#removing outliers
propRHSulf.out <- propRHSulf
propRHSulf.out  <- propRHSulf.out [propRHSulf.out $SampleID != "2-10",]
propRHSulf.out  <- propRHSulf.out [propRHSulf.out $SampleID != "3-1",]
propRHSulf.out  <- propRHSulf.out [propRHSulf.out $SampleID != "3-2",]
propRHSulf.out  <- propRHSulf.out [propRHSulf.out $SampleID != "3-11",]
propRHSulf.out  <- propRHSulf.out [propRHSulf.out $SampleID != "3-12",]
str(propRHSulf.out )
droplevels(propRHSulf.out )
hist(propRHSulf.out $sulfide.um)

#without outliers
propRHSulfOUT.regclass <- ggplot(propRHSulf.out, aes(x= propotion.root.hair.area.mm, y= sulfide.um, color= poptreat)) + geom_point(size=3) +
  scale_color_manual(values=c("orangered2","skyblue2", "deeppink","cornflowerblue")) +
  ylab(expression(paste(F[0], "proportion of root hairs vs sulfide (um)"))) +
  geom_smooth(method=lm) 
propRHSulfOUT.regclass 

##ARRANGING TWO TGP PLOTS FOR A:B
grid.arrange(propRHSulf.regclass, propRHSulfOUT.regclass, ncol=2)

###########################
##BOXPLOT OF PRODUCTIVITY##
##########################
prod <- meta %>% select (SampleID:TankID, productivity, poptreat)

set.seed(5)
prod.plot <- ggplot(prod, aes(x=PlantOrigin, y=productivity, group=PlantOrigin)) + 
  geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
prod.plot

#ANALYSIS
prod.2way  <- lmer(productivity ~ PlantOrigin*SedOrigin + (1|TankID), data=prod, na.action = na.exclude)
prod.1way <- lmer(productivity ~ PlantOrigin  + SedOrigin + (1|TankID), data=prod, na.action = na.exclude )
model.sel(prod.2way, prod.1way) 
#prod.1way is best fit

prod.full <- lmer(productivity ~ PlantOrigin  + SedOrigin + (1|TankID), data=prod, na.action = na.exclude )
plot(prod.full)
qqnorm(resid(prod.full))
qqline(resid(prod.full))
shapiro.test(resid(prod.full)) #p=0.1057
summary(prod.full) 
Anova(prod.full) #no differences in productivity p >0.1

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
#setting seed so that code is run from same random path
set.seed(5)

avgSs.plot <- ggplot(avgSs, aes(x=PlantOrigin, y=avgSsSize, group=PlantOrigin)) + 
  geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
avgSs.plot

#ANALYSIS
avgSs.2way  <- lmer(avgSsSize ~ PlantOrigin*SedOrigin + (1|TankID), data=avgSs, na.action = na.exclude)
avgSs.1way <- lmer(avgSsSize ~ PlantOrigin  + SedOrigin + (1|TankID), data=avgSs, na.action = na.exclude )
model.sel(avgSs.2way, avgSs.1way) 
#avgSs.1way is best fit

avgSs.full <- gls(avgSsSize ~ PlantOrigin  + SedOrigin, data=avgSs, na.action = na.exclude )
#dropped random effect bc there was a singular boundary fit error
plot(avgSs.full)
qqnorm(resid(avgSs.full))
qqline(resid(avgSs.full))
shapiro.test(resid(avgSs.full)) #p<0.05
summary(avgSs.full) 
Anova(avgSs.full) 
#both plant and sed origin have effect of p<0.05

#looking at contrasts for 2-way interaction
avgSs.emm <- emmeans(avgSs.full, ~ PlantOrigin * SedOrigin)
contrast(avgSs.emm, "consec", simple = "each", combine = FALSE, adjust = "mvt")

#############
##SULFIDE BOXPLOT PER TREATMENT
##############
sul <- meta %>% select(SampleID:TankID, poptreat, sulfide.um)
set.seed(5)

#with outliers
sul.plot <- ggplot(sul, aes(x=PlantOrigin, y=sulfide.um, group=PlantOrigin)) + 
  geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
sul.plot

#without outliers
#removing outliers
#removing outliers
sul.out <- sul
sul.out <- sul.out[sul.out$SampleID != "2-10",]
sul.out <- sul.out[sul.out$SampleID != "3-1",]
sul.out <- sul.out[sul.out$SampleID != "3-2",]
sul.out <- sul.out[sul.out$SampleID != "3-11",]
sul.out <- sul.out[sul.out$SampleID != "3-12",]
sulOut.plot <- ggplot(sul.out, aes(x=PlantOrigin, y=sulfide.um, group=PlantOrigin)) + 
  geom_boxplot(color="gray") + facet_grid(. ~ SedOrigin) + geom_jitter(aes(color=PlantOrigin)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
sulOut.plot

#ANALYSIS WITHOUT OUTLIERS
#ANALYSIS
sulOut.2way  <- lmer(sulfide.um ~ PlantOrigin*SedOrigin + (1|TankID), data=sul.out, na.action = na.exclude)
sulOut.1way <- lmer(sulfide.um ~ PlantOrigin  + SedOrigin + (1|TankID), data=sul.out, na.action = na.exclude )
model.sel(sulOut.2way, sulOut.1way) 
#sulOut.2way is best fit

sulOut.full <- lmer(sulfide.um ~ PlantOrigin*SedOrigin + (1|TankID), data=sul.out, na.action = na.exclude)
#dropped random effect bc there was a singular boundary fit error
plot(sulOut.full)
qqnorm(resid(sulOut.full))
qqline(resid(sulOut.full))
shapiro.test(resid(sulOut.full)) #p<0.05
summary(sulOut.full) 
Anova(sulOut.full) 
#plant and sed origin have effect and interaction of p<0.05

#looking at contrasts for 2-way interaction
sulOut.emm <- emmeans(sulOut.full, ~ PlantOrigin * SedOrigin)
contrast(sulOut.emm, "consec", simple = "each", combine = FALSE, adjust = "mvt")

#ANALYSIS WITH OUTLIERS
sul.2way  <- lmer(sulfide.um ~ PlantOrigin*SedOrigin + (1|TankID), data=sul, na.action = na.exclude)
sul.1way <- lmer(sulfide.um ~ PlantOrigin  + SedOrigin + (1|TankID), data=sul, na.action = na.exclude )
model.sel(sul.2way, sul.1way) 
#sul.2way is best fit

sul.full <- lmer(sulfide.um ~ PlantOrigin*SedOrigin + (1|TankID), data=sul, na.action = na.exclude)
#dropped random effect bc there was a singular boundary fit error
plot(sul.full)
qqnorm(resid(sul.full))
qqline(resid(sul.full))
shapiro.test(resid(sul.full)) #p<0.05
summary(sul.full) 
Anova(sul.full) 
#plant and sed origin have effect and interaction of p<0.05

#looking at contrasts for 2-way interaction
sulOut.emm <- emmeans(sulOut.full, ~ PlantOrigin * SedOrigin)
contrast(sulOut.emm, "consec", simple = "each", combine = FALSE, adjust = "mvt")

##########
##PCA ATTEMPTS
##############
##loading libraries needed for analysis
library(missMDA)
library(FactoMineR)
library(naniar)
library(VIM)
library(factoextra)
library(naniar)

############
###PREPARING DF FOR MANIPULATION
###########
#subset meta to include just vars of intrest
reduced.meta <- meta %>% select(ShootCount, ShootLength_mm, LongestRoot_mm:No_RootBundles, proportion.root.area.mm, 
                                propotion.root.hair.area.mm:productivity)

#groups
groups <- meta %>% select(PlantOrigin, SedOrigin, poptreat, SampleID)

#################
##PCA ANALYSIS WITH DATA
#################
#pull only variables
vars <- reduced.meta
#calculating how many NAs are in df after removing dead individuals
gg_miss_var(vars)
#estimate number of dimentions for PCA
nb= estim_ncpPCA(vars,ncp.max=5) 
#imputes missing values in dataset
res.comp <- imputePCA(vars, ncp = nb$ncp)
res.comp
#looking at the imputed dataset
res.comp$completeObs[1:3,]
#merginging imputed dataset with class
imp <- cbind.data.frame(res.comp$completeObs,groups)
str(imp)
head(imp, 3)
#compute PCA without class
vars.pca <- PCA(imp[,-20:-23], graph = FALSE)
vars.pca
#coloring PCA indiviudals by class
fviz_pca_ind(vars.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = imp$poptreat, # color by groups
             palette = c("orangered2","skyblue2", "deeppink","cornflowerblue"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
#biplot overlayed on top of pca
fviz_pca_biplot(vars.pca, 
                col.ind = imp$poptreat, palette = c("orangered2","skyblue2", "deeppink","cornflowerblue"), 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Groups") 

##########################
###SUBSETTING DATA INTO PLANT AND SED ORIGINS FOR PCA
##########################
mpPlant <- meta %>% filter((PlantOrigin == "MILL"))

blPlant <- meta %>% filter((PlantOrigin == "BLAK"))

mpSed <- meta %>% filter((SedOrigin == "MILL"))

blSed <- meta %>% filter((SedOrigin == "BLAK"))
###############################
##PCA OF JUST MP PLANTS##
##############################
#subset break.2021 to include just vars of intrest
reduced.mpPlant <- mpPlant %>% select(PlantOrigin:SedOrigin, ShootLength_mm, LongestRoot_mm:No_RootBundles,
                                      sumAbove:ab, totalGrowth:productivity)
#removing dead individuals
processed.mpPlant <- reduced.mpPlant %>% filter(if_any(ShootLength_mm:productivity))

#groups
groups.mpPlant <- processed.mpPlant %>% select(PlantOrigin, SedOrigin)

#pull only variables
vars.mpPlant <- processed.mpPlant %>% select(ShootLength_mm:productivity)

#calculating how many NAs are in df after removing dead individuals
gg_miss_var(vars.mpPlant)
#estimate number of dimentions for PCA
nb.mpPlant= estim_ncpPCA(vars.mpPlant,ncp.max=5) 
#imputes missing values in dataset
res.comp.mpPlant<- imputePCA(vars.mpPlant, ncp = nb.mpPlant$ncp)
res.comp.mpPlant
#looking at the imputed dataset
res.comp.mpPlant$completeObs[1:3,]
#merginging imputed dataset with class
imp.mpPlant <- cbind.data.frame(res.comp.mpPlant$completeObs,groups.mpPlant)
str(imp.mpPlant)
head(imp.mpPlant, 3)
#compute PCA without class
vars.pca.mpPlant<- PCA(imp.mpPlant[,-11:-12], graph = FALSE)
vars.pca.mpPlant
#coloring PCA indiviudals by class
fviz_pca_ind(vars.pca.mpPlant,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = imp.mpPlant$SedOrigin, # color by groups
             palette = c("skyblue2","orangered2", "cornflowerblue","deeppink", "darkblue","coral2", "deepskyblue4", "brown2"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

###############################
##PCA OF JUST BL PLANTS##
##############################
#subset break.2021 to include just vars of intrest
reduced.blPlant <- blPlant %>% select(PlantOrigin:SedOrigin, ShootLength_mm, LongestRoot_mm:No_RootBundles,
                                      sumAbove:ab, totalGrowth:productivity)
#removing dead individuals
processed.blPlant <- reduced.blPlant %>% filter(if_any(ShootLength_mm:productivity))

#groups
groups.blPlant <- processed.blPlant %>% select(PlantOrigin, SedOrigin)

#pull only variables
vars.blPlant <- processed.blPlant %>% select(ShootLength_mm:productivity)

#calculating how many NAs are in df after removing dead individuals
gg_miss_var(vars.blPlant)
#estimate number of dimentions for PCA
nb.blPlant= estim_ncpPCA(vars.blPlant,ncp.max=5) 
#imputes missing values in dataset
res.comp.blPlant<- imputePCA(vars.blPlant, ncp = nb.blPlant$ncp)
res.comp.blPlant
#looking at the imputed dataset
res.comp.blPlant$completeObs[1:3,]
#merginging imputed dataset with class
imp.blPlant <- cbind.data.frame(res.comp.blPlant$completeObs,groups.blPlant)
str(imp.blPlant)
head(imp.blPlant, 3)
#compute PCA without class
vars.pca.blPlant<- PCA(imp.blPlant[,-11:-12], graph = FALSE)
vars.pca.blPlant
#coloring PCA indiviudals by class
fviz_pca_ind(vars.pca.blPlant,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = imp.blPlant$SedOrigin, # color by groups
             palette = c("skyblue2","orangered2", "cornflowerblue","deeppink", "darkblue","coral2", "deepskyblue4", "brown2"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

###############################
##PCA OF JUST MP SED##
##############################
#subset break.2021 to include just vars of intrest
reduced.mpSed <- mpSed %>% select(PlantOrigin:SedOrigin, ShootLength_mm, LongestRoot_mm:No_RootBundles,
                                      sumAbove:ab, totalGrowth:productivity)
#removing dead individuals
processed.mpSed <- reduced.mpSed %>% filter(if_any(ShootLength_mm:productivity))

#groups
groups.mpSed <- processed.mpSed %>% select(PlantOrigin, SedOrigin)

#pull only variables
vars.mpSed <- processed.mpSed %>% select(ShootLength_mm:productivity)

#calculating how many NAs are in df after removing dead individuals
gg_miss_var(vars.mpSed)
#estimate number of dimentions for PCA
nb.mpSed= estim_ncpPCA(vars.mpSed,ncp.max=5) 
#imputes missing values in dataset
res.comp.mpSed<- imputePCA(vars.mpSed, ncp = nb.mpSed$ncp)
res.comp.mpSed
#looking at the imputed dataset
res.comp.mpSed$completeObs[1:3,]
#merginging imputed dataset with class
imp.mpSed <- cbind.data.frame(res.comp.mpSed$completeObs,groups.mpSed)
str(imp.mpSed)
head(imp.mpSed, 3)
#compute PCA without class
vars.pca.mpSed<- PCA(imp.mpSed[,-11:-12], graph = FALSE)
vars.pca.mpSed
#coloring PCA indiviudals by class
fviz_pca_ind(vars.pca.mpSed,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = imp.mpSed$PlantOrigin, # color by groups
             palette = c("skyblue2","orangered2", "cornflowerblue","deeppink", "darkblue","coral2", "deepskyblue4", "brown2"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

###############################
##PCA OF JUST BL SED##
##############################
#subset break.2021 to include just vars of intrest
reduced.blSed <- blSed %>% select(PlantOrigin:SedOrigin, ShootLength_mm, LongestRoot_mm:No_RootBundles,
                                      sumAbove:ab, totalGrowth:productivity)
#removing dead individuals
processed.blSed  <- reduced.blSed  %>% filter(if_any(ShootLength_mm:productivity))

#groups
groups.blSed  <- processed.blSed  %>% select(PlantOrigin, SedOrigin)

#pull only variables
vars.blSed  <- processed.blSed  %>% select(ShootLength_mm:productivity)

#calculating how many NAs are in df after removing dead individuals
gg_miss_var(vars.blSed )
#estimate number of dimentions for PCA
nb.blSed = estim_ncpPCA(vars.blSed ,ncp.max=5) 
#imputes missing values in dataset
res.comp.blSed <- imputePCA(vars.blSed , ncp = nb.blSed $ncp)
res.comp.blSed 
#looking at the imputed dataset
res.comp.blSed$completeObs[1:3,]
#merginging imputed dataset with class
imp.blSed  <- cbind.data.frame(res.comp.blSed $completeObs,groups.blSed )
str(imp.blSed )
head(imp.blSed , 3)
#compute PCA without class
vars.pca.blSed<- PCA(imp.blSed [,-11:-12], graph = FALSE)
vars.pca.blSed 
#coloring PCA indiviudals by class
fviz_pca_ind(vars.pca.blSed ,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = imp.blSed $PlantOrigin, # color by groups
             palette = c("skyblue2","orangered2", "cornflowerblue","deeppink", "darkblue","coral2", "deepskyblue4", "brown2"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)




