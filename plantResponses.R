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

meta$productivity <- as.numeric(meta$totalGrowth*meta$ShootWidth_mm/meta$DaysGrowing) #all NAs remained NAs

#average ss side
meta <- meta %>% mutate(avgSsSize = side.g/(ShootCount-1))
meta$avgSsSize[meta$avgSsSize == 0] <- NA
meta$avgSsSize[meta$avgSsSize == Inf] <- NA
meta$avgSsSize[meta$avgSsSize == -Inf] <- NA
meta$avgSsSize[meta$avgSsSize == NaN] <- NA



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
##PROPORTION ROOT HAIR AREA##
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
##BOXPLOT OF PRODUCTIVITY##
##########################
prod <- meta %>% select (SampleID:TankID, productivity, poptreat)

prod.plot <- ggplot(prod, aes(x=PlantOrigin, y=productivity, fill=SedOrigin)) +
  geom_boxplot()
prod.plot

###########################
##BOXPLOT OF AB##
##########################
aboveBelow <- meta %>% select (SampleID:TankID, ab, poptreat)
#setting seed so that code is run from same random path
set.seed(5)
ab.plot <- ggplot(aboveBelow, aes(x=PlantOrigin, y=ab, fill=SedOrigin)) +
  geom_boxplot() + 
  geom_jitter()
ab.plot

###########################
##BOXPLOT OF Avg SS Size##
##########################
avgSs <- meta %>% select (SampleID:TankID, avgSsSize, poptreat)
#setting seed so that code is run from same random path
set.seed(5)
avgSs.plot <- ggplot(avgSs, aes(x=PlantOrigin, y=avgSsSize, fill=SedOrigin)) +
  geom_boxplot() + 
  geom_jitter()
avgSs.plot


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
