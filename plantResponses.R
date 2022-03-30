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
