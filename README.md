# **Plasticity and adaptation of eelgrass and its microbiome in response to sediment conditions**

**Karolina Zabinski<sup>1,2</sup>, Claire Murphy<sup>1,2</sup>, Katherine DuBois<sup>3</sup>, John J. Stachowicz<sup>1</sup>**    

<sup>1</sup> Department of Evolution and Ecology, University of California Davis, Davis, California, USA  
<sup>2</sup> Bodega Marine Laboratory, University of California Davis, Bodega Bay, California, USA  
<sup>3</sup> School of Ocean Sciences, Bangor University, Isle of Anglesey LL59 5AB, UK  

# Abstract from Manuscript:  
1.	Interactions between local adaptation and plasticity shape population differentiation in response to stressors throughout a species’ range. In addition, there is gaining interest in understanding how microbes may play a role in mitigating organismal performance by extending phenotypes. 
2.	Here, we use a common garden mesocosm experiment to assess eelgrass (Zostera marina) plastic and fixed adaptive responses to sediment grain size, considering not only variation in plant traits but also variation in the plant-associated microbiome that may contribute to eelgrass growth and reproductive success. 
3.	We found consistent fixed differences in plant traits from the two populations, as plants from the clay-dominated, fine-grained site maintained greater above and below ground biomass regardless of sediment type. 
4.	We also found plasticity in both populations. Clay-dominated sediment increased root length and clonal side shoot size compared to the larger grained, sand-dominated sediment. 
5.	Plant by sediment origin interactions indicated some measure of home site advantage with respect to sediment conditions. Specifically, plants from both populations reduced porewater sulfide to low levels in fine-grained sediment, but only plants native to sand-dominated sediment decreased porewater sulfide in the coarse-grained treatment.
6.	This home site advantage may be mediated by microbiome differences between populations as plants native to clay-dominated sediment had approximately half the relative abundance of a dominant sulfur oxidizer, Sulfurimonadaceae, compared to plants native to sand-dominated sediment when grown in sediment from the coarse-grained site. 
7.	Synthesis: These results support that sediment partially mediates home site advantage in eelgrass populations and suggest differential population responses may be mediated by the associated microbial community.

# Notes on Scripts:  
1. General note: Working directories are set at the start of each .RMD or .R. All downstream paths are relative to that set working directory.
2. GSE-preprocessing.r : running DADA2 pipeline and preparing phyloseq object for downstream analyses
3. plant-responses-&-analyses.rmd : analyses and figures corresponding to 1-3, S1-2
4. community-&-distances-analyses.rmd : analyses and figures corresponding to 4A-C
5. alpha-diversity-beta-dispersion-analyses-&-figures.rmd : analyses and figures corresponding to Fig. 2B and S5A-D
6. DAA-figures-&-analyses.rmd : analyses and figures corresponding to Fig. 5 and S4-6
7. GSE-temperature.rmd : analyses of logger data to caluclate temperature means

# Notes on Program Version:  
1. R: 4.4.1

# Notes on Package Versions:  
1. dada2 (1.26 or later) for preprocessing amplicon data
2. tdyr (1.3.1), dplyr (1.1.4), tidyverse (2.0.0), reshape2 (1.4.4) for easier data manipulation
3. doBy (4.6.22) for groupwise calculations
4. ggplot2(3.5.1), gridExtra (2.3), ggpubr (0.6.0), cowplot (1.1.3) for figure creation and manipulation
5. formatR (1.14), usethis(3.0.0) for github package integration and rmd formating
6. lme4 (1.1.35.5), car(3.1.2), emmeans(1.10.3) for linear mixed model analyses
7. phyloseq (1.48.0), speedyseq (0.5.3.9021) for manipulation of microbiome (16s rRNA amplicon) dataset
8. corncob (0.4.1) for differential abundance tests of microbiome data
9. lubridate (1.9.3) manipulate time series data/logger data

