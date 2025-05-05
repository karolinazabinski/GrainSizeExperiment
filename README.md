# **Plasticity and adaptation of northern California eelgrass in response to sediment conditions**

**Karolina L. Zabinski<sup>1,2</sup>, Claire E. Murphy<sup>1,2</sup>, Katherine DuBois<sup>3</sup>, John J. Stachowicz<sup>1</sup>**    

<sup>1</sup> Department of Evolution and Ecology, University of California Davis, Davis, California, USA  
<sup>2</sup> Bodega Marine Laboratory, University of California Davis, Bodega Bay, California, USA  
<sup>3</sup> School of Ocean Sciences, Bangor University, Isle of Anglesey LL59 5AB, UK  

# Abstract from Manuscript:  
Considerable research describes the interactions between seagrasses and their sedimentary environment, but there is little information on how populations differ in their innate versus plastic responses to these differences. Here, we test whether sediment contributes to eelgrass population differentiation and the nature of plastic responses to different sediment environments. We do this via a 15-week, fully crossed common garden experiment with two populations and their native sediment types. Plants from the warmer-temperature, clay-dominated site (90% silt + clay, 10% sand) consistently maintained greater biomass than plants from the cooler, sand-dominated site (60% sand, 40% silt + clay). Plants from both populations were highly plastic for root length and clonal shoot size, with both increasing when planted in clay-dominated compared to sand-dominated sediment. Plants from the clay-dominated site grew longer rhizomes in foreign sediment while plants from the sand-dominated site had no change in this plant trait, indicating some measure of home site advantage with respect to sediment conditions. Porewater sulfide also exhibited this pattern where concentrations were very low in clay-dominated sediment for all plants, but in the sand-dominated treatment, only plants native to sand-dominated sediment maintained porewater sulfide concentrations below toxic levels. These patterns may be mediated by microbiome differences between populations as roots from plants native to clay-dominated sediment had more fixed microbiomes between treatments compared to plants native to sand-dominated sediment. These results support that sediment type partially mediates home site advantage in eelgrass populations and suggest differential population responses may be mediated by the associated microbiome.

Keywords: seagrass, root microbiomes, sediment, plasticity, local adaptation


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

