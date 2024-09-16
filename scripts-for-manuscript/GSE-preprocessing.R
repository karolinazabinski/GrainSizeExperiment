###################
##DESCRIPTION#####
##################
##1. download rawfastq.gz files from box folder in this repo
##2. Trim sequences at quality score of around 30 and save file 
##3. Filter sequences to get rid of low quality scores and save file
##4. Dereplicate sequences by collapsing all reads with the same sequences and save file
##5. Merge R1 & R2 reads to remove pairs that don't match and save files
##6. Remove chimeras and save file
##7. Create sequence table as final output and save file

############
##Logistics##
############
##installing packages needed for analysis
if (!requireNamespace("BiocManager", quietly = TRUE))#bc this package needs an older version of R, call the compatible version
  install.packages("BiocManager")
BiocManager::install(version='3.18') # switching from R version 4.3.1 to 3.18 w/in this session
BiocManager::install("dada2", version = "3.18") 
library("dada2"); packageVersion("dada2") #check to see if downloaded
  # version 1.30.0
#if (!requireNamespace("BiocManager", quietly = TRUE)) #bc this package needs an older version of R, call the compatible version
#  + +     install.packages("BiocManager")
#BiocManager::install("phyloseq")
library("phyloseq"); packageVersion("phyloseq") #check to see if downloaded
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("vegan"); packageVersion("vegan")
library("grid"); packageVersion("grid")
library("knitr"); packageVersion("knitr")
#if (!requireNamespace("BiocManager", quietly = TRUE)) #bc this package needs an older version of R, call the compatible version
# + +     install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2"); packageVersion("DESeq2") #check to see if downloaded
#if (!requireNamespace("BiocManager", quietly = TRUE)) #bc this package needs an older version of R, call the compatible version
 #+ +     install.packages("BiocManager")
BiocManager::install("decontam")
library("decontam"); packageVersion("decontam") 
  # version 1.22.0
library(usethis) 
library(dplyr)

#############
##Setting WDs####
#############
rm(list=ls())
##setting main wd throughout session
setwd("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment")

##############
##Manipulating raw sequences via fastq files##
###############
##defining path to raw sequences##
path.raw <- "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/raw-FASTQ"
#checking to see if files are present
head(list.files(path.raw))

##reading names of fastq.gz files and creating list of matched forward and reverse reads##
##Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz
##creating list of forward reads.
fnFs <- sort(list.files(path.raw, pattern="_R1.fastq.gz", full.names = TRUE))
head(fnFs)
tail(fnFs)
#creating list of reverse reads. Based on Guillaumes's demultiplexing code, files ending in 1 are the Rvs reads
fnRs <- sort(list.files(path.raw, pattern="_R2.fastq.gz", full.names = TRUE))
head(fnRs)
tail(fnRs)

# Extract sample names, where filenames have format: SAMPLENAME_XXX.fastq
# split a string vector by starting with fnFs df and then splitting when you get to "_" 
sample.names <- sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1)
sample.names

# plotting Q of forward and reverse reads
plotQualityProfile(fnFs) #looks good can trim at 15ish and 240ish
plotQualityProfile(fnRs) #looks pretty good up until 200 bp

# Define the paths to filtered files we are going to create
filtFs <- file.path(path.raw, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.raw, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Perform filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=2, 
                     trimLeft=c(20, 20), # REPLACE XXX/YYY with proper parameter choices
                     truncLen=c(245, 200)) # REPLACE XXX/YYY with proper parameter choices
head(out)
tail(out)
#less than 30% of reads were lost which is good!
out <- write.csv(out, file = "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/out.csv")

#################
##Looking at error rates##
################
# Learn the error model from the filtered data.
errF <- learnErrors(filtFs, multi=TRUE) 
# 123149475 total bases in 547331 reads from 6 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multi=TRUE)
# 118629720 total bases in 659054 reads from 7 samples will be used for learning the error rates.

# Visualize the error model. Points are observations, black line is fitted error model.`
plotErrors(errF) #passes
plotErrors(errR) #passes

#################
##Dereplication##
################
## combining all identical sequences into "unique sequences" while keeping Q information and abundance = # of reads with that sequence
# saves computational time by collapsing dataset & average Q scores inform error models for denoising step
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

## Naming the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Run the DADA2 method using the fitted error model.
# pseudo bc there are more samples so saving computation time
ddF <- dada(derepFs, errF, pool= "pseudo" , multi=TRUE) 
ddR <- dada(derepRs, errR, pool="pseudo", multi=TRUE) 

## checking output of dadaFs
ddF[[1]]
  # 2546 sequence variants were inferred from 38651 input unique sequences.
## checking output of dadaRs
ddR[[1]]
  # 2240 sequence variants were inferred from 28606 input unique sequences.

# Merge the denoised forward and reverse reads together
# justConcatenate = TRUE bc less than 25 bp of overlap after filter and trim if region is 460 bp long (based on google)
mm <- mergePairs(ddF, derepFs, ddR, derepRs, verbose=TRUE, justConcatenate=TRUE)

# Construct a sequence table: rows are samples, columns are ASVs, values are abundances.
# ASVs (columns) are ordered by their summed abundance across all samples.
sta <- makeSequenceTable(mm)
dim(sta)
  # 26 samples by 98260 unique ASVs

# Remove chimeric ASVs and construct a new chimera-free sequence table.
st <- removeBimeraDenovo(sta, multi=TRUE, verbose=TRUE)
  # Identified 76464 bimeras out of 98260 input sequences.
sum(st)/sum(sta)
  # 0.9146952 of ASVs are retained
saveRDS(st,"~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/seqtab_nochim.rds" )

####################################################################
######  Inspect the number of reads passing through each step ######
######  THIS IS THE SINGLE MOST IMPORTANT SANITY CHECK!!      ######
####################################################################
# Code derived from the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(ddF, getN), sapply(ddR, getN), sapply(mm, getN), rowSums(st))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- basename(fnFs)
head(track)
write.csv(track, file = "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/track.csv" )

#############
##Assigning Taxonomy##
#############
#load in sequence table
seqtab.nochim <- readRDS("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/seqtab_nochim.rds")
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim,  "~/Library/CloudStorage/Box-Box/ReynoldsMicrobes/reynoldsMicrobes/rawseq/guillaume_rawfastq/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
head(unname(taxa))
saveRDS(taxa, "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/taxa.rds")

#########
##Handoff to Phyloseq##
#############
library("tidyverse")
##loading in metadata
metadata <- read.table("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/data/GSE-meta.txt",sep = "\t", header = TRUE)
str(metadata)
metadata$tank.id <- as.factor(metadata$tank.id)
rownames(metadata)
# convert rownames to column named sample.no
library(tibble)
metadata <- tibble::rownames_to_column(metadata, "sample.no")
# add KZ. to beginning of sample.no
metadata$sample.no <- paste("KZ.", metadata$sample.no, "a", sep = "")
rownames(metadata) <- metadata$sample.no

# load in ASV table
seqtab.nochim <- readRDS("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/seqtab_nochim.rds")
# load in taxa table
taxa <- readRDS("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/taxa.rds")

# check if row names of taxa are the same as column names of seqtab.nochim
idx <- match(rownames(taxa), colnames(seqtab.nochim))
plot(idx)
  # yes they match

# to be safe, let's filter seqtab.nochim to only include the subst of column indices found in the idx vector. 
seqtab.nochim <- seqtab.nochim[,idx]
ncol(seqtab.nochim)

# create df of replacement ASV numbers for sequences and save as a df so you can reference the sequences in the future
ASVseqs <- data.frame("asv" = paste0("ASV", seq(from = 1, to = 21796, by = 1)), "sequence" = rownames(taxa))
saveRDS(ASVseqs, "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/ASVseq.rds")

# rename seqtab.nochim and taxa df so they are easier to interpret
colnames(seqtab.nochim) <- ASVseqs$asv
rownames(taxa) <- ASVseqs$asv

#creating phyloseq object using ASV and taxa table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(taxa))
ps
sample_sums(ps) #a few samples have <10,000 reads
rownames(sample_data(ps)) # nice
rownames(otu_table(ps)) # nice

#####################
##Cleaning the data##
#####################
##removing mitochondia and chloroplast sequences which are not bacteria that we care about for this study
ps_nomitochloro <-subset_taxa(ps, Family !="Mitochondria" & Order !="Chloroplast") 
nomito_filtrationSummary <- as.data.frame(sample_sums(ps_nomitochloro))
saveRDS(nomito_filtrationSummary, file ="~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/sequences/nomito_filtrationSummary.txt" )

#save RDS so you can close out code and load into new file for analysis
saveRDS(ps_nomitochloro, file= "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/data/ps_nomitochloro")

###################
##Decontaminating samples##
###################
# source code from Ben himself: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# setting up
library(phyloseq); packageVersion("phyloseq") # 1.46.0
library(decontam); packageVersion("decontam") # 1.22.0
library(ggplot2); packageVersion("ggplot2") # 3.4.4

ps <- readRDS("~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/data/ps_nomitochloro")

# calculate number of singletons
length(which(taxa_sums(ps) == 1))
  # 3233

head(sample_data(ps))
  # nice

# inspect library sizes
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()
  # looks good where controls have the lowest sample size
colnames(df)
# Identify contaminants- frequency
contamdf.freq <- isContaminant(ps, method="frequency", conc="total.dna.ng.ul.")
head(contamdf.freq)
# This calculation has returned a data.frame with several columns, 
# the most important being $p which containts the probability that was used for 
# classifying contaminants, and $contaminant which contains TRUE/FALSE classification 
# values with TRUE indicating that the statistical evidence that the associated sequence
# feature is a contaminant exceeds the user-settable threshold. 
# As we did not specify the threshold, the default value of threshold = 0.1 was used,
# and $contaminant=TRUE if $p < 0.1.

table(contamdf.freq$contaminant)
# 77 out of 15201 ASVs are classified as contaminants

tail(which(contamdf.freq$contaminant))

# looking at sequence 1 (not contaminant) vs 643 (contaminant)
plot_frequency(ps, taxa_names(ps)[c(1,10102)], conc="total.dna.ng.ul.") + 
  xlab("DNA Concentration (Qubit intensity)")
# sequence 643 fits red line more than 1 which models a contaminant strain

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),8)], conc="total.dna.ng.ul.") +
  xlab("DNA Concentration (Qubit intensity)")

# pruning taxa of contaminants based on frequency
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
ps.noncontam

##### Another way to identify contaminants
# identifying contaminants based on prevalence
sample_data(ps)$is.neg <- sample_data(ps)$sample_or_control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
  # no taxa identified as being a contaminant

otu.table <- data.frame(otu_table(ps))
otu.table.pruned <- data.frame(otu_table(ps.noncontam))

saveRDS(ps.noncontam, file= "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/data/ps_nocontam")

# prune decontamed ps object of NC samples
ps.nomitochloroDecontamed <- subset_samples(ps.noncontam, sample_data(sample_or_control != "control"))
  
  library(phyloseq)

# Assuming 'your_phyloseq_object' is your phyloseq object
# Replace 'sample_to_remove1' and 'sample_to_remove2' with the actual sample names you want to remove

samples_to_remove <- c("KZ.26a", "KZ.13a")

ps.final <- subset_samples(ps.noncontam, !sample_names(ps.noncontam) %in% samples_to_remove)
saveRDS(ps.final, file = "~/Library/CloudStorage/Box-Box/GSE/GrainSizeExperiment/data/GSE-ps-final.rds")
