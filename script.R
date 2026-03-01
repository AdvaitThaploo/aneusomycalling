
# Note: Tested and functional on R 4.4.3 as of 5/21/2025. Choose R version in Tools > Global Options > General 
### Note: You must install RTools as well. https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html

# Installation - run once 

install.packages("BiocManager")

library(BiocManager)

install.packages("devtools")

library(devtools)

install.packages("remotes")

library(remotes)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", ask = FALSE, force = TRUE)
BiocManager::install("Rsamtools", ask = FALSE, force = TRUE)
BiocManager::install("GenomicRanges", ask = FALSE, force = TRUE)
BiocManager::install("bamsignals", ask = FALSE, force = TRUE)
BiocManager::install("DNAcopy", ask = FALSE, force = TRUE)
BiocManager::install("AneuFinderData")
BiocManager::install("AneuFinder")

install.packages(c("foreach", "doParallel", "reshape2", "ggdendro", "mclust"))


install.packages(c(
    "ecp", "tidyverse", "ggsignif", "ggplot2", "dplyr", "ggpubr",
    "scales", "tidyr", "openxlsx", "writexl", "Rtools"))

remotes::install_local("C:/Users/lampsonlab-adm/Downloads/JunMaAneu/ReorderCluster") # Local File

install.packages("C:/Users/lampsonlab-adm/Downloads/aneufinder-developer", 
                 repos = NULL, 
                 type = "source")

# Loading - run every time 

library(BSgenome.Mmusculus.UCSC.mm10)
library(ecp)
library(tidyverse)
library(GenomicRanges)
library(bamsignals)
library(openxlsx)
library(ReorderCluster)
library(writexl)
library(AneuFinderDev)

# Analysis

chromosomes <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
                'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                'chr17', 'chr18', 'chr19', 'chrX')

# Skip if using premade euploid reference Step 1: Run Aneufinder naively in order to find 3 representative euploid cells for each of the 2xB and 5xUb groups to construct euploid references

AneuFinderDev::Aneufinder(inputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xB', outputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/NaiveAneufinderResults2xB', reuse.existing.files = FALSE,
                       pairedEndReads = TRUE, assembly = "mm10", chromosomes = chromosomes,
                       use.bamsignals = TRUE, correction.method = 'GC',
                       GC.BSgenome = "BSgenome.Mmusculus.UCSC.mm10",
                       method = "edivisive", min.ground.ploidy = 1.5, max.ground.ploidy = 2.5)

AneuFinderDev::Aneufinder('C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUb', 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/NaiveAneufinderResults5xUb', reuse.existing.files = FALSE,
                       pairedEndReads = TRUE, assembly = "mm10", chromosomes = chromosomes,
                       use.bamsignals = TRUE, correction.method = 'GC',
                       GC.BSgenome = "BSgenome.Mmusculus.UCSC.mm10",
                       method = "edivisive")

# This will produce two output folders. The next step is to cluster the cells by quality, and 
# then manually inspect the profiles plot to select for 3 euploid references from the first cluster 
# with high reads and good appearance to merge into the euploid reference.

files <- list.files('C:/Users/lampsonlab-adm/Downloads/JunMaAneu/NaiveAneufinderResults2xB/MODELS/method-edivisive/', full.names = TRUE)
Clusters <- clusterByQuality(files)
Cluster1 <- Clusters$classification[[1]]
files <- sapply(strsplit(basename(Cluster1), '_binsize'), '[[', 1)

files  

# Set files to merge (these are the ones I used to produce my 2xB reference)
filesmerge <- c("B95.bam", "B96.bam", "B99.bam")

# Merge BAMs
mergeBam(filesmerge, destination = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xBeuploidreference.bam')


files <- list.files('C:/Users/lampsonlab-adm/Downloads/JunMaAneu/NaiveAneufinderResults5xUb/MODELS/method-edivisive/', full.names = TRUE)
Clusters <- clusterByQuality(files)
Cluster1 <- Clusters$classification[[1]]
files <- sapply(strsplit(basename(Cluster1), '_binsize'), '[[', 1)

files  

# Set files to merge (these are the ones I used to produce my 2xB reference)

filesmerge <- c("B10.bam", "B12.bam", "B13.bam")

# Merge BAMs
mergeBam(filesmerge, destination = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUbeuploidreference.bam')


# Bin reads
bins <- binReads("C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xBeuploidreference.bam", assembly = 'mm10', binsizes = 500e3, 
                 chromosomes = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9', 'chr10', 
                                 'chr11','chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'))[[1]]

# Blacklist generation
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.999)
blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist <- GenomicRanges::reduce(blacklist)

blacklist.file <- "C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xBblacklist.bed.gz"
exportGRanges(blacklist, filename = blacklist.file, header = FALSE, chromosome.format = 'UCSC')



# Bin reads
bins <- binReads("C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUbeuploidreference.bam", assembly = 'mm10', binsizes = 500e3, 
                 chromosomes = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9', 'chr10', 
                                 'chr11','chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'))[[1]]

# Blacklist generation
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.999)
blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist <- GenomicRanges::reduce(blacklist)

blacklist.file <- "C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUbblacklist.bed.gz"
exportGRanges(blacklist, filename = blacklist.file, header = FALSE, chromosome.format = 'UCSC')


# Run Aneufinder - skip to here if using premae euploid reference and blacklist
chromosomes <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9', 'chr10', 
                 'chr11','chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX')


AneuFinderDev::Aneufinder(
  inputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xB',
  outputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/FinalAneufinderResults2xB',
  reuse.existing.files = FALSE,
  binsizes = 2e+06,
  variable.width.reference = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xBeuploidreference.bam',
  pairedEndReads = TRUE,
  assembly = "mm10",
  chromosomes = chromosomes,
  blacklist = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/2xBblacklist.bed.gz',
  use.bamsignals = TRUE,
  correction.method = 'GC',
  GC.BSgenome = "BSgenome.Mmusculus.UCSC.mm10",
  method = "edivisive",
  min.ground.ploidy = 1.5,
  max.ground.ploidy = 2.5
)

AneuFinderDev::Aneufinder(
  inputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUb',
  outputfolder = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/FinalAneufinderResults5xUb',
  reuse.existing.files = FALSE,
  binsizes = 2e+06,
  variable.width.reference = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUbeuploidreference.bam',
  pairedEndReads = TRUE,
  assembly = "mm10",
  chromosomes = chromosomes,
  blacklist = 'C:/Users/lampsonlab-adm/Downloads/JunMaAneu/5xUbblacklist.bed.gz',
  use.bamsignals = TRUE,
  correction.method = 'GC',
  GC.BSgenome = "BSgenome.Mmusculus.UCSC.mm10",
  method = "edivisive",
  min.ground.ploidy = 1.5,
  max.ground.ploidy = 2.5
)
