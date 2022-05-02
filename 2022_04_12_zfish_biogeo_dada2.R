# Title and author information --------------------------------------------
#!/usr/bin/R

#################################################
#                                               #
#          2022_04_12_zfish_biogeo_dada2.R     #
#                                               #
#################################################

#Title: Zebrafish microbiome biogeography
#
#Copyright (C) 2022-2023  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

#Purpose: The purpose of this script is to generate ASV profiles for
# zebrafish biogeography analysis

# PACKAGE INSTALLATION ----------------------------------------------------

# If you haven't already, install the required software. The commands in this
# section are run only once and should be commented after you have finished running
# them

# Installing Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# Installing dada2
BiocManager::install("dada2")

# Installing other packages needed
install.packages("ggplot2")
install.packages("vegan")

# Unpacking and sorting files ---------------------------------------------

# The following commands can be run in a terminal window (Mac) or in Rstudio console
# window. The tar archive should be unpacked in a suitable location
# (i.e., not the desktop). Preferably you have a robust file naming system for
# your research. For example, my research analysis are all kept in
# ~/Documents/research/<name_of_analysis_project>. Each project has the following
# subdirectories: data/, analysis, scripts/. Scripts contains the R project. Other
# directories can be added as needed (e.g., manuscript, reports, ppt, etc.). This
# compartmentalizes the projects and helps keep things organized


#decompress the and untar the archive
#tar -xvf zebrafish_biogeography_filtered.tar.gz

#NOT RUN

# Used this to move all empty files, and those with vanishingly small number of reads
# to a dump dir
#find . -maxdepth 1 -type f -size -100 -exec mv {} dumps/ \;


# SET ENVIRONMENT ---------------------------------------------------------

#load required libraries

library(dada2)
library(ggplot2)
library(vegan)

options(stringsAsFactors = F)


# Functions ---------------------------------------------------------------

###
#              Function filter_df                #
###

filter <- function(df) {
  df <- df[which(rowSums(df) > 0), which(colSums(df) > 0)]
  return(df)
}

###
#             Function normalize                 #
###

#normalize counts either by relative abundance of rarefying
#depends on vegan

normalize <- function(df, method="rel", depth=depth){
  #default method = relative abundance
  if(method == "rare"){
    if( is.null(depth)){
      depth <- min(rowSums(df))
    }
    ndf <- rrarefy(df, depth)
    ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }else{
    ndf <- sweep(df,
                 1,
                 rowSums(df),
                 `/`)
    #  ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }
  return(ndf)
}


# IMPORT DATA -------------------------------------------------------------

path <- "/Users/cgaulke/unsynced_projects/raw_data/2022_04_12_zebrafish_biogeo/"
filt.path <- "/Users/cgaulke/Documents/research/zfish_biogeography/data/filtered_data/" #filtered file directory make sure to update

#git ignore to ignore this
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# # ANALYSIS: QC ------------------------------------------------------------

#make sure the lengths are the same
length(fnFs) == length(fnRs)

#get sample names, in our case we have to leave some extra on the end
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

#preview
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

#aggregate all data together now
fnFs_qual.plot <- plotQualityProfile(fnFs,aggregate = T)
fnRs_qual.plot <- plotQualityProfile(fnRs,aggregate = T)

#set up for filtering
filtFs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_R_filt.fastq.gz"))

#make these named
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim
#by looking at the data we see a rapid drop in quality around 250bp R1 (200 R2). Since the
#average drops below ~30 around 250 we will truncate at 200 for the reverse. The
#forward looks better (this is usual) so we will truncate around 250. We will
#also take off about 10 bases on the left as these bases are highly skewed.

#note the original tutorial uses generic variable names

#stopped here

filter.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = 10,
                     compress=TRUE, multithread=TRUE)
#take a look
View(filter.out)

colMeans(filter.out) #mean number of reads in and out of filtering
mean(1-(filter.out[,2]/filter.out[,1])) #mean % filtered reads

#lets look at these numbers in a little more detail
fivenum(1-(filter.out[,2]/filter.out[,1])) #five number report for % filtered
hist(1-(filter.out[,2]/filter.out[,1]))  # hist #mean % filtered reads

# since some libraries look like there is a higher level of filtration
# than others lets take a closer look at this

sort(1-(filter.out[,2]/filter.out[,1]))

#let's keep this in mind moving forward

# # ANALYSIS: ERROR ---------------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# # ANALYSIS: MERGE AND FILTER -----------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#how much data was wrapped up in chimeras
sum(seqtab.nochim)/sum(seqtab)

#good to have full path
dir.create("/Users/cgaulke/Documents/research/zfish_biogeography/data/dada2/")

write.table(seqtab.nochim,
            file = "../../data/dada2/seqtab_nochim.txt",
            quote = FALSE,
            sep = "\t"
            )

# # ANALYSIS: TRACK READS ---------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(filter.out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# # ANALYSIS: ADD TAX -------------------------------------------------------

#***Here is where you will need to go and download the silva databases.
#***Be sure to get the right ones (the names are the same as the ones below)
#***These files can be downloaded here:https://zenodo.org/record/4587955#.YSlzKC1h1hA


taxa <- assignTaxonomy(seqtab.nochim,
          "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_nr99_v138.1_train_set.fa",
          multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_species_assignment_v138.1.fa")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
