library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(pheatmap)
library(ggplot2)
library("reticulate")
Sys.setenv(RETICULATE_PYTHON = "./anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

# INPUT: A bam file and a bai file for each single cell 
# GOAL: Generate count tables for mt sites, which are used for mitoclone input

# Directions to run on cluster:
# screen
# request job
# conda activate mitoclone
# open R
# Run this file

rm(list=ls())

# This is how we would generate count tables from individual bam files
# Note: you can skip this step because we provide its output: "countTablesAllForMitoclone.rds"
# We had to run this in four sets of ten due to memory constraints
# First set of ten:
bamDirFirstTen <- './firstTenSamples_singleCell_BAM/'

# Get bam files, exclude bai files (which are in same dir), but make sure each bam has a bai
bamfiles <- list.files(bamDirFirstTen, pattern = "bam", full.names = T)
baifiles <- list.files(bamDirFirstTen, pattern = 'bai', full.names = T)
bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
baifilesCut <- gsub(".bai\\.*","",baifiles)
bamfiles <- bamfiles[bamfiles %in% baifilesCut]

# Generate count table (list of matrices with one matrix for each bam file, so one for each cell)
countTables_firstTen <- baseCountsFromBamList(bamfiles,  sites = "MT:1-16569", 
                                              ncores = 12) #This will take a while
# Get cell names for this set of ten (from single cell bam file names)
namesFrom_bamFiles_firstTen <- sub(".*singleCell_BAM/", "", bamfiles)

# Repeat previous steps for second, third, and fourth set of ten samples, then:
# Merge the countTables and save the larger RDS file (to be used as mitoclone input)
countTables_allSamples <- c(countTables_firstTen, countTables_secondTen, countTables_thirdTen, countTables_fourthTen)
# Add cell names for each cell to the countTables
namesFrom_bamFiles <- c(namesFrom_bamFiles_firstTen, namesFrom_bamFiles_secondTen, namesFrom_bamFiles_thirdTen, namesFrom_bamFiles_fourthTen)
names(countTables_allSamples) <- namesFrom_bamFiles

# Save output (used as input to mitoclone)
saveRDS(countTables_allSamples, './countTablesAllForMitoclone.rds')
