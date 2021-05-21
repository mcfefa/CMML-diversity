library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(Rmagic)
library(gdata)

rm(list = ls())

#Import the seurat object with all of our samples assigned to clusters
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Try the nature paper hsc list 
EppertGeneList <- read.xls('/Users/4472241/scCode/eppert/natureEppert_HSC+LSC_signatures.xls', header = T, skip = 1, sheet = 8)
genes_hscEp <- EppertGeneList$Gene.Symbol[c(1:53, 57, 59:105)]
genes_hscEp[2] <- 'ABCB4'
genes_hscEp[19] <- 'CFH'
genes_hscEp[48] <- 'GNL1'
genes_hscEp[55] <- 'HIST2H2AA3'
genes_hscEp[56] <- "HLA.DDRB4"
genes_hscEp[71] <- "SERINC2"

#Check which eppert genes are in the seurat object
genesKeep <- toupper(genes_hscEp)[toupper(genes_hscEp) %in% allSeurat@assays[['RNA']]@data@Dimnames[[1]]]

#write to csv
write.csv(genesKeep, '/Users/4472241/scCode/eppert/eppertGeneList.csv')

#Add module score based on eppert gene signature, then rename it
allSeurat <- AddModuleScore(allSeurat, list(genesKeep))
allSeurat$EppertHSC <- allSeurat$Cluster1

#Feature plot of HSC score
pdf('/Users/4472241/scCode/eppert/eppertHSC_featurePlot.pdf')
FeaturePlot(allSeurat, features = "EppertHSC")
dev.off()

df <- data.frame("Sample" = allSeurat$orig.ident, "EppertScore" = allSeurat$EppertHSC)

#Split by sample
eppert.split <- split(df, df$Sample)

out.df <- data.frame("Sample" = NA, "HSC_Score" = NA)
for (i in 1:length(eppert.split)){
  sample <- eppert.split[[i]]$Sample[1]
  score <- mean(as.numeric(eppert.split[[i]]$EppertScore))
  out.df[i,] <- c(sample, score)
}

#Make into dataframe
out.df <- data.frame("Sample" = NA, "HSC_Score" = NA)
for (i in 1:length(eppert.split)){
  sample <- eppert.split[[i]]$Sample[1]
  score <- mean(eppert.split[[i]]$eppertHSCScore)
  out.df[i,] <- c(sample, score)
} 

normalsOnly <- out.df$HSC_Score[out.df$Sample %in% c('HuaPt2', "HuaPt4", "CD34", "HuaPt1", "HuaPt3",
                                                     "SettyPt1", "SettyPt2", "SettyPt3")]

out.df$Z_HSC_Score <- (as.numeric(out.df$HSC_Score)-mean(as.numeric(normalsOnly)))/sd(as.numeric(normalsOnly))

write.csv(out.df, '/Users/4472241/scCode/eppert/eppertForHeatmap_normalizedByNormals_withHuaPt3.csv')


#Read in csv with score
eppertScore <- read.csv('/Users/4472241/scCode/eppert/eppertHSCScore.csv')

eppert.split <- split(eppertScore, eppertScore$Sample)

out.df <- data.frame("Sample" = NA, "HSC_Score" = NA)
for (i in 1:length(eppert.split)){
  sample <- eppert.split[[i]]$Sample[1]
  score <- mean(eppert.split[[i]]$eppertHSCScore)
  out.df[i,] <- c(sample, score)
} 

normalsOnly <- out.df$HSC_Score[out.df$Sample %in% c('HuaPt2', "HuaPt4", "CD34", "HuaPt1", "HuaPt3",
                                                          "SettyPt1", "SettyPt2", "SettyPt3")]

out.df$Z_HSC_Score <- (as.numeric(out.df$HSC_Score)-mean(as.numeric(normalsOnly)))/sd(as.numeric(normalsOnly))

write.csv(out.df, '/Users/4472241/scCode/eppert/eppertForHeatmap_normalizedByNormals.csv')

hist(eppert.split[["4K001"]]$eppertHSCScore)
