library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(raster)

rm(list = ls())

#Read in seurat object with updated clustering
allSeurat <- readRDS('/Users/Brian/Downloads/seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds')

############## Make txt files for COMET analysis ##############################
#Specify the genes we need to use (flow receptors we have profiled)
panelA <- c("IL3RA", "KIT", "MPL", "HAVCR2", "CSF2RA", 'CCR2', 'IL2RG', 'IL5RA',
            'IL15RA', 'TLR2', 'FLT3', 'CSF1R', 'CXCR4', 'CSF3R', 'CXCR2',
            'CXCR1', 'IL18R1', 'IFNGR1', 'TLR4', 'TNFRSF1A', 'TNFRSF1B', 'IL6R',
            'CD34', 'CD38', 'PTPRC', 'IL3RA', 'THY1', 'ITGA6')

dna_seq_genes <- read.csv("~/scCode/gene_names_dna_seq.csv")

for (i in 1:length(dna_seq_genes[,1])) {
  if (! dna_seq_genes$Gene.Name[i] %in% allSeurat@assays[["RNA"]]@data@Dimnames[[1]]) {
    print(dna_seq_genes$Gene.Name[i])
  }
}

all_genes <- unique(c(panelA, dna_seq_genes$Gene.Name))


#Get cluster 2 data and visualization (umap)
clus2.data <- allSeurat@assays[['RNA']]@data[,(allSeurat@meta.data[["RNA_snn_res.0.05"]] == "2")]
clus2.umap <- allSeurat@reductions[["umap"]]@cell.embeddings[allSeurat@meta.data[["RNA_snn_res.0.05"]] == "2",]

#Get non cluster 2 data and visualization
nonClus2.data <- allSeurat@assays[['RNA']]@data[,(allSeurat@meta.data[["RNA_snn_res.0.05"]] != "2")]
nonClus2.umap <- allSeurat@reductions[["umap"]]@cell.embeddings[allSeurat@meta.data[["RNA_snn_res.0.05"]] != "2",]

#Max number of cells for COMET tool is 65k, will have to downsample nonclus2 data to be under that
nSamples <- 64999-dim(clus2.data)[2]

#Assign cluster arbitrarily, clus 2 becomes 1. All others become 0.
nonClus2.cluster <- data.frame("X" = rep(0, nSamples))
clus2.cluster <- data.frame("X" = rep(1, dim(clus2.data)[2]))

#Convert to df's and downsample nonClus2 data
clus2.data.df <- data.frame(clus2.data[all_genes,])
nonClus2.data.df <- data.frame(nonClus2.data[all_genes,])
cellsKeep <- sample(length(nonClus2.data.df), nSamples)
nonClus2.data.DS <- nonClus2.data.df[,cellsKeep]
nonClus2.umap.DS <- data.frame(nonClus2.umap[cellsKeep,])

#Check the names of cells to make sure they match, Note: they don't but we know why ("X" added in front)
#colnames(nonClus3.data.DS)[!(colnames(nonClus3.data.DS) %in% row.names(nonClus3.umap.DS))]
#row.names(nonClus3.umap.DS)[!row.names(nonClus3.umap.DS) %in% colnames(nonClus3.data.DS)]

markers.df <- cbind(clus2.data.df, nonClus2.data.DS)
umap.df <- rbind(clus2.umap, nonClus2.umap.DS)
clusters.df <- rbind(clus2.cluster, nonClus2.cluster)
row.names(clusters.df) <- row.names(umap.df)
#Get rid of the weird "X" in front of the names introduced by R for ones that start with number
colnames(markers.df) <- row.names(umap.df)

#Write to outfiles, specifying the format needed for input to COMET
outdir <- "/Users/Brian/scCode/Comet_04252022/"
write.table(markers.df,paste0(outdir, "res=0.05_markers.txt"),sep="\t",row.names=T, quote = F)
write.table(umap.df,paste0(outdir, "res=0.05_umap.txt"),sep="\t",row.names=T, col.names =F, quote = F)
write.table(clusters.df, paste0(outdir, "res=0.05_clusters.txt"),sep="\t",row.names=T, col.names =F, quote = F)
write.table(all_genes, paste0(outdir, "res=0.05_genes.txt"),sep="\t",row.names=F, col.names =F, quote = F)
