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

allSeurat <- readRDS('/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringOldUMAP_04-30-2021_all_46_samples_umap_layered_noHarmony.rds')
split.list <- SplitObject(allSeurat, split.by = "orig.ident")
split.list_ordered <- list()
split.list_ordered[["Normal_1"]] <- split.list[["SettyPt1"]]
split.list_ordered[["Normal_2"]] <- split.list[["SettyPt2"]]
split.list_ordered[["Normal_3"]] <- split.list[["SettyPt3"]]
split.list_ordered[["Normal_4"]] <- split.list[["CD34"]]
split.list_ordered[["Normal_5"]] <- split.list[["HuaPt1"]]
split.list_ordered[["Normal_6"]] <- split.list[["HuaPt2"]]
split.list_ordered[["Normal_7"]] <- split.list[["HuaPt4"]]
split.list_ordered[["CMML_1"]] <- split.list[["LTB3966"]]
split.list_ordered[["CMML_2"]] <- split.list[["LTB4121"]]
split.list_ordered[["CMML_3"]] <- split.list[["LTB6169"]]
split.list_ordered[["CMML_4"]] <- split.list[["4J003"]]
split.list_ordered[["CMML_5"]] <- split.list[["4K001"]]
split.list_ordered[["CMML_6"]] <- split.list[["4Q001"]]
split.list_ordered[["CMML_7"]] <- split.list[["5E001"]]
split.list_ordered[["CMML_8"]] <- split.list[["5H001"]]
split.list_ordered[["CMML_9"]] <- split.list[["SF100109106293"]]
split.list_ordered[["CMML_10"]] <- split.list[["SF100109111451"]]
split.list_ordered[["CMML_11"]] <- split.list[["SF100109110236"]]
split.list_ordered[["CMML_12"]] <- split.list[["SF14040100158"]]
split.list_ordered[["CMML_13"]] <- split.list[["SF14060200025"]]
split.list_ordered[["CMML_14"]] <- split.list[["SF12062800475"]]
split.list_ordered[["CMML_15"]] <- split.list[["SF14072200012"]]
split.list_ordered[["CMML_16"]] <- split.list[["SF13061200056"]]
split.list_ordered[["CMML_17"]] <- split.list[["4S001"]]
split.list_ordered[["CMML_18"]] <- split.list[["2V001"]]
split.list_ordered[["CMML_19"]] <- split.list[["SF14101000049"]]
split.list_ordered[["CMML_20"]] <- split.list[["SF16112900158"]]
split.list_ordered[["CMML_21"]] <- split.list[["6AE001"]]
split.list_ordered[["CMML_22"]] <- split.list[["6AC001"]]
split.list_ordered[["CMML_23"]] <- split.list[["SF100109101914"]]
split.list_ordered[["CMML_24"]] <- split.list[["SF12042500035"]]
split.list_ordered[["CMML_25"]] <- split.list[["6AD001"]]
split.list_ordered[["CMML_26"]] <- split.list[["SF12092600014"]]
split.list_ordered[["CMML_27"]] <- split.list[["SF14031800065"]]
split.list_ordered[["CMML_28"]] <- split.list[["SF14050700419"]]
split.list_ordered[["CMML_29"]] <- split.list[["SF16026800045"]]
split.list_ordered[["CMML_30"]] <- split.list[["SF16072200003"]]
split.list_ordered[["CMML_31"]] <- split.list[["SF16112300029"]]
split.list_ordered[["CMML_32"]] <- split.list[["SF13032800016"]]
split.list_ordered[["CMML_33"]] <- split.list[["SF14110400108"]]
split.list_ordered[["CMML_34"]] <- split.list[["SF14111400033"]]
split.list_ordered[["CMML_35"]] <- split.list[["SF14092500135"]]
split.list_ordered[["CMML_36"]] <- split.list[["SF14061300036"]]
split.list_ordered[["CMML_37"]] <- split.list[["SF14080400065"]]
split.list_ordered[["CMML_38"]] <- split.list[["SF15010200008"]]
split.list_ordered[["CMML_39"]] <- split.list[["SF13070900171"]]

names(split.list_ordered)

gatingGenes <- c("PTPRC", "CD34", "CD38", "IL3RA", "THY1", "ITGA6")
miscGenes <- c("MPO", "JUN", "CD79A", "PLAUR", "GATA1", "IRF8", "EBF1", "CLEC12A", "IL10",
               "IL1A", "IL1B", "IL1RN", "FOS", "ATF1", "JDP2", "SEMA4D")
panelA <- c('CCR2', 'CSF2RA', 'HAVCR2', 'TLR2', "TNFRSF1B", 'CSF1R', 'IFNGR1', 'TLR4', 'CSF3R',
            "IL3RA", "KIT", "MPL", 'IL2RG', 'IL5RA',
            'IL15RA', 'FLT3', 'CXCR4', 'CXCR2',
            'CXCR1', 'IL18R1', 'TNFRSF1A', 'IL6R')

#Make csv of the whole thing
#allSeurat <- FindVariableFeatures(allSeurat)
#Include variable genes as well as receptor genes (panel A) and genes used to gate receptors
#Gating Genes: PTPRC, CD34, CD38, IL3RA, THY1, ITGA6
#allSeurat.counts <- t(allSeurat@assays[['RNA']]@counts)[,unique(c(VariableFeatures(allSeurat),gatingGenes, 
#                                                        miscGenes, panelA))]
#allSeurat.matrix <- as.matrix(allSeurat.counts)
#write.csv(allSeurat.matrix, paste0('/Users/4472241/scCode/runPalantir/countData/cmml_All_Counts.csv'))

#namesInClus3_9 <- row.names(allSeurat.matrix)[c(as.numeric(as.character(allSeurat@meta.data[["clusters_noHarmony_res.0.05_30Neighbors"]]))==3,
#                                           as.numeric(as.character(allSeurat@meta.data[["clusters_noHarmony_res.0.05_30Neighbors"]]))==9)]
#write.csv(namesInClus3_9, paste0('/Users/4472241/scCode/runPalantir/countData/cmml_All_namesInClus3_9.csv'))


for (i in 1:39){
  cmml <- split.list[[7+i]]
  
  cmml <- FindVariableFeatures(cmml)
  #Include variable genes as well as receptor genes (panel A) and genes used to gate receptors
  #Gating Genes: PTPRC, CD34, CD38, IL3RA, THY1, ITGA6
  cmml.counts <- t(cmml@assays[['RNA']]@counts)[,unique(c(VariableFeatures(cmml),gatingGenes, 
                                                            miscGenes, panelA))]
  cmml.matrix <- as.matrix(cmml.counts)
  write.csv(cmml.matrix, paste0('/Users/4472241/scCode/runPalantir/countData/cmml',i,'Counts.csv'))
  
  namesInClus2 <- row.names(cmml.matrix)[as.numeric(as.character(cmml@meta.data[["RNA_snn_res.0.05"]]))==2]
  write.csv(namesInClus2, paste0('/Users/4472241/scCode/runPalantir/countData/cmml',i,'_namesInClus2.csv'))
}

for (i in 1:7){
  normal <- split.list[[i]]
  normal <- FindVariableFeatures(normal)
  normal.counts <- t(normal@assays[['RNA']]@counts)[,unique(c(VariableFeatures(normal),gatingGenes, 
                                                              miscGenes, panelA))]
  normal.matrix <- as.matrix(normal.counts)
  write.csv(normal.matrix, paste0('/Users/4472241/scCode/runPalantir/countData/normal',i,'Counts.csv'))
  
  namesInClus2 <- row.names(cmml.matrix)[as.numeric(as.character(cmml@meta.data[["RNA_snn_res.0.05"]]))==2]
  write.csv(namesInClus2, paste0('/Users/4472241/scCode/runPalantir/countData/normal',i,'_namesInClus2.csv'))
}

check <- allSeurat@assays[['RNA']]@counts[c(gatingGenes, miscGenes, panelA),]





############## Make txt files for COMET analysis ##############################
#split_by_cluster <- SplitObject(allSeurat, split.by = "RNA_snn_res.0.05")

panelA <- c("IL3RA", "KIT", "MPL", "HAVCR2", "CSF2RA", 'CCR2', 'IL2RG', 'IL5RA',
            'IL15RA', 'TLR2', 'FLT3', 'CSF1R', 'CXCR4', 'CSF3R', 'CXCR2',
            'CXCR1', 'IL18R1', 'IFNGR1', 'TLR4', 'TNFRSF1A', 'TNFRSF1B', 'IL6R',
            'CD34', 'CD38', 'PTPRC', 'IL3RA', 'THY1', 'ITGA6')

clus2.data <- allSeurat@assays[['RNA']]@data[,(allSeurat@meta.data[["RNA_snn_res.0.05"]] == "2")]
clus2.umap <- allSeurat@reductions[["umap"]]@cell.embeddings[allSeurat@meta.data[["RNA_snn_res.0.05"]] == "2",]

nonClus2.data <- allSeurat@assays[['RNA']]@data[,(allSeurat@meta.data[["RNA_snn_res.0.05"]] != "2")]

nonClus2.umap <- allSeurat@reductions[["umap"]]@cell.embeddings[allSeurat@meta.data[["RNA_snn_res.0.05"]] != "2",]

nSamples <- 64999-dim(clus2.data)[2]

nonClus2.cluster <- data.frame("X" = rep(0, nSamples))
clus2.cluster <- data.frame("X" = rep(1, dim(clus2.data)[2]))

clus2.data.df <- data.frame(clus2.data[panelA,])
nonClus2.data.df <- data.frame(nonClus2.data[panelA,])
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


write.table(markers.df,"/Users/4472241/scCode/CometNewClustering/isolateClus2/res=0.05_markers.txt",sep="\t",row.names=T, quote = F)
write.table(umap.df,"/Users/4472241/scCode/CometNewClustering/isolateClus2/res=0.05_umap.txt",sep="\t",row.names=T, col.names =F, quote = F)
write.table(clusters.df, "/Users/4472241/scCode/CometNewClustering/isolateClus2/res=0.05_clusters.txt",sep="\t",row.names=T, col.names =F, quote = F)
write.table(panelA, "/Users/4472241/scCode/CometNewClustering/isolateClus2/res=0.05_genes.txt",sep="\t",row.names=F, col.names =F, quote = F)


