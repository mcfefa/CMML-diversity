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

#Load the seurat object and split by sample
allSeurat <- readRDS('/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringOldUMAP_04-30-2021_all_46_samples_umap_layered_noHarmony.rds')
split.list <- SplitObject(allSeurat, split.by = "orig.ident")

#Ugly but safest way to make sure the patients are aligned in the correct order
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

#Include the genes that we want to visualize
gatingGenes <- c("PTPRC", "CD34", "CD38", "IL3RA", "THY1", "ITGA6")
miscGenes <- c("MPO", "JUN", "CD79A", "PLAUR", "GATA1", "IRF8", "EBF1", "CLEC12A", "IL10",
               "IL1A", "IL1B", "IL1RN", "FOS", "ATF1", "JDP2", "SEMA4D")
panelA <- c('CCR2', 'CSF2RA', 'HAVCR2', 'TLR2', "TNFRSF1B", 'CSF1R', 'IFNGR1', 'TLR4', 'CSF3R',
            "IL3RA", "KIT", "MPL", 'IL2RG', 'IL5RA',
            'IL15RA', 'FLT3', 'CXCR4', 'CXCR2',
            'CXCR1', 'IL18R1', 'TNFRSF1A', 'IL6R')

for (i in 1:39){
  #Get cmml seurat object for individual sample
  cmml <- split.list[[7+i]]
  
  #Find variable features
  cmml <- FindVariableFeatures(cmml)
  #Include variable genes as well as receptor genes (panel A) and genes used to gate receptors
  #Gating Genes: PTPRC, CD34, CD38, IL3RA, THY1, ITGA6
  cmml.counts <- t(cmml@assays[['RNA']]@counts)[,unique(c(VariableFeatures(cmml),gatingGenes, 
                                                            miscGenes, panelA))]
  cmml.matrix <- as.matrix(cmml.counts)
  #Write count matrix to csv
  write.csv(cmml.matrix, paste0('/Users/4472241/scCode/runPalantir/countData/cmml',i,'Counts.csv'))
  
  #Write cell names of cells in cluster 2 to csv
  namesInClus2 <- row.names(cmml.matrix)[as.numeric(as.character(cmml@meta.data[["RNA_snn_res.0.05"]]))==2]
  write.csv(namesInClus2, paste0('/Users/4472241/scCode/runPalantir/countData/cmml',i,'_namesInClus2.csv'))
}

for (i in 1:7){
  #Get seurat object of individual sample, find var. features and add genes we want to visualize
  normal <- split.list[[i]]
  normal <- FindVariableFeatures(normal)
  normal.counts <- t(normal@assays[['RNA']]@counts)[,unique(c(VariableFeatures(normal),gatingGenes, 
                                                              miscGenes, panelA))]
  normal.matrix <- as.matrix(normal.counts)
  #Write count matrix to csv
  write.csv(normal.matrix, paste0('/Users/4472241/scCode/runPalantir/countData/normal',i,'Counts.csv'))
  
  #Get cluster 2 cells and write their cell names to csv
  namesInClus2 <- row.names(cmml.matrix)[as.numeric(as.character(cmml@meta.data[["RNA_snn_res.0.05"]]))==2]
  write.csv(namesInClus2, paste0('/Users/4472241/scCode/runPalantir/countData/normal',i,'_namesInClus2.csv'))
}
