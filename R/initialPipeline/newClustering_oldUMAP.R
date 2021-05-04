library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(uwot)
library(data.table)
library(matrixStats)
library(clustree)

#Clear workspace
rm(list=ls())

#Set directory to save to
filedir <- "/Users/4472241/scCode/finalLayerUMAP_Code+Output/allSamples/"
tag <- '_noHarmony_30Neighbors_AllSamples'

#Open the fully umapped file
allSeurat <- readRDS(paste0(filedir, 'all_46_samples_umap_layered', tag, '.rds'))

#Remove old clustering
allSeurat@meta.data$RNA_snn_res.0 <- NULL
allSeurat@meta.data$RNA_snn_res.0.01 <- NULL
allSeurat@meta.data$RNA_snn_res.0.025 <- NULL
allSeurat@meta.data$RNA_snn_res.0.05 <- NULL
allSeurat@meta.data$RNA_snn_res.0.1 <- NULL
allSeurat@meta.data$RNA_snn_res.0.15 <- NULL
allSeurat@meta.data$RNA_snn_res.0.2 <- NULL
allSeurat@meta.data$RNA_snn_res.0.25 <- NULL
allSeurat@meta.data$RNA_snn_res.0.3 <- NULL
allSeurat@meta.data$RNA_snn_res.0.35 <- NULL
allSeurat@meta.data$RNA_snn_res.0.4 <- NULL
allSeurat@meta.data$RNA_snn_res.0.45 <- NULL
allSeurat@meta.data$RNA_snn_res.0.5 <- NULL
allSeurat@meta.data$RNA_snn_res.0.55 <- NULL
allSeurat@meta.data$RNA_snn_res.0.6 <- NULL

#Set new clustering at various resolutions
allSeurat <- FindNeighbors(allSeurat)
res_vec <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15)
for (res in res_vec){
  allSeurat <- FindClusters(allSeurat, resolution = res)
  pdf(paste0('/Users/4472241/scCode/setNewClustering+keepOldUMAP/res=', res, '_Seuratv4_defaults.pdf'), width=10, height=6)
  print(DimPlot(allSeurat, reduction = 'umap', group.by = 'seurat_clusters'))
  dev.off()
}

## Grouping Cells as CMML, HMA, or Normal
allSeurat$sampleType <- plyr::mapvalues(
  x = allSeurat$orig.ident, 
  from = c("SettyPt1", "SettyPt2", "SettyPt3", "CD34","HuaPt1","HuaPt2","HuaPt4","LTB3966","LTB4121","LTB6169","4J003","4K001","4Q001","5E001","5H001","SF100109106293","SF100109111451","SF100109110236","SF14040100158","SF14060200025","SF12062800475","SF14072200012","SF13061200056","4S001","2V001","SF14101000049","SF16112900158","6AE001","6AC001","6AD001","SF100109101914","SF12042500035","SF12092600014","SF14031800065","SF14050700419","SF16026800045","SF16072200003","SF16112300029",'SF15010200008', 'SF14111400033','SF14061300036','SF13032800016', 'SF13070900171', 'SF14092500135', 'SF14110400108', 'SF14080400065'), 
  to = c("Normal", "Normal", "Normal", "Normal","Normal","Normal","Normal","CMML","CMML","CMML","Normal","HMA","HMA","HMA","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","HMA","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML", 'RUX', 'RUX', 'Chemo', 'HMA', 'HMA', 'HMA', 'HMA', 'HMA')
)

## Grouping Cells as Normal or Malignant
allSeurat$HMA_Normal_Other <- plyr::mapvalues(
  x = allSeurat$sampleType, 
  from = c("Normal", 'HMA', 'CMML', 'Chemo', 'RUX'), 
  to = c("Normal", 'HMA', 'CMML', 'CMML', 'CMML')
)

#Plot UMAP showing HMA, Normal, CMML/Other(Rux, Chemo)
pdf(paste0('/Users/4472241/scCode/setNewClustering+keepOldUMAP/HMA_or_Norm_or_Other_coloredUMAP.pdf'), width=9, height=7)
print(DimPlot(allSeurat, reduction = 'umap', group.by = 'HMA_Normal_Other', cols = c('grey', 'red', 'black')))
dev.off()

#Plot UMAP of resolution = 0.05 clustering
pdf(paste0('/Users/4472241/scCode/setNewClustering+keepOldUMAP/oldUmap_newClustering_res=0.05.pdf'), width=9, height=7)
print(DimPlot(allSeurat, reduction = 'umap', group.by = 'RNA_snn_res.0.05'))
dev.off()

#Plot clustree
pdf('/Users/4472241/scCode/setNewClustering+keepOldUMAP/newClustering_all46_ClustreePlot.pdf', width = 9, height = 7)
clustree(allSeurat, prefix = "RNA_snn_res.")
dev.off()

#Save the seurat file with the new clustering
saveRDS(allSeurat, '/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringOldUMAP_04-30-2021_all_46_samples_umap_layered_noHarmony.rds')
#allSeurat <- readRDS('/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringOldUMAP_04-30-2021_all_46_samples_umap_layered_noHarmony.rds')

#Extract data from clustree and save as csv
graph <- clustree(allSeurat, prefix = "RNA_snn_res.", return = "graph")
edges <- graph %>%
  activate("edges") %>%
  as.data.frame()
colnames(edges)
write.csv(graph,'/Users/4472241/scCode/setNewClustering+keepOldUMAP/graphClustreeData.csv')

#Now re-run UMAP, if desired (looks almost the same)
allSeurat <- RunUMAP(allSeurat, dims = 1:100)

#Plot new umap with new clustering, notice umap doesn't change much
pdf('/Users/4472241/scCode/setNewClustering+keepOldUMAP/newUMAP_newClustering_all46_res.0.05.pdf', width = 9, height = 7)
DimPlot(allSeurat, reduction = 'umap', group.by = 'RNA_snn_res.0.05')
dev.off()

#Save seurat object as rds with new UMAP
saveRDS(allSeurat, '/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringNewUMAP_04-30-2021_all_46_samples_noHarmony.rds')

