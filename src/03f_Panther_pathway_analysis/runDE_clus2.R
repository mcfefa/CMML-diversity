library(Seurat)

#Run differential expression on allSeurat object for cluster 2 and export to csv

rm(list = ls())

#Import the seurat object with all of our samples assigned to clusters
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Map values to a separate cluster 2 metadata entry
allSeurat$clus2 <- plyr::mapvalues(
  x = allSeurat$RNA_snn_res.0.05, 
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
  to = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

#Run differential expression. group.by specifies which metadata, ident.1 specifies the group within that metadata.
#In this way, we run DE between cluster 2 and all other cells (ident.2 blank means run against all other cells). 
clus2.de.markers <- FindMarkers(allSeurat, group.by = "clus2", ident.1 = 1)

write.csv(clus2.de.markers, '/Users/4472241/scCode/cluster2_DE/clus2.de.markers.csv')
