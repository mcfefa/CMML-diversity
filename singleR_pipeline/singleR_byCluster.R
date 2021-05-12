library(SingleR)
#browseVignettes("SingleR") #To see html page with examples/info about reference data
library(Seurat)
library(dplyr)

rm(list = ls())

#Look at datasets, Novershtern is most specific to human hematopoiesis
novershtern <- NovershternHematopoieticData()
hpca <- HumanPrimaryCellAtlasData()
blueprintEncode <- BlueprintEncodeData()

#Load the clustered object
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
seuratClustered <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Run singleR
data <- seuratClustered@assays[["RNA"]]@data
clusters <- seuratClustered@meta.data[["RNA_snn_res.0.05"]]
pred.cluster <- SingleR(test = data, method = "cluster", ref = novershtern, labels = novershtern$label.main, 
                        clusters = clusters)

