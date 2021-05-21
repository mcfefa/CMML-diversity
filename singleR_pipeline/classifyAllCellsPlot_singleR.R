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
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Run singleR on all cells
data <- allSeurat@assays[["RNA"]]@data
pred.all <- SingleR(test = data, ref = novershtern, labels = novershtern$label.fine)

#Convert output into df and write to csv
allCells.df <- data.frame("CellNames" = pred.all@rownames, "CellType" = pred.all@listData[["first.labels"]])
write.csv(allCells.df, '/Users/4472241/scCode/singleR/allCellOutputNovershtern.csv')


#Plot the umap by cell type
## Grouping Cells as CMML, HMA, or Normal
allSeurat$cellTypeIndividualRes <- pred.all@listData[["first.labels"]]

pdf('/Users/4472241/scCode/singleR/UMAP_res=0.05_standard_pipeline_withHarmony_coloredByCellTypeSingleCellResNovershtern.pdf', width = 19, height = 5)
DimPlot(allSeurat, reduction = "umap", group.by = "cellTypeIndividualRes", shuffle = T)
dev.off()


