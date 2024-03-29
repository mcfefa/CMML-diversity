library(SingleR)
#browseVignettes("SingleR") #To see html page with examples/info about reference data
library(Seurat)
library(dplyr)
library(celldex)

rm(list = ls())

#Look at datasets, Novershtern is most specific to human hematopoiesis
novershtern <- NovershternHematopoieticData()
hpca <- HumanPrimaryCellAtlasData()
blueprintEncode <- BlueprintEncodeData()

#Load the clustered object
dir <- '/Users/Brian/Downloads/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

############# HPCA #############
#Run singleR on all cells, restricting the reference to progenitors
labels <- hpca$label.main[c(36:59, 66:73, 116:121, 227:230, 308:320, 479:480, 516:534, 
                          625:634, 637:638, 643:646, 693:697)]
ref <- hpca@assays@data@listData[["logcounts"]][,c(36:59, 66:73, 116:121, 227:230, 308:320, 479:480, 516:534, 
                                                   625:634, 637:638, 643:646, 693:697)]
data <- allSeurat@assays[["RNA"]]@data
pred.all.hpca <- SingleR(test = data, ref = ref, labels = labels)
#Write to csv
allCells.hpca <- data.frame("CellNames" = pred.all.hpca@rownames, "CellType" = pred.all.hpca@listData[["first.labels"]])
write.csv(allCells.hpca, '/Users/Brian/scCode/singleR/allCellOutput_hpca_ProgOnly_08-03-2021.csv')

########## Novershtern #########
labels <- as.vector(novershtern$label.main[c(7:30, 79:82, 96:109, 122:130, 154:211)])
ref <- novershtern@assays@data@listData[['logcounts']][,c(7:30, 79:82, 96:109, 122:130, 154:211)]
data <- allSeurat@assays[["RNA"]]@data
pred.all.nov <- SingleR(test = data, ref = ref, labels = labels)
#Write to csv
allCells.nov <- data.frame("CellNames" = pred.all.nov@rownames, "CellType" = pred.all.nov@listData[["first.labels"]])
write.csv(allCells.nov, '/Users/Brian/scCode/singleR/allCellOutput_nov_ProgOnly_08-03-2021.csv')

########### Blueprint Encode ##########
labels <- as.vector(blueprintEncode$label.fine[c(4, 12:17, 20:21, 26:29, 32, 35, 
                                                 43:54, 57:59, 62:64, 66:67, 72:73, 77, 
                                                 79:82, 84:92, 102, 139:140, 158)])
ref <- blueprintEncode@assays@data@listData[['logcounts']][,c(4, 12:17, 20:21, 26:29, 32, 35, 
                                                              43:54, 57:59, 62:64, 66:67, 72:73, 77, 
                                                              79:82, 84:92, 102, 139:140, 158)]
data <- allSeurat@assays[["RNA"]]@data
pred.all.blu <- SingleR(test = data, ref = ref, labels = labels)
#Write to csv
allCells.blu <- data.frame("CellNames" = pred.all.blu@rownames, "CellType" = pred.all.blu@listData[["first.labels"]])
write.csv(allCells.blu, '/Users/Brian/scCode/singleR/allCellOutput_blueprintEncode_ProgOnly_08-03-2021.csv')
