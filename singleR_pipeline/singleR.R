library(SingleR)
browseVignettes("SingleR") #To see html page with examples/info about reference data
library(Seurat)
library(dplyr)

#Use Novershtern as the dataset (flow-sorted bulk RNAseq dataset)
hem.ref <- NovershternHematopoieticData()
?NovershternHematopoieticData
hem.ref

#Set dir to read in QC'ed data
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/'
#Read in seurat object with all samples
allSeurat <- readRDS(paste0(dir, 'postQC_CMML39+healthy8_05-06-2021.rds'))
allSeurat <- NormalizeData(allSeurat)

#Test to make sure it works
test.data <- allSeurat@assays[['RNA']]@data[,1:1000]
pred.hesc <- SingleR(test = test.data, ref = hem.ref, labels = hem.ref$label.main)
pred.hesc

#Run for each sample independently and make into df (df.all)
namesVec <- levels(allSeurat)
for (name in namesVec){
  sample <- allSeurat@assays[['RNA']]@data[,allSeurat@meta.data[["orig.ident"]] %in% name]
  pred.sample <- SingleR(test = sample, ref = hem.ref, labels = hem.ref$label.fine)
  df.row <- data.frame(rbind(table(pred.sample@listData[["first.labels"]])/
                    length(pred.sample@listData[["first.labels"]])))
  if (name %in% namesVec[1]){
    df.all <- df.row
  }else{
    df.all <- dplyr::bind_rows(df.all, df.row)
  }
}

row.names(df.all) <- namesVec
df.all[is.na(df.all)] <- 0
write.csv(df.all, '/Users/4472241/scCode/df_all_cellType.csv')
