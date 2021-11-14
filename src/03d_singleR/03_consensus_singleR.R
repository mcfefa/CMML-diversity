#Load all of the cell output files from each ref for prog only
library(SingleR)
#browseVignettes("SingleR") #To see html page with examples/info about reference data
library(Seurat)
library(dplyr)
library(celldex)
library(SummarizedExperiment)
library(ggplot2)
library(plyr)


rm(list = ls())

novershtern <- NovershternHematopoieticData()

classified_nov <- read.csv('/Users/Brian/scCode/singleR/allCellOutput_nov_ProgOnly_08-03-2021.csv')
classified_hpca <- read.csv('/Users/Brian/scCode/singleR/allCellOutput_hpca_ProgOnly_08-03-2021.csv')
classified_blu <- read.csv('/Users/Brian/scCode/singleR/allCellOutput_blueprintEncode_ProgOnly_08-03-2021.csv')
classified_rap <- read.csv('/Users/Brian/scCode/singleR/allCells_cellType_RapinProgOnly.csv')

unique(classified_nov$CellType)
classified_nov$CellType <- gsub("s", "", classified_nov$CellType)
classified_nov$CellType <- gsub("B cell", "CLP", classified_nov$CellType)
classified_nov$CellType <- gsub("CD8\\+ T cell", "CLP", classified_nov$CellType)
classified_nov$CellType <- gsub("CD4\\+ T cell", "CLP", classified_nov$CellType)
classified_nov$CellType <- gsub("NK T cell", "CLP", classified_nov$CellType)
unique(classified_nov$CellType)

unique(classified_hpca$CellType)
classified_hpca$CellType <- gsub("HSC_CD34\\+", "HSC", classified_hpca$CellType)
classified_hpca$CellType <- gsub("Pro-B\\_cell\\_CD34\\+", "CLP", classified_hpca$CellType)
classified_hpca$CellType <- gsub("T\\_cells", "CLP", classified_hpca$CellType)
classified_hpca$CellType <- gsub("Erythroblast", "MEP", classified_hpca$CellType)
classified_hpca$CellType <- gsub("Pro-Myelocyte", "GMP", classified_hpca$CellType)
classified_hpca$CellType <- gsub("NK\\_cell", "CLP", classified_hpca$CellType)
unique(classified_hpca$CellType)

unique(classified_blu$CellType)
classified_blu$CellType <- gsub("naive B-cells", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("Memory B-cells", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD8\\+ Tem", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD8\\+ T-cells", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD8\\+ Tcm", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD4\\+ Tem", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD4\\+ Tcm", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("Tregs", "CLP", classified_blu$CellType)
classified_blu$CellType <- gsub("CD4\\+ T-cells", "CLP", classified_blu$CellType)
unique(classified_blu$CellType)

unique(classified_rap$Cell.Type)
classified_rap$Cell.Type <- gsub("Granulocyte Macrophage Progenitors", "GMP", classified_rap$Cell.Type)
classified_rap$Cell.Type <- gsub("Common myeloid progenitors", "CMP", classified_rap$Cell.Type)
classified_rap$Cell.Type <- gsub("Megakaryocyte-Erythroid Progenitors", "MEP", classified_rap$Cell.Type)
classified_rap$Cell.Type <- gsub("Hematopoietic Stem Cells", "HSC", classified_rap$Cell.Type)
classified_rap$Cell.Type <- gsub("Multipotent Progenitors", "MPP", classified_rap$Cell.Type)
unique(classified_rap$Cell.Type)

df_combined <- data.frame("CellName" = classified_nov$CellNames, "cellType_nov" = classified_nov$CellType, 
                          "cellType_hpca" = classified_hpca$CellType, "cellType_blu" = classified_blu$CellType,
                          "cellType_rap" = classified_rap$Cell.Type)
#Define function used next
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

out <- c()
for (i in 1:length(df_combined[,1])){
  #Blueprint encode has only CLP, so trust it for lymphoid cells
  #if(df_combined$cellType_blu[i] == "CLP"){
  #  out[i] <- "CLP"
  #  next
  #}
  
  #Other cells, get mode
  vec <- as.character(df_combined[i,2:5])
  uniqv <- unique(vec)
  if (length(uniqv) == 1){
    out[i] <- uniqv[1] 
  }else if (length(uniqv) == 2){
    modes <- unique(as.character(Modes(vec)))
    if (length(modes) == 2){
      out[i] <- "No Consensus"
    }else{
      out[i] <- modes
    }
  }else if (length(uniqv) == 3) {
    out[i] <- as.character(Modes(vec))
  }else if (length(uniqv) == 4){
    out[i] <- "No Consensus"
  }
}

#Check output
table(out)
out.df <- data.frame("CellName" = df_combined$CellName, "ConsensusType" = out)
write.csv(out.df, '/Users/Brian/scCode/singleR/singleR_consensus_08-15-2021.csv')

#Set dir to read in QC'ed and plotted/clustered data
dir <- '/Users/Brian/Downloads/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))
allSeurat$consensusSingleR <- out
saveRDS(allSeurat, '/Users/Brian/scCode/singleR/seuratObject_allSamples_withsingleRType.rds')

#Subset to cluster 2
clus2 <- subset(allSeurat, subset = clusterResolution_0.05 == 2)
table(clus2$consensusSingleR)/dim(clus2)[2]

pdf('/Users/Brian/scCode/singleR/umap_singleR_consensus_08-22-2021.pdf', width = 8, height = 6)
DimPlot(allSeurat, group.by = "consensusSingleR", 
        cols = c("cyan", "black", "blue", "pink", "red", "green", "grey"),
        raster = F, shuffle = T, pt.size = 0.01)
dev.off()

#Get outfile with cell type fraction for each sample:
df <- data.frame("sample" = allSeurat$orig.ident, "cellType" = allSeurat$consensusSingleR)
split.df <- split(df, df$sample)

for ( i in 1:length(split.df)){
  this_step <- as.data.frame(table(split.df[[i]]))$Freq/sum(as.data.frame(table(split.df[[i]]))$Freq)
  this_step.df <- as.data.frame(t(as.data.frame(this_step)))
  colnames(this_step.df) <- as.data.frame(table(split.df[[i]]))$cellType
  row.names(this_step.df) <- unique(as.data.frame(table(split.df[[i]]))$sample)
  if (i > 1){
    all.df <- rbind.fill(all.df, this_step.df)
  }else{
    all.df <- as.data.frame(this_step.df)
  }
}

all.df[is.na(all.df)] <- 0
row.names(all.df) <- names(split.df)

#Write HSC number for each sample to outfile
write.csv(all.df, '/Users/Brian/scCode/singleR/singleR_consensus_sampleSummary_08-22-2021.csv')

#Representative case
table(allSeurat$orig.ident)
# Subset on a value in the object meta data
cmml15 <- subset(allSeurat, subset = orig.ident == "SF14072200012")
hua3 <- subset(allSeurat, subset = orig.ident == "HuaPt3")


DimPlot(cmml15, group.by = "consensusSingleR", 
        cols = c("cyan", "black", "blue", "pink", "red", "dark blue", "grey"),
        raster = F, shuffle = T)
DimPlot(hua3, group.by = "consensusSingleR", 
        cols = c("cyan", "black", "blue", "pink", "red", "dark blue", "grey"),
        raster = F, shuffle = T)

seurat_rep_HSC_depletion <- hua3 <- subset(allSeurat, subset = orig.ident %in% c("HuaPt3", "SF14072200012"))

DimPlot(seurat_rep_HSC_depletion, group.by = "consensusSingleR", split.by = "orig.ident",
        cols = c("cyan", "black", "blue", "pink", "red", "grey"),
        raster = F, shuffle = T)
