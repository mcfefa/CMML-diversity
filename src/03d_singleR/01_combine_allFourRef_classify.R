library(SingleR)
#browseVignettes("SingleR") #To see html page with examples/info about reference data
library(Seurat)
library(dplyr)
library(affyio)
library(SummarizedExperiment)
library(GEOquery)
library(biomaRt)
library(BiocParallel)
library(ArrayExpress)

rm(list = ls())

#Get dataset from GEO accession number
gse = getGEO("GSE42519")[[1]]

#Convert to "SummarizedExperiment" format
gse = as(gse, "SummarizedExperiment")

#Get expression matrix
gse.mat <- gse@assays@data@listData[["exprs"]]

#Use biomaRt to convert from Affy to HGNC symbol to match our seurat object
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl",mart=mart)
genes_affy <- row.names(gse.mat)
#Run the getBM function to convert
genes_symbol <- getBM(c("affy_hg_u133_plus_2", "hgnc_symbol"), c("affy_hg_u133_plus_2"), genes_affy, mart)

# Now go from matched affy and hgnc to hgnc
outNames <- c()
for (i in 1:length(row.names(gse.mat))){
  newName <- genes_symbol[genes_symbol$affy_hg_u133_plus_2 %in% row.names(gse.mat)[i],"hgnc_symbol"]
  if(length(newName)>0){
    outNames[i] <- newName
  }else{
    outNames[i] <- ""
  }
}
#Set row names to be the HGNC symbols associated with the affy ID
row.names(gse.mat) <- outNames

#Look at datasets, Novershtern is most specific to human hematopoiesis
novershtern <- NovershternHematopoieticData()
hpca <- HumanPrimaryCellAtlasData()
blueprintEncode <- BlueprintEncodeData()

#Load the clustered object
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Run singleR on all cells, restricting the reference to progenitors and combining ref's
labels_hpca <- hpca$label.main[c(231:235, 257:284, 302:307, 416:431, 592:599, 608:616, 625:634, 637:638, 676:683)]
labels_hpca <- labels_hpca[64:84]
labels_hpca <- c(rep("MSC", 9), rep("HSC", 6), labels_hpca[16:21])
ref_hpca <- hpca@assays@data@listData[["logcounts"]][,c(159:168, 231:235, 257:284, 302:307, 416:431, 592:599, 608:616, 625:634, 637:638, 676:683)]
ref_hpca <- ref_hpca[,64:84]
labels_nov <- as.vector(novershtern$label.main[c(27:30, 79:82, 96:109, 122:130)])
labels_nov <- gsub("s", "", labels_nov)
ref_nov <- novershtern@assays@data@listData[['logcounts']][,c(27:30, 79:82, 96:109, 122:130)]
labels_blue <- as.vector(blueprintEncode$label.fine[c(4, 26:29, 43:54, 58, 62:64, 66:67, 72, 80:82, 
                                                 87:91, 158)])
ref_blue <- blueprintEncode@assays@data@listData[['logcounts']][,c(4, 26:29, 43:54, 58, 62:64, 66:67, 72, 80:82, 
                                                                   87:91, 158)]
labels_rapin <- gse@colData@listData[["source_name_ch1"]][1:16]
labels_rapin <- c(rep("HSC", 4), rep("MPP", 2), rep("CMP", 3), rep("GMP", 5), rep("MEP", 2))
ref_rapin <- gse.mat[,1:16]
data <- allSeurat@assays[["RNA"]]@data
pred.all <- SingleR(test = data, ref = list(ref_hpca, ref_nov, ref_blue, ref_rapin),
                    labels = list(labels_hpca, labels_nov, labels_blue, labels_rapin))

#Convert output into df and write to csv
allCells.df <- data.frame("CellNames" = pred.all@rownames, "CellType" = pred.all@listData[["first.labels"]])
write.csv(allCells.df, '/Users/4472241/scCode/singleR/allCellOutput_allFourRefs_ProgOnly.csv')


#Plot the umap by cell type
## Grouping Cells as CMML, HMA, or Normal
allSeurat$cellTypeIndividualRes <- pred.all@listData[["first.labels"]]

pdf('/Users/4472241/scCode/singleR/UMAP_res=0.05_standard_pipeline_withHarmony_coloredByCellTypeSingleCellRes_allFourRefs_ProgOnly_noHSCGCSF.pdf', 
    width = 8, height = 5)
DimPlot(allSeurat, reduction = "umap", group.by = "cellTypeIndividualRes", shuffle = T)
dev.off()

