library(Seurat)
BiocManager::install("affyio")
library(affyio)
library(SummarizedExperiment)
library(SingleR)
BiocManager::install("GEOquery")
library(GEOquery)
BiocManager::install("biomaRt")
library("biomaRt")
library(BiocParallel)
BiocManager::install("ArrayExpress")
library(ArrayExpress)
library(celldex)
library(ggplot2)

rm(list = ls())

#Set dir to read in QC'ed and plotted/clustered data
dir <- '/Users/Brian/Downloads/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

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

#Run singleR using this dataset as a reference
all.data <- allSeurat@assays[['RNA']]@data
#Set to first 16 bulk samples to only look at the CD34+ progenitors in the dataset
pred.sample <- SingleR(test = all.data, ref = gse.mat[,1:16], labels = gse@colData@listData[["source_name_ch1"]][1:16], BPPARAM = MulticoreParam(workers=4))

#Save output as csv so we can re-load the cell type assignment at any time
out.csv <- data.frame("Cell Name" = pred.sample@rownames, "Cell Type" = pred.sample@listData[["first.labels"]])
write.csv(out.csv, '/Users/Brian/scCode/singleR/allCells_cellType_RapinProgOnly.csv', row.names = F)
