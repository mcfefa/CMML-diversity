library(Seurat)
library(affyio)
library(SummarizedExperiment)
library(SingleR)
library(GEOquery)
library("biomaRt")
library(BiocParallel)
BiocManager::install("ArrayExpress")
library(ArrayExpress)

rm(list = ls())

#Set dir to read in QC'ed and plotted/clustered data
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
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
pred.sample <- SingleR(test = all.data, ref = gse.mat, labels = gse@colData@listData[["source_name_ch1"]], BPPARAM = MulticoreParam(workers=4))


table(pred.sample@listData[["first.labels"]])/length(pred.sample@listData[["first.labels"]])

allSeurat$cellTypeRapinIndividualRes <- pred.sample@listData[["first.labels"]]

#Plot, showing cell type
pdf('/Users/4472241/scCode/singleR/standard_pipeline_withHarmony_coloredByCellTypeSingleCellRapin.pdf', width = 8, height = 5)
DimPlot(allSeurat, reduction = "umap", group.by = "cellTypeRapinIndividualRes", 
        cols = c("dark grey", "cyan", "orange", "red", "pink", "green", "light grey", "yellow",
                 "purple", "black", "brown"), shuffle = T)
dev.off()

AEDownloadBulk <- function(accession, type = "processed", out = getwd()) {
  # create output dir
  dir.create(paste0(out,"/ArrayExpress"), showWarnings = F)
  setwd(paste0(out,"/ArrayExpress"))
  # run getAE() in loop for accession numbers
  for(i in accession) {
    dir.create(i, showWarnings = F)
    getAE(i, type = type, path = i, extract = T)
    zip <- list.files(i, full.names = T)[grep(".zip",list.files(i))]
    if(file.exists(zip[1])) file.remove(zip)
  }
  setwd("../")
}
accession = c("E-MTAB-5456")
x <- AEDownload(accession)


AEDownload <- function(accession, type = "processed", out = getwd(), import = T) {
  if(length(accession)>1){stop("length(accession) > 1. Please provide only a single accession number. ")}
  # create output dir
  dir.create(paste0(out,"/ArrayExpress"), showWarnings = F)
  setwd(paste0(out,"/ArrayExpress"))
  
  # run getAE() in loop for accession numbers
  i = accession
  dir.create(i, showWarnings = F)
  AE <- getAE(i, type = type, path = i, extract = T)
  zip <- list.files(i, full.names = T)[grep(".zip",list.files(i))]
  if(file.exists(zip[1])) file.remove(zip)
  
  # Import data in R list object
  if(import == T) {
    ls <- list()
    sdrf <- list.files(i, full.names = T)[grep(".sdrf.txt",list.files(i))]
    idf <- list.files(i, full.names = T)[grep(".idf.txt",list.files(i))]
    ls[["sdrf"]] <- read.delim(file = sdrf, row.names = 1)
    ls[["idf"]] <- read.delim(file = idf)
    for(k in AE$processedFiles) {
      ls[[k]] <- read.table(file = paste(i,k,sep = "/"),row.names = 1)
    }
    ls[["info"]] <- AE
  }
  setwd("../")
  if(import == T) return(ls)
}


