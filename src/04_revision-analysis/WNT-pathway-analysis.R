library(Seurat)

#Import the seurat object with all of our samples assigned to clusters
datadir <- '/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/figures/R data (Seurat objects)/'
allSeurat <- readRDS(paste0(datadir, 'allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021.rds'))

dir <- '~/GitHub/CMML-diversity/src/04_revision-analysis/'

# read in GSEA Geneset
wntsig <- read.delim(paste0(dir,"GSEA_KEGG-WNT-SIGNALING-Geneset.txt",sep=""),header=TRUE,sep="\t",dec=",",skip=1)

for (j in 1:length(orig.ident.vec)){
  
  cmml <- subset(allSeurat, orig.ident == orig.ident.vec[j], features = allGenes)
  
  #Make into dataframes (separate by gene list then combine)
  genesKeep_VG <- genes_hscVG[genes_hscVG %in% cmml@assays[['RNA']]@data@Dimnames[[1]]]
  cmmlVG.df <- data.frame(cmml@assays[["RNA"]]@data[genesKeep_VG,])

  cmml.df <- do.call("rbind", list(cmmlWu.df, cmmlVG.df, cmmlEp.df))

  #Downsample the number of cells for equal comparison/easier plotting
  n_samplesCMML <- 858
  n_samplesNormal <- 1000
  cmmlDS <- cmml.df[,sample(ncol(cmml.df), n_samplesCMML)]
  normalsDS <- normals.df[,sample(ncol(normals.df), n_samplesNormal)]
  dfAll <- cbind(cmmlDS, normalsDS)
  
  #Scale by z-score of the included samples, set max at +/- 3 for color scale
  minV <- -3
  maxV <- 3
  for (i in 1:dim(dfAll)[[1]]){
    dfAll[i,] <- (dfAll[i,]-rowMeans(dfAll)[i])/sd(dfAll[i,])
    dfAll[i,] <- sapply(dfAll[i,], function(y) min(max(y,minV),maxV))
  }
  
  #convert to matrix
  matrixNormal <- as.matrix(dfAll)
  
  #Find where to draw lines based on last gene in Van Galen and Wu
  genesKeep <- allGenes[allGenes %in% row.names(cmml.df)]
  colBreak1 <- which(allGenes %in% lastWu)
  colBreak2 <- which(allGenes %in% lastVG)
  
  #Define colormap 
  coul <- colorRampPalette(c('blue','white', 'red'))
  
  #Plot and save
  png(paste0(dir,'Wnt-signaling-score.png'))#, width = 8, height = 6)
  heatmap.2(t(matrixNormal), rowsep = n_samplesCMML,
            colsep = c(colBreak1, colBreak2),
            sepwidth=c(.3,3),
            sepcolor="black", trace="none",
            Rowv=F,Colv=F, scale="none", dendrogram="none",key=F, 
            lhei = c(0.01,.5),margins=c(1,8), col = coul)
  dev.off()
}




