library(reshape2)# For melt function
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library("MASS")

#Clear workspace
rm(list = ls())

#Set dir to save to
filedir <- '/Users/Brian/scCode/palantirProjection_11012021/'

#Import the reference and convert to geyser 1 dataframe
ref_seurat <- readRDS(paste0(filedir, 'reps123Integrated.rds'))
geyser1 <- data.frame('duration' = ref_seurat@meta.data[["Palantir Embeddings"]][["X"]], 
                      'waiting' = ref_seurat@meta.data[["Palantir Embeddings"]][["Y"]])

rdsFileList <- c('CMML1', 'CMML2', 'CMML3', 'CMML4_postBMT',
                 'CMML5_HMA', 'CMML6_HMA', 'CMML7_HMA','CMML8', 'CMML9', 'CMML10',
                 'CMML11', 'CMML12', 'CMML13','CMML14', 'CMML15', 'CMML16',
                 'CMML17','CMML18', 'CMML19','CMML20', 'CMML21', 'CMML22',
                 'CMML24','CMML25', 'CMML23+26','CMML27','CMML28', 'CMML29',
                 'CMML30', 'CMML31', 'CMML32','CMML33', 'CMML34', 'CMML35',
                 'CMML36', 'CMML37', 'CMML38','CMML39', 'CMML40')

for (i in 1:39){
  #Set sample ID to read in the right file
  sampleID <- paste0('cmml',i)
  
  #Read in the csv file containing the embeddings and convert to dataframe geyser2
  query <- read.csv(paste0(filedir,sampleID,'PalantirEmbeddings.csv'))
  geyser2 <- data.frame('duration' = query$X1, 'waiting' = query$X2)
  
  #Set manually the common x and y range for geyser1 and geyser2
  xrng = c(-32, 28)
  yrng = c(-31, 28)
  
  # Calculate the 2d density estimate over the common range
  d1 = kde2d(geyser1$duration, geyser1$waiting, lims=c(xrng, yrng), n=200, h=3)
  d2 = kde2d(geyser2$duration, geyser2$waiting, lims=c(xrng, yrng), n=200, h=3)
  
  # Confirm that the grid points for each density estimate are identical
  identical(d1$x, d2$x) # TRUE
  identical(d1$y, d2$y) # TRUE
  
  # Calculate the difference between the 2d density estimates
  diff12 = d1 
  diff12$z = log(d2$z*1000+1) - log(d1$z*1000+1)
  
  # First, add row and column names (x and y grid values) to the z-value matrix
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  
  # Now melt it to long format
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  names(diff12.m) = c("Duration","Waiting","DensityDiff")
  
  # Plot and save difference between geyser2 and geyser1 density
  pdf(paste0(filedir,sampleID,'_densityPlotRedBlue_11012021.pdf'),width = 8, height = 6)
  print(ggplot()+
          geom_tile(diff12.m, mapping = aes(Duration, Waiting, fill=DensityDiff)) +
          stat_contour(aes(colour=..level..), binwidth=0.001) +
          scale_fill_gradient2(low="blue",mid="white", high="red", midpoint=0) +
          coord_cartesian(xlim=xrng, ylim=yrng)+
          guides(colour=FALSE)+
          geom_point(geyser1, mapping = aes(x=duration, y=waiting), alpha = .1, size = .0005)+
          theme(panel.background = element_rect(fill = "white")))
  dev.off()
}
