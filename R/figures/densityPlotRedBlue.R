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
filedir <- '/Users/4472241/scCode/density_over-under_rep1/'

#Import the reference and convert to geyser 1 dataframe
ref_seurat <- readRDS('/Users/4472241/scCode/palantirEmbeddings/mapToRep1/reps123Integrated.rds')
geyser1 <- data.frame('duration' = ref_seurat@meta.data[["Palantir Embeddings"]][["X"]], 
                      'waiting' = ref_seurat@meta.data[["Palantir Embeddings"]][["Y"]])

rdsFileList <- c('CMML1_LTB3966_2020-10-15', 'CMML2_LTB4121_2020-10-15', 'CMML3_LTB5109_2020-10-15', 'CMML4_4-J-003_postBMT_2020-10-15',
                 'CMML5_4-K-001_HMA_2020-10-15', 'CMML6_4-Q-001_HMA_2020-10-15', 'CMML7_5-E-001_HMA_2020-10-15',
                 'CMML8_5-H-001_2020-10-15', 'CMML9_SF-100109-106293_2020-10-15', 'CMML10_SF-100109-111451_2020-10-15',
                 'CMML11_SF-100109-110236_2020-10-15', 'CMML12_SF-140401-00158_2020-10-15', 'CMML13_SF-140602-00025_2020-10-15',
                 'CMML14_SF-120628-00475_2020-10-15', 'CMML15_SF-140722-00012_2020-10-15', 'CMML16_SF-130612-00056_2020-10-15',
                 'CMML17_4-S-001_HMA_2020-10-15','CMML18_2-V-001_2020-10-15', 'CMML19_SF-141010-00049_2020-10-15',
                 'CMML20_SF-161129-00158_2020-10-15', 'CMML21_6-AE-001_2020-10-15', 'CMML22_6-AC-001_2020-10-15',
                 'CMML24_SF-100109-101914_2020-10-15','CMML25_SF-120425-00035_2020-10-15', 'CMML23+26_6-AD-001_2020-10-15',
                 'CMML27_SF-120926-00014_2020-10-15','CMML28_SF-140318-00065_2020-10-15', 'CMML29_SF-140507-00419_2020-10-15',
                 'CMML30_SF-160268-00045_2020-10-15', 'CMML31_SF-160722-00003_2020-10-15', 'CMML32_SF-161123-00029_2020-10-15',
                 'CMML33_SF-130328-00016_2020-10-15', 'CMML34_SF-141104-00108_2020-10-15', 'CMML35_SF-141114-00033_2020-10-15',
                 'CMML36_SF-140925-00135_2020-10-15', 'CMML37_SF-140613-00036_2020-10-15', 'CMML38_SF-140804-00065_2020-10-15',
                 'CMML39_SF-150102-00008_2020-10-15', 'CMML40_SF-130709-00171_2020-10-15')

for (i in 1:39){
  
  #Set sample ID to read in the right file
  sampleID <- paste0('cmml',i)
  
  #Read in the csv file containing the embeddings and convert to dataframe geyser2
  query <- read.csv(paste0('/Users/4472241/scCode/palantirEmbeddings/mapToRep1/',sampleID,'PalantirEmbeddings.csv'))
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
  pdf(paste0(filedir,sampleID,'_densityPlotRedBlueNEW.pdf'),width = 8, height = 6)
  print(ggplot()+
    geom_tile(diff12.m, mapping = aes(Duration, Waiting, fill=DensityDiff)) +
    stat_contour(aes(colour=..level..), binwidth=0.001) +
    scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
    coord_cartesian(xlim=xrng, ylim=yrng)+
    guides(colour=FALSE)+
    geom_point(geyser1, mapping = aes(x=duration, y=waiting), alpha = .1, size = .0005)+
    theme(panel.background = element_rect(fill = "white")))
  dev.off()
}

