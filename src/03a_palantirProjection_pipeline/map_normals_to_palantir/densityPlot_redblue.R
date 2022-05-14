densityPlot_redblue <- function(query_embeddings, reference_embeddings, sampleID, outdir) {
  library(reshape2)# For melt function
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(BiocNeighbors)
  library("MASS")
  
  #sampleID <- paste0("CMML_", cmml_number)
  
  geyser1 <- data.frame('duration' = reference_embeddings[,1], 
                        'waiting' = reference_embeddings[,2])
  geyser2 <- data.frame('duration' = query_embeddings[,1], 
                        'waiting' = query_embeddings[,2])
  
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
  
  pdf(paste0(outdir,sampleID,'_densityPlotRedBlue_palantir_03142022.pdf'),width = 8, height = 6)
  print(ggplot()+
          geom_tile(diff12.m, mapping = aes(Duration, Waiting, fill=DensityDiff)) +
          stat_contour(aes(colour=..level..), binwidth=0.05) +
          scale_fill_gradient2(low="blue",mid="white", high="red", midpoint=0) +
          coord_cartesian(xlim=xrng, ylim=yrng)+
          guides(colour=FALSE)+xlab("Pal_TSNE_1")+ylab("Pal_TSNE_2")+
          geom_point(geyser2, mapping = aes(x=duration, y=waiting), alpha = .1, size = .0005)+
          geom_point(geyser1, mapping = aes(x=duration, y=waiting), alpha = .1, size = .0005)+
          theme(panel.background = element_rect(fill = "white")))
  dev.off()
}
