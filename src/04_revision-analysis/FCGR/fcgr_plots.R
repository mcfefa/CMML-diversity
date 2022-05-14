# FCGR
library(ggplot2)
library(Seurat)

rm(list = ls())

allSeurat <- readRDS("~/scCode/allSeurat_withsingleR_withBias_03092022.rds")

p <- FeaturePlot(allSeurat, features = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A"), combine = F)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()
}

pdf("~/scCode/FCGR_feature_plots_for_supp.pdf")
cowplot::plot_grid(plotlist = p)
dev.off()

