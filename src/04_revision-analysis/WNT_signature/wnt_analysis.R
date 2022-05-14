# WNT pathway figure
library(ggplot2)
library(Seurat)
library(gridExtra)

rm(list = ls())

allSeurat <- readRDS("~/scCode/allSeurat_withsingleR_withBias_03092022.rds")

# Restrict to cluster 0 and 2 GMPs
GMP_Seurat1 <- subset(allSeurat, subset = consensusSingleR == "GMP")
GMP_Seurat <- subset(GMP_Seurat1, subset = clusterResolution_0.05 %in% c("0", "2"))

wnt <- read.table("~/scCode/wnt_signatures/GSEA_KEGG-WNT-SIGNALING-Geneset.txt", sep = "\t", skip = 2)$V1

# Add module score to seurat object
wnt_list <- list(wnt1)
GMP_Seurat <- AddModuleScore(GMP_Seurat, wnt_list, name = "WNT_signature", search = F)
# SKP1P2 not found, continue without (could not find alternative name that was in our data)

pdf("~/scCode/wnt_signatures/wnt_plot.pdf")
FeaturePlot(GMP_Seurat, features = c("WNT_signature1"), order = T, min.cutoff = 0.0) + NoAxes()
dev.off()

# Make feature plots of CTNNB1, IRF8, and the WNT score we just added
p1 <- FeaturePlot(GMP_Seurat, features = c("CTNNB1"), order = T, min.cutoff = 0.0) + NoAxes()
p2 <- FeaturePlot(GMP_Seurat, features = c("IRF8"), order = T, min.cutoff = 0.0) + NoAxes()
p3 <- FeaturePlot(GMP_Seurat, features = c("WNT_signature1"), order = T, min.cutoff = 0.0) + NoAxes() + ggtitle("WNT Signature")

# Plot feature plots together
pdf("~/scCode/wnt_signatures/supp_plot.pdf", width = 12/1.5, height = 5/1.5)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()


# Make violin plot comparing WNT across clus2 and not clus2
pdf("~/scCode/vlnPlot_WNT.pdf", width = 4, height = 5)
VlnPlot(GMP_Seurat, features = "WNT_signature1", group.by = "clusterResolution_0.05", pt.size = 0)+
  ggtitle("WNT signature score in GMPs") + scale_fill_manual(values = c("black", "#31CEFF")) + xlab("Cluster") +
  guides(fill=guide_legend(title="Cluster"))+
  theme(legend.text = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold")) + NoLegend()
dev.off()




