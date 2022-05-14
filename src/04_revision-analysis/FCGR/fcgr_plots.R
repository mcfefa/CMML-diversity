# FCGR

library(Seurat)

rm(list = ls())

allSeurat <- readRDS("~/scCode/allSeurat_withsingleR_withBias_03092022.rds")

pdf("~/scCode/FCGR_feature_plots_for_supp.pdf")
FeaturePlot(allSeurat, features = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A"), NoAxes)
dev.off()


p <- FeaturePlot(allSeurat, features = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A"), combine = F)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()
}

pdf("~/scCode/FCGR_feature_plots_for_supp.pdf")
cowplot::plot_grid(plotlist = p)
dev.off()

#VlnPlot(allSeurat, features = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A"), group.by = "bias")

# Make a signature score and compare across GMPs (like in WNT)
GMP_Seurat1 <- subset(allSeurat, subset = consensusSingleR == "GMP")
GMP_Seurat <- subset(GMP_Seurat1, subset = clusterResolution_0.05 %in% c("0", "2"))

wnt1 <- read.table("~/scCode/wnt_signatures/wnt1_GSEA_KEGG-WNT-SIGNALING-Geneset.txt", sep = "\t", skip = 2)$V1
wnt2 <- read.csv("~/scCode/wnt_signatures/wnt2.txt", header = F)[1,]

wnt_1_2 <- list(wnt1, wnt2)
GMP_Seurat <- AddModuleScore(GMP_Seurat, list(c(c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B"))), name = "FCGR_signature", search = F)

pdf("~/scCode/vlnPlot_FCGR.pdf", width = 5, height = 5)
VlnPlot(GMP_Seurat, features = "FCGR_signature1", group.by = "clusterResolution_0.05", pt.size = 0)+
  ggtitle("Fc gamma receptor score in GMPs") + scale_fill_manual(values = c("black", "#31CEFF")) + xlab("Cluster") +
  guides(fill=guide_legend(title="Cluster"))+
  theme(legend.text = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold")) + NoLegend()
dev.off()
