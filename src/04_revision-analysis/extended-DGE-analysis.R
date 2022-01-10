library(Seurat)

#Import the seurat object with all of our samples assigned to clusters
dir <- '/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/figures/R data (Seurat objects)/'
allSeurat <- readRDS(paste0(dir, 'allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021.rds'))

## check FcgR expression per Reviewer comment
feature_vec <- c("FCGR2A")
FeaturePlot(allSeurat, features = feature_vec)

## differential gene expression 

# set clustering as the default identity
allSeurat <- SetIdent(object=allSeurat, value = allSeurat$RNA_snn_res.0.05)
levels(allSeurat)

# differential gene expression
options(future.globals.maxSize = 16000 * 1024^5)
mrks <- FindMarkers(allSeurat, ident.1=2, ident.2=0)
write.csv(mrks, './src/04_revision-analysis/clus2-vs-clus0.de.markers.csv')
### ran out of memory on laptop, moved to HiPerGator

# sample names are stored in
unique(allSeurat@meta.data$orig.ident)
allSeurat <- SetIdent(object=allSeurat, value = allSeurat@meta.data$orig.ident)

outdir <- "./src/04_revision-analysis"

mrks.ptHMA.A <- FindMarkers(allSeurat, ident.1="LTB3966", ident.2="SF13032800016")
write.csv(mrks.ptHMA.A,paste(outdir,"/HMA-TX_A1_LTB3966-vs-SF13032800016.de.markers.csv",sep=""))

mrks.ptHMA.A2 <- FindMarkers(allSeurat, ident.1="LTB3966", ident.2="SF13070900171")
write.csv(mrks.ptHMA.A2,paste(outdir,"/HMA-TX_A2_LTB3966-vs-SF13070900171.de.markers.csv",sep=""))

mrks.ptHMA.B <- FindMarkers(allSeurat, ident.1="SF14060200025", ident.2="SF14092500135")
write.csv(mrks.ptHMA.B,paste(outdir,"/HMA-TX_B_SF14060200025-vs-SF14092500135.de.markers.csv",sep=""))

mrks.ptHMA.C <- FindMarkers(allSeurat, ident.1="6AD001", ident.2="SF14110400108")
write.csv(mrks.ptHMA.C,paste(outdir,"/HMA-TX_C_6AD001-vs-SF14110400108.de.markers.csv",sep=""))

mrks.ptHMA.D <- FindMarkers(allSeurat, ident.1="SF14031800065", ident.2="SF14080400065")
write.csv(mrks.ptHMA.D,paste(outdir,"/HMA-TX_D_SF14031800065-vs-SF14080400065.de.markers.csv",sep=""))

mrks.ptRux.E <- FindMarkers(allSeurat, ident.1="SF14050700419", ident.2="SF15010200008")
write.csv(mrks.ptRux.E,paste(outdir,"/Rux-TX_E_SF14050700419-vs-SF15010200008.de.markers.csv",sep=""))

mrks.ptRux.F <- FindMarkers(allSeurat, ident.1="SF14072200012", ident.2="SF14111400033")
write.csv(mrks.ptRux.F,paste(outdir,"/Rux-TX_F_SF14072200012-vs-SF14111400033.de.markers.csv",sep=""))

mrks.ptChemo.G <- FindMarkers(allSeurat, ident.1="SF14040100158", ident.2="SF14061300036")
write.csv(mrks.ptChemo.G,paste(outdir,"/Chemo-TX_G_SF14040100158-vs-SF14061300036.de.markers.csv",sep=""))












