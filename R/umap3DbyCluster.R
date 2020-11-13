#Make a nice 3D UMAP
library(Seurat)
library(plotly)

#Clear workspace variables
rm(list=ls())

#Set file directory to save to
filedir <- "/Users/4472241/scCode/layerUMAP+ClusteringNoInit/"

#Import the object you want to plot
integratedObject <- readRDS('/Users/4472241/scCode/layerUMAP+ClusteringNoInit/31+7_samples_umap_layered_excludeInitial_with3DUMAP.rds')

# Re-run UMAPs if there's no 3D UMAP reduction in the object
#integratedObject <- RunUMAP(integratedObject, dims = 1:50, n.components = 3L, reduction = "harmony",
                            #reduction.key = "UMAP3D_", reduction.name = "umap3D")

# Extract UMAPinformation from Seurat Object
umap_1 <- integratedObject[["umap3D"]]@cell.embeddings[,1]
umap_2 <- integratedObject[["umap3D"]]@cell.embeddings[,2]
umap_3 <- integratedObject[["umap3D"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = integratedObject, reduction = "umap3D")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = integratedObject, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3", "clusters_noHarmony_res.0.05_30Neighbors"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))



# Plot your data, in this example my Seurat object had 21 clusters (0-20)
p <- plot_ly(data = plot.data, 
             x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, 
             color = ~clusters_noHarmony_res.0.05_30Neighbors, 
             #Optionally set your own color scheme
             #colors = c("lightseagreen",
             #          "gray50",
             #         "darkgreen",
             #        "red4",
             #       "red",
             #      "turquoise4",
             #     "black",
             #    "yellow4",
             #   "royalblue1",
             #  "lightcyan3",
             # "peachpuff3",
             #"khaki3",
             #"gray20",
             #"orange2",
             #"royalblue4",
             #"yellow3",
             #"gray80"),
             #"darkorchid1",
             #"lawngreen",
             #"plum2",
             #"darkmagenta"),
             type = "scatter3d",
             mode = "markers", 
             marker = list(size = 1, width=0.5), # controls size of points
             text=~label, #This is that extra column we made earlier for which we will use for cell ID
             hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

#Save html plot for viewing
htmlwidgets::saveWidget(as_widget(p), paste0(filedir, "UMAP3D_NoHarmony_LayeredOnMeghansCluster_colorbyCluster_NoInit.html"))

#Save seurat object with umap reduction
saveRDS(integratedObject, paste0(filedir,'31+7_samples_umap_layered_excludeInitial_with3DUMAP.rds'))
