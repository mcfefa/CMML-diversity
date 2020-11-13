#Make a nice 3D UMAP
library(Seurat)
library(plotly)

#Clear workspace variables
rm(list=ls())

#Set file directory to save to
filedir <- "/Users/4472241/scCode/layerUMAP+ClusteringNoInit/"

#Import the object you want to plot
integratedObject <- readRDS('/Users/4472241/scCode/layerUMAP+ClusteringNoInit/31+7_samples_umap_layered_excludeInitial.rds')

# Re-run UMAPs if there's no 3D UMAP reduction in the object
integratedObject <- RunUMAP(integratedObject, dims = 1:50, n.components = 3L, reduction = "pca",
                reduction.key = "UMAP3D_", reduction.name = "umap3D")

## Grouping Cells as CMML, HMA, or Normal in metadata to plot by condition
integratedObject$sampletype <- plyr::mapvalues(
        x = integratedObject$orig.ident, 
        from = c("SettyPt1", "SettyPt2", "SettyPt3", "CD34","HuaPt1","HuaPt2","HuaPt4","LTB3966","LTB4121","LTB6169","4J003","4K001","4Q001","5E001","5H001","SF100109106293","SF100109111451","SF100109110236","SF14040100158","SF14060200025","SF12062800475","SF14072200012","SF13061200056","4S001","2V001","SF14101000049","SF16112900158","6AE001","6AC001","6AD001","SF100109101914","SF12042500035","SF12092600014","SF14031800065","SF14050700419","SF16026800045","SF16072200003","SF16112300029",'SF15010200008', 'SF14111400033','SF14061300036','SF13032800016', 'SF13070900171', 'SF14092500135', 'SF14110400108', 'SF14080400065'), 
        to = c("Normal", "Normal", "Normal", "Normal","Normal","Normal","Normal","CMML","CMML","CMML","Normal","HMA","HMA","HMA","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","HMA","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML","CMML", 'RUX', 'RUX', 'Chemo', 'HMA', 'HMA', 'HMA', 'HMA', 'HMA')
)
#Ignore warning if using NoInit object, as some samples are supposed to be missing

# Extract UMAPinformation from Seurat Object
umap_1 <- integratedObject[["umap3D"]]@cell.embeddings[,1]
umap_2 <- integratedObject[["umap3D"]]@cell.embeddings[,2]
umap_3 <- integratedObject[["umap3D"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = integratedObject, reduction = "umap3D")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = integratedObject, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3", "sampletype"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))



# Plot your data, in this example my Seurat object had 21 clusters (0-20)
p <- plot_ly(data = plot.data, 
        x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, 
        color = ~sampletype, 
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
htmlwidgets::saveWidget(as_widget(p), paste0(filedir, "UMAP3D_NoHarmony_LayeredOnMeghansCluster_colorbyCondition_NoInit.html"))

#Save seurat object with umap reduction
saveRDS(integratedObject, paste0(filedir,'31+7_samples_umap_layered_excludeInitial_with3DUMAP.rds'))
