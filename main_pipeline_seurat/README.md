The first file, "QC+Plots_normal+CMML.R" contains the quality controlling for both the normal samples and CMML samples. Input required is the unprocessed data as an rds file (or 10X directory for the case of Hua 3)
The second file, "seurat_standard_pipeline+harmony.R", takes the output of the first file and runs many of the standard operations (scaling, dimension reduction etc.)
