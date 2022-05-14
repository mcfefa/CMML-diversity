# Install comet / activate env where comet is installed 

#Navigate to the directory where we have input files
cd /Users/Brian/scCode/Comet_04252022

#Run comet, adjust the -K to set how many combinations (1=singlet, 2=pair+singlet, etc.)
Comet res=0.05_markers.txt res=0.05_umap.txt res=0.05_clusters.txt res=0.05_output/ -g res=0.05_genes.txt -C 2 -K 3
