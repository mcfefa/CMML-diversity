## Install devtools from CRAN
#install.packages("devtools")
## Use devtools to install hdf5r and loomR from GitHub
#devtools::install_github(repo = "hhoeflin/hdf5r")
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(loomR)

pt1file <- "/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/analysis/scDNAseq/Tapestri_Output_Files/5-M-001.cells.loom"

lfile <- connect(filename = pt1file, mode = "r+")
