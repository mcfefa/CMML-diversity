## Using COMET to find markers for flow cytometry

We use [COMET](http://www.cometsc.com/index), a tool previously published by [Delaney et al. in Mol. Sys. Bio (PMID: 31657111)] (https://pubmed.ncbi.nlm.nih.gov/31657111/) to leverage scRNA-Seq to identify interesting populations in flow cytometry. We specifically use COMET to identify cluster 2 cells in our parallel flow cytometry dataset.

Given our seurat object, we use `01_makeComet_txtFiles.R` to make txt files matching the criteria for input to COMET. We distinguish our cells based on whether they appear in cluster 2, and downsample the non-cluster 2 cells to account for the maximum of the COMET package. 

We run COMET in the command line using the simple bash script `02_comet_commandLine.sh`
