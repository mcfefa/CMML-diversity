## Cell Type Assignment

Here, we used [SingleR](https://github.com/dviraran/SingleR) for cell type recognition of single-cell RNA sequencing data, leveraging reference trancriptomic datasets of pure cell types to infer the cell of origin. More information about the methodology and implementation of SingleR is available from the manuscript [Aran et al. Nature Immunology, 2019](https://www.nature.com/articles/s41590-018-0276-y) ([PMID: 30643263](https://pubmed.ncbi.nlm.nih.gov/30643263/)). 

`01_singleR_threeRefs.R` generates cell type estimates using three built-in references from singleR.

`02_singleR_Rapin.R` uses a reference from publicly available data to esimate cell type using singleR.

`03_consensus_singleR.R` takes the four independent cell type assignments and assigns a consensus cell type. 
