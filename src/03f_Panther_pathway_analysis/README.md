## Differential Gene Expression Pathway Analysis

Previously, we found the mono-biased patients were enriched for cluster 2 and those patients had a statistically significant worse prognosis. To determine, how cluster 2 cells were different from the rest of the cells, we used the `FindMarkers` feature in `Seurat` to identify the genes that were differentially expressed in cluster 2 compared to all other cells. From that list, we split the genes into those that were significantly up-regulated in cluster 2 compared to others and those that were significantly down-regulated in cluster 2 compared to others. These two gene lists were then uploaded to [Enrichr](https://maayanlab.cloud/Enrichr/) a gene enrichment analysis tool to identify what pathways were up-regulated and down-regulated in cluster 2. 

More information about Enrichr is available from these publications: 
+ [Chen et al. BMC Bioinformatics, 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128) ([PMID: 23586463](https://pubmed.ncbi.nlm.nih.gov/23586463/))
+ [Kuleshov et al. Nucleic Acids Research, 2016](https://academic.oup.com/nar/article/44/W1/W90/2499357) ([PMID: 27141961](https://pubmed.ncbi.nlm.nih.gov/27141961/))
+ [Xie et al. Current Protocols, 2021](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.90) ([PMID: 33780170](https://pubmed.ncbi.nlm.nih.gov/33780170/))
