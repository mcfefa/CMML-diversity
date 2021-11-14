## Clonal tracing using mitochondrial variants

We use [mitoClone](https://github.com/veltenlab/mitoClone), a package previously published by [Velten et al. in Nature Communications (PMID: 33649320)](https://pubmed.ncbi.nlm.nih.gov/33649320/), to establish the clonal distribution of each of our samples. We focus on samples with sequential timepoints in order to quantify how the clonal construction shifts over time and in response to therapy. 

Given the bam and bai files from each single cell, we use `01_mitoclone_generateCountTables_11142021.R` to construct count tables, which denote the reads at each site in the mitochondrial genome. The count tables are provide in [our code ocean capsule](https://codeocean.com/capsule/1315403/tree), but this first step is shown for completeness.

Then, we take the count tables and use `02_runMitoclone_11142021.R` to reconstruct clonal relationships of the samples. 
