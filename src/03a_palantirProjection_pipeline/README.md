## Differentiation Trajectory Analysis

We leverage the Differentiation Trajectory Analysis on CD34+ cells enriched from normal bone marrow mononuclear cells previously published by [Setty et al in Nature Biotechnology in 2019](https://www.nature.com/articles/s41587-019-0068-4) ([PMID: 30899105](https://pubmed.ncbi.nlm.nih.gov/30899105/)).  We projected each patient sample in our cohort onto the Rep 1 reference CD34+ normal single-cell RNA sequencing data to evaluate in which lineages CMML patients have over-density and under-density of cells. 

`01_map_to_rep1_use123.R` is the code for mapping each individual patients onto the normal/reference (Rep1) from Setty et al.
`02_densityPlot_11012021.R` is the code for visualizing the over-densities and under-densities along the differentiation trajectories.
