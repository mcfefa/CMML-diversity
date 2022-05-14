## Gene Expression Signature Analysis 

This analysis contains two distinct parts. 

(1) We compare previously gene signatures across different hematopoietic stem and progenitor populations on our pseudo-bulk aggregation of single-cell data to produce a heatmap with distinct samples on one axis and genes from each signature on the other. For this, we leverage the gene signatures published by [Wu et al in Blood Advances in 2020](https://ashpublications.org/bloodadvances/article/4/12/2702/460973/Sequencing-of-RNA-in-single-cells-reveals-a) ([PMID: 32556286](https://pubmed.ncbi.nlm.nih.gov/32556286/)) and focused on the HSC, GMP, and MEP populations based on the over-densities and under-densities we found in our Differentiation Trajectory Analysis (analysis [03a](https://github.com/mcfefa/CMML-diversity/tree/main/src/03a_palantirProjection_pipeline)). The code for this analysis is described in `Wu_pseudobulk_heatmap`.

(2) We used hematopoietic stem cell (HSC) gene signatures to produce single-cell resolution heatmaps comparing normal HSC expression of HSC markers from three independent HSC signatures to construct patient-specific heatmaps with genes as one axis and individual cells as the other. The HSC gene signatures used were from: 
+ Wu: [Wu et al.  Blood Advances, 2020](https://ashpublications.org/bloodadvances/article/4/12/2702/460973/Sequencing-of-RNA-in-single-cells-reveals-a) ([PMID: 32556286](https://pubmed.ncbi.nlm.nih.gov/32556286/))
+ VanGalen: [Van Galen et al. Cell Reports, 2018](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31597-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124718315973%3Fshowall%3Dtrue) ([PMID: 30380403](https://pubmed.ncbi.nlm.nih.gov/30380403/))
+ Eppert: [Eppert et al. Nature Medicine, 2011](https://www.nature.com/articles/nm.2415) ([PMID: 21873988](https://pubmed.ncbi.nlm.nih.gov/21873988/))

The code for this analysis is described in `HSC_depletion_scRes_eppert_wu_vg`.

(3) `WNT_signature` details the signature scoring and plots shown in Supp. Figure S10.  Signature derives from:
"Herault A, Binnewies M, Leong S, Calero-Nieto FJ, Zhang SY, Kang YA, et al. Myeloid progenitor cluster formation drives emergency and leukaemic myelopoiesis. Nature 2017;544(7648):53-8 doi 10.1038/nature21693."



