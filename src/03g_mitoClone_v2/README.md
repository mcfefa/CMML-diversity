"01_mitoClone_generateCountTables.R" takes bam and bai files for each single cell as input and constructs count tables, which contain the base counts at each site for each cell. This will be used as input to mitoClone's next steps. We provide count tables as an rds file in CodeOcean, so this step is only necessary if running on a different dataset.

"02_mitoClone_mutationCalls_and_cluster.R" runs the mitoclone functions and saves the output mitoclone objects. Mitoclone functions include filtering sites to identify which sites will be used for clonal assignment and clustering to assign clones. This script will call mutationCallsFromBlacklist_singleCoreFn.R because I needed to manually edit one of the mitoClone functions to run only on a single core (multicore function was not working).

"umap_embeddings_pre_mito_02022022.rds" aggregates all information from seurat object, including palantir embeddings, umap embeddings, cluster assignment, and singleR cell type assignment.

"mitoclone_in_palantir_plotAll.R" takes all mitoclone objects computed in "02_mitoClone_mutationCalls_and_cluster.R" and plots them with position of palantir embeddings and color of mitoClone assigned clone. This is shown in figure 4A of the paper for three representative examples. 
