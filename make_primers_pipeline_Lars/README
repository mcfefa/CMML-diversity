This folder contains all necessary data and scripts to generate primers to amplify mutations from
3' scRNAseq from 10x. The pipeline consists in computing the distance from the mutation to the polyA as well as
the raw counts per cell and then design primers to amplify the manually selected mutations of interest.

I included a yml file which contains the info to create a conda environment with all the available packages.
The pipeline only requires R packages specified in the Rscripts so it can also be run outside of conda.

-------------------------------------------------------------------------------------------------------------

1.- Compute distance to polyA

The 1st script distance_polyA.sh computes the distance of the mutation to the end of the coding region.
This is key in deciding which mutations to target. Mutations which are further away from 1.5 kb to the polyA
will be very difficult to amplify in Illumina sequencers (they have a strong preference towards shorter
fragments, particularly HiSeq and MiSeq). The script requires the following files:

  - variants.csv: csv file with the variants of interest found in the sample. It must contain the columns
                  symbol, CHROM, POS indicating the gene symbol, chromosome and position of the variant.

  - subset.bam: bam file from the sample subsetted for the genes with variants. The complete BAM file
                from cellranger can be included but it will take longer to run. I generally create a bed file
                with the coordinates of the genes in variants.csv and subset the BAM file with samtools.

  - gtf_file: gtf file with information about exons, genes... I included one for the GRCh38 in the data folder.

  - seurat_object.rds: seurat object with raw counts for the sample of interest. This is used to estimate the
                      raw average counts per cell for each of the genes with mutations. This is key when deciding
                      what mutations to amplify as the lower the expression the less cells will be covered.

The output file annotated_variants.csv contains two additional columns distance_3_end and counts_cells for
each of the mutated sites in variants.csv. Mutations for which primers will be designed have to be selected
manually. As mentioned above, we generally only consider mutations which are <1.5-2kb from the polyA of the gene.
Also it is important to make sure that the gene is expressed in our sample.

It is always good practice to check the inferred polyAs in a genome browser. The pipeline generates bed files
with the polyAs and mutation coordinates which can be loaded into IGV for example. I attached a screenshot
of one example (distance_polyA.png).

Before executing the design_mutation_primers.sh script a txt file (targeted_mutations.txt) with the selected mutations has to be
saved in the primers folder. The format must be symbol_position (e.g. NPM1_171410539), one mutation per line.

-------------------------------------------------------------------------------------------------------------

2.- Design primers for the mutations of interest

An outer, middle and 4 staggered inner primers will be designed for each mutation. In this script
is key to specify the read length of read2. We recommend to use at least 75bp so that there is enough sequence
for the inner primers to anneal.

The script generates two csv files, primers_details.csv with primer information provided by Primer3 and
primers_staggered.csv with the corresponding sequences. Finally a bed file with the genomic coordinates
of the primers is also exported.

I usually also check that the order of the primer makes sense and that they fall into exons by loading
them in IGV. You can see an example of this in 'browser_images/primers_igv.png'.
