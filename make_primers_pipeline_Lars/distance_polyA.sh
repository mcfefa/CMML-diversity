# script to compute the distance from the mutation to the polyA of the gene.
# It requires:
# - csv file with variant information
# - subsetted BAM file for the genes of interest (it can also be the complete file but it will take longer if it is very large)
# - gtf file
# - seurat object with raw counts (ideally from the same sample)

# specifications to run in UGE cluster (it can also be run locally)
#$ -cwd
#$ -e log_files/stderr.txt
#$ -o log_files/stdout.txt
#$ -N create_primers
#$ -pe smp 8
#$ -l h_rt=2:00:00
#$ -l virtual_free=40G
#$ -q short-sl7

# conda env create -f scripts/conda_primers.yml -n primers
#conda activate primers ;
Rscript scripts/distance_polyA.R -i data/variants.csv -b data/subset.bam -g data/Homo_sapiens.GRCh38.100.chr.gtf -c data/seurat_object.rds -d genomic_files/ -o annotated_variants.csv -e FALSE
