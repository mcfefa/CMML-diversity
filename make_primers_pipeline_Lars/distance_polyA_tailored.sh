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
cd /Users/4472241/scCode/make_primers/
for ((i = 1 ; i <= 31 ; i++)); do
    dir=/Users/4472241/scCode/make_primers/CMML_${i}
    if [[ -d "${dir}" ]]
    then
        Rscript ./distance_polyA_tailored.R -i ./CMML_${i}/CMML_${i}_variants.csv -b ./CMML_${i}/CMML_${i}.subset.bam -g data/Homo_sapiens.GRCh38.100.chr.gtf -c ./CMML_${i}/CMML_${i}.rds -d CMML_${i}/ -o annotated_variants_CMML_${i}.csv -e FALSE
    fi
done