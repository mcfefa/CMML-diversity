#!/bin/bash
#First step automated approach for calculating the polyA distance for eventual primer generation
#Request job: qsub -I -l nodes=1:ppn=12 -l walltime=4:00:00

cd /share/lab_altrock/MeghanCluster/BrianCluster/make_primers/

#Define input csv file with CMML_#, gene symbol, mutation POS and CHROM
maf_input='/share/lab_altrock/MeghanCluster/BrianCluster/make_primers/maf_all_GRCh38.csv'

for ((i = 1 ; i <= 31 ; i++)); do
    rOutput=($(Rscript --vanilla polyA_fileGeneration.R ${maf_input} ${i}))
    name=${rOutput[0]}
    echo ${name} 
    bedfile=${rOutput[1]}
    echo ${bedfile} 
    csv_outfile=${rOutput[2]}
    echo ${csv_outfile} 
    pathToRDS=${rOutput[3]}
    echo ${pathToRDS}

    #Move rds file to appropriate directory
    scp ${pathToRDS} ./${name}/${name}.rds

done

#Now need to subset bam files (see github)
#Then the rest will be run locally.
