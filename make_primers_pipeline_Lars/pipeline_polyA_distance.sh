#!/bin/bash
#Full automated approach for calculating the polyA distance for eventual primer generation
#qsub -I -l nodes=1:ppn=12 -l walltime=4:00:00

cd /share/lab_altrock/MeghanCluster/BrianCluster/make_primers/

#Define input csv file with CMML_#, gene symbol, mutation POS and CHROM
maf_input='/share/lab_altrock/MeghanCluster/BrianCluster/make_primers/maf_all_GRCh38.csv'

for ((i = 1 ; i <= 21 ; i++)); do
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

    #samtools view -b -L ${bedfile} /share/lab_altrock/MeghanCluster/BrianCluster/BAM/${name}/${name}.bam > ./${name}/${name}.subset.bam
    #Create job to do the bam subsetting
    #echo "samtools view -b -L ${bedfile} /share/lab_altrock/MeghanCluster/BrianCluster/BAM/${name}/${name}.bam > ./${name}/${name}.subset.bam" \
    #| qsub -l nodes=1:ppn=4 -l walltime=1:00:00 -j oe -V
    #samtools index ${name}.subset.bam > ${name}.subset.bam.bai

done

#Now need to create a job and manually subset bam files (see github)
#Then the rest will be run locally.