#!/bin/bash

# Subset all the bam files and make bai index files for them (for 21 samples with mutation position)

#Request job (this one needs bigMem to run faster): qsub -I -q bigmemQ -l nodes=1:ppn=12 -l walltime=4:00:00

cd /share/lab_altrock/MeghanCluster/BrianCluster
module load samtools/1.9

for ((i = 1 ; i <= 31 ; i++)); do
    dir=/share/lab_altrock/MeghanCluster/BrianCluster/make_primers/CMML_${i}
    if [[ -d "${dir}" ]]
    then
        bedfile=./make_primers/CMML_${i}/CMML_${i}.bed
        name=CMML_${i}
        samtools view -b -L ${bedfile} /share/lab_altrock/MeghanCluster/BrianCluster/BAM/${name}/${name}.bam > \
        /share/lab_altrock/MeghanCluster/BrianCluster/make_primers/${name}/${name}.subset.bam
        samtools index /share/lab_altrock/MeghanCluster/BrianCluster/make_primers/${name}/${name}.subset.bam > \
        /share/lab_altrock/MeghanCluster/BrianCluster/make_primers/${name}/${name}.subset.bam.bai
    fi
done
