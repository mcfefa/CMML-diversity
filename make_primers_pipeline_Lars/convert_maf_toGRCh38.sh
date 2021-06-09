#Convert the MAF file for all samples to get correct position

#screen
#Request job: qsub -I -q bigmemQ -l nodes=1:ppn=12 -l walltime=104:00:00
#cd /share/lab_altrock/MeghanCluster/BrianCluster

conda activate crossMap

# Run Crossmap (see details: http://crossmap.sourceforge.net/#convert-vcf-format-files)
CrossMap.py  maf  GRCh37_to_GRCh38.chain.gz  ./make_primers/MAF_All_GRCh37.maf  genome.fa  GRCh38 ./make_primers/MAF_ALL_GRCh38.maf
