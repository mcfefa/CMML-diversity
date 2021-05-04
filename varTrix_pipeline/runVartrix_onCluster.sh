#screen
#Request job: qsub -I -q bigmemQ -l nodes=1:ppn=12 -l walltime=104:00:00
#cd /share/lab_altrock/MeghanCluster/BrianCluster

conda activate liftOver

CrossMap.py  vcf  GRCh37_to_GRCh38.chain.gz  ./varTrix_CMML_20/145249.IWG-219_S36.smCounter.anno.vcf  genome.fa  ./varTrix_CMML_20/out.hg38.145249.IWG-219_S36.smCounter.anno.vcf

#Then manually convert "chr#" to "#" in hg38 vcf file (if needed, depends on vcf file)

samtools faidx genome.fa
./vartrix_linux -v ./varTrix_CMML_20/out.hg38.145249.IWG-219_S36.smCounter.anno.vcf -b ./varTrix_CMML_20/CMML_20.bam -f genome.fa -c ./varTrix_CMML_20/CMML_20.tsv -o ./varTrix_CMML_20/outputMatrix.mtx

module load bcftools/1.9
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%ID][\t%FILTER]\n' "./varTrix_CMML_20/out.hg38.145249.IWG-219_S36.smCounter.anno.vcf" > "./varTrix_CMML_20/SNV.loci.withID.txt"
sed -i 's/\s/:/g' ./varTrix_CMML_20/SNV.loci.withID.txt

