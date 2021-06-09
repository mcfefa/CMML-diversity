# script to design primers to amplify selected mutations.
# Variants must be selected manually from the output of distance_polyA.sh
# the input file must have the following format: genesymbol_position with one gene per line
# -r refers to the read length of read2. It is crucial to specify the correct one.
# -p refers to the sample name

# specifications to run in UGE cluster (it can also be run locally)
#$ -cwd
#$ -e log_files/stderr.txt
#$ -o log_files/stdout.txt
#$ -N create_primers
#$ -pe smp 1
#$ -l h_rt=2:00:00
#$ -l virtual_free=40G
#$ -q short-sl7

#conda activate primers ;
Rscript scripts/design_mutation_primers.R -i annotated_variants.csv -o primers/ -r 120 -g data/Homo_sapiens.GRCh38.100.chr.gtf -b genomic_files/polyA_final_selection.bed -p P100 -l primers/targeted_mutations.txt
