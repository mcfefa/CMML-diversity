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
cd /Users/4472241/scCode/make_primers

for ((i =  1; i <= 31 ; i++)); do
    dir=/Users/4472241/scCode/make_primers/CMML_${i}
    if [[ -d "${dir}" ]]
    then
        name=CMML_${i}
        Rscript design_mutation_primers_tailored.R -i ${name}/annotated_variants_${name}.csv -o ${name}/ -r 120 -g data/Homo_sapiens.GRCh38.100.chr.gtf -b ${name}/polyA_final_selection.bed -p P100 -l ${name}/targeted_mutations.txt
    fi
done