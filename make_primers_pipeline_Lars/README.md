We start with a list of annotated variants which we want to focus on in MAF format (`MAF_All_GRCh37.maf`, derived from Eric's annotated csv `tru_variants_scRNA.csv`). This file contains info on the chromosome, symbol, position (as well as a bunch of other things) for all of the IWG samples. Note that for us, these are GRCh37, not GRCh38.

In addition to this file, we also need: bam and bai files for each sample (mapped to GRCh38), a seurat object for the sample, a chain.gz file for converting GRCh37 to GRCh38 `GRCh37_to_GRCh38.chain.gz`, and a gtf file for GRCh 38 (CHR).


1. The first step is to convert the MAF files to GRCh38. For this we use crossMap. The info for creating a conda environment with this tool is detailed in `crossMap_crossPlatform.yml`. The bash script for converting is given by `convert_maf_toGRCh38.sh`. Make sure it has the correct path to `GRCh37_to_GRCh38.chain.gz`

2. Then, we take the GRCh38 MAF file and use it as input for the shell script `pipeline_polyA_fileGeneration.sh`. This will call an R script, `polyA_fileGeneration.R` which:
  * makes a bed file for each sample containing the positions of the genes we look at in that sample. This will be used for subsetting the bam files. This uses biomaRt and the server can sometimes be inaccessible, so it might take a few runs to get bed files for all the samples.
  * Makes a csv with the required variant information, CHROM, POS, symbol. Note that POS is now the GRCh38 position.
  * Gives the path to the seurat object for the given sample. This won't be useful once the platform changes, but can be altered for the new workflow.
  
3. Now we are ready to subset the bam files. The script `subset_bamFiles.sh` will do this for all of the bam files which we have DNA sequencing for. Note: Change the directory for the updated workflow so it puts the subsetted bam file in the right spot. 

4. At this point, I moved off the cluster. I copied all folders from the cluster to my local computer. Everything should be ready for input into Lar's pipeline (see his README at `README_Lars`). We must build the "environment" to run this. First, install all of the R packages in `distance_polyA_tailored.R` and `design_mutation_primers_tailored.R`. Then, Download and install primer3_core (use Primer3 v.2.5.0 from github): (https://github.com/primer3-org/primer3/releases instructions: https://vcru.wisc.edu/simonlab/bioinformatics/programs/primer3/primer3_manual.htm#installLinux backup download: http://sourceforge.net/projects/primer3/). Move primer3_core into your PATH (/home/4472241/.conda/envs/exome_v2/bin/primer3_core on Moffitt Cluster, /usr/local/bin/primer3_core on my mac). Install makeblastdb and blastn tools using instructions here (use blastn v.2.6.0): (https://www.ncbi.nlm.nih.gov/books/NBK52640/ Source: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Move blastn and makeblastdb executables into your PATH variable as well. 
  * We have to make changes to Lar's pipeline to accomodate the notational differences between his samples and ours. Most notably, our bam files list chromosome only by number, instead of his which are "chr#". This will throw an error in his pipeline, but has been fixed in ours.
  * Run `distance_polyA_tailored.sh` which calls `distance_polyA_tailored.R`
  * Manually go through the `annotated_variants_CMML_#.csv` produced for each CMML sample to elect the mutations we want to make primers for. Make a .txt file for each sample in the format genesymbol_POS (KRAS_128798340), with a new mutation on each line. These are the mutations you elect to make primers for. Call this `targeted_mutations.txt` and place it in the `CMML_#` folder.
  * Run `design_primers_tailored.sh` which calls the R script `design_mutation_primers_tailored.R`
  * Load bed files and bam (with bai) in IGV to double check and make sure everything is properly placed. Be sure to change the genome to GRCh38
