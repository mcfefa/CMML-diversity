#Start screen to be able to go back and check work
screen 

#Request job, big memory, four hours is more than enough
qsub -I -q bigmemQ -l nodes=1:ppn=12 -l walltime=4:00:00

#cd to 4472241 home dir on cluster, activate env
cd ~
source ./Comet_env/bin/activate

#Navigate to the directory
cd /share/lab_altrock/MeghanCluster/BrianCluster/COMET/isolateClus2_newClustering

#Run comet, adjust the -K to set how many combinations (1=singlet, 2=pair+singlet, etc.)
Comet res=0.05_markers.txt res=0.05_umap.txt res=0.05_clusters.txt res=0.05_output/ -g res=0.05_genes.txt -C 2 -K 3
