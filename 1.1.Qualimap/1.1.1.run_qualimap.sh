#!/usr/bin/env bash
#
#SBATCH -J quali # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 #<= this may depend on your resources
#SBATCH --mem 60G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/quali.%A_%a.out # Standard output
#SBATCH -p bluemoon
#SBATCH --array=1

##-20

JAVAMEM=59G
CPU=5
SAMP=${SLURM_ARRAY_TASK_ID}
QUAL=60

qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap
root=/gpfs2/scratch/jcnunez/Dsu.prelim.data/bam_files_final

####
mkdir all_KY_quals
cd all_KY_quals

mkdir KY${SAMP}_qual

$qualimap bamqc \
-bam $root/KY${SAMP}.dedup.bam  \
-outdir  KY${SAMP}_qual \
-nt $CPU \
--java-mem-size=$JAVAMEM

echo "done"
