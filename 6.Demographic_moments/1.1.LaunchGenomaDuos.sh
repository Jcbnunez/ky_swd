#!/bin/sh
#
#SBATCH -J genoDuo
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=3:00:00
#SBATCH -o ./slurmOut/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-68

module load Rgeospatial

guide=VA_KY_pairs.txt

Rscript --vanilla 1.Genomalicious.SWD.R \
${SLURM_ARRAY_TASK_ID} \
$guide
