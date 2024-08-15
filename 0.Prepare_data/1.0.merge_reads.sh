#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/merge.%A_%a.out # Standard output
#SBATCH --array=1-20


## sanity check post load
my_job_header
####

#### user inputs
metadata="/gpfs2/scratch/jcnunez/Dsu.prelim.data/SampleId_metadata.txt"
seq_id=$( cat $metadata  | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

####### Locate directories
#######
dir_1="/netfiles/nunezlab/D_suzukii_resources/reads/KY_NTeets_2020_2021"
dir_2="/netfiles/nunezlab/D_suzukii_resources/reads/KY_round3"

##### Merge forward reads

F_a=$dir_1/${seq_id}_R1.*
echo $F_a
F_b=$dir_2/${seq_id}_*_R1_*
echo $F_b

cat $F_a $F_b > ${seq_id}.F.merged.fastq.gz

##### Merge reverse reads

R_a=$dir_1/${seq_id}_R2.*
echo $R_a
R_b=$dir_2/${seq_id}_*_R2_*
echo $R_b

cat $R_a $R_b > ${seq_id}.R.merged.fastq.gz


echo "done"
date