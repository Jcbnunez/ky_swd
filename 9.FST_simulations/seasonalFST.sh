#!/bin/sh

# Specify a partition
#SBATCH --account=pi-jcnunez
#SBATCH --partition=general
# Request nodes
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Request memory
#SBATCH --mem=2G
# Run for some time
#SBATCH --time=1:00:00 
# Name job:
#SBATCH --job-name=AF_10-80
# Name output file
#SBATCH --output=%x_%j.out
#SBATCH --array=1-9

#change to directory in which to submit code
cd /users/a/r/armccrac/suzukii

## load modules if any
module load slim/4.3

## print job header
my_job_header

####
# Working directory
outdir="/users/a/r/armccrac/suzukii/8.27.25_AF_0.1-0.8_out"
wd="/users/a/r/armccrac/suzukii"

cd $outdir

# Parameter file
paramFile=${wd}/par.txt

# Extract constants from parameter file in epoch1
BaseName=$(cat $paramFile | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d"  )
wt=$(cat $paramFile | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d"  )



## get filepath to output files for epoch2_par


## sanity check
echo 'winter temp: ' $wt
echo 'BaseName: ' $BaseName

#### Run SLiM

for replicate in {1..10}; do

    # Progress
    echo 'Replicate:' $replicate
    slim \
    -d BaseName="'${BaseName}'" \
    -d replicatei=$replicate \
    -d wt=$wt \
    ${wd}/8.12.25_AF_FST.slim
done 





    
    