#!/usr/bin/env bash
#
#SBATCH -J clean # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/clean.%A_%a.out # Standard output
#SBATCH --array=1-20


## sanity check post load
my_job_header
####
NTHREADS=$SLURM_CPUS_ON_NODE
echo $NTHREADS

#### user inputs
metadata="/gpfs2/scratch/jcnunez/Dsu.prelim.data/SampleId_metadata.txt"
seq_id=$( cat $metadata  | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

####Locations
dir_data=/gpfs2/scratch/jcnunez/Dsu.prelim.data/Reads_Merged
json_folder=/gpfs2/scratch/jcnunez/Dsu.prelim.data/res_clark_jsons
out_dir=/gpfs2/scratch/jcnunez/Dsu.prelim.data/Clean_Reads

mkdir -p $out_dir

### Softwares
fastp=/netfiles/nunezlab/Shared_Resources/Software/fastp_latest/fastp

#cleanning seqeunce
$fastp \
-i $dir_data/${seq_id}.F.merged.fastq.gz \
-I $dir_data/${seq_id}.R.merged.fastq.gz \
-o $out_dir/${seq_id}.clean.F.gz \
-O $out_dir/${seq_id}.clean.R.gz \
-w $NTHREADS \
-h $json_folder/${QUERY_NAME}.html \
-j $json_folder/${QUERY_NAME}.json

#rm ${QUERY_NAME}.json
#rm $dir_data/${QUERY_NAME}.fq_1.gz $dir_data/${QUERY_NAME}.fq_2.gz




